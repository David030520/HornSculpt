#include "bem_hf_test_oracle.hpp"
#include <algorithm>
#include <array>
#include <cuda/galerkin_pair_moments.cuh>
#include <cuda/regular_pair_contraction.cuh>
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace {
using namespace hornsculpt;
using C = std::complex<double>;
using Z = cuDoubleComplex;
using Rule = std::vector<galerkin::TriangleQuadraturePoint>;
constexpr C imaginary{0, 1};
constexpr double budget = 1.25e-6;
void require(bool condition, const char *message) {
    if (!condition)
        throw std::runtime_error(message);
}
void check(cudaError_t status) {
    if (status != cudaSuccess)
        throw std::runtime_error(cudaGetErrorString(status));
}
template <class T> class Device {
  public:
    T *data{};
    explicit Device(std::size_t count) {
        check(cudaMalloc(reinterpret_cast<void **>(&data), count * sizeof(T)));
    }
    ~Device() {
        cudaFree(data);
    }
    Device(const Device &) = delete;
    Device &operator=(const Device &) = delete;
};
__global__ void actualRegular(const cuda_detail::Panel *panels,
                              const galerkin::TriangleQuadraturePoint *rule, unsigned count,
                              double k, Z *result, int *error) {
    auto moments = cuda_detail::integrateRegularPairThread<true>(panels[0], panels[1], rule, count,
                                                                 k, error, threadIdx.x, 32);
    cuda_detail::reduceWarp(moments);
    if (threadIdx.x == 0) {
        for (unsigned i = 0; i < 3; ++i) {
            result[9 + i] = make_cuDoubleComplex(moments.rhs[i].real, moments.rhs[i].imag);
            for (unsigned j = 0; j < 3; ++j) {
                const auto value =
                    cuda_detail::pressureCoefficient(moments, panels[0], panels[1], i, j, k);
                result[3 * i + j] = make_cuDoubleComplex(value.x, value.y);
            }
        }
    }
}
cuda_detail::Panel packed(const galerkin::Panel &panel) {
    cuda_detail::Panel result{};
    for (unsigned i = 0; i < 3; ++i) {
        result.vertices[i] = panel.vertices[i];
        result.nodes[i] = panel.nodes[i];
        result.curls[i] = panel.curls[i];
    }
    result.centroid = panel.centroid;
    result.normal = panel.normal;
    result.area = panel.area;
    result.edge = panel.maximumEdgeLength;
    return result;
}

__global__ void actualContracted(const cuda_detail::Panel *panels,
                                 const galerkin::TriangleQuadraturePoint *rule, unsigned count,
                                 double k, Z *result, int *error) {
    __shared__ cuda_detail::RegularPairScratch scratch;
    const auto value =
        cuda_detail::contractRegularPairWarp(panels[0], panels[1], rule, count, k, scratch, error);
    if (threadIdx.x < 12)
        result[threadIdx.x] = value;
}
Rule independentRule(unsigned count) {
    Rule result;
    const auto gauss = test::wideband_oracle::gauss(count);
    for (const auto &[u, wu] : gauss)
        for (const auto &[v, wv] : gauss)
            result.push_back({{1 - u, u * (1 - v), u * v}, wu * wv * u});
    return result;
}
std::array<C, 12> originalWeak(const galerkin::Panel &xPanel, const galerkin::Panel &yPanel,
                               double k, const Rule &rule) {
    using namespace test::wideband_oracle;
    std::array<C, 12> result{};
    for (const auto &x : rule)
        for (const auto &y : rule) {
            const auto delta = difference(xPanel, x.barycentric, yPanel, y.barycentric);
            const double r = std::sqrt(dot(delta, delta));
            const double weight = 4 * xPanel.area * yPanel.area * x.weight * y.weight;
            const C green = std::exp(-imaginary * k * r) / (4 * std::numbers::pi * r);
            const C radial = green * (1. + imaginary * k * r) / (r * r);
            const C sourceDerivative = radial * dot(delta, yPanel.normal);
            const C targetDerivative = -radial * dot(delta, xPanel.normal);
            for (unsigned i = 0; i < 3; ++i) {
                result[9 + i] +=
                    weight * x.barycentric[i] * (-green + imaginary * targetDerivative / k);
                for (unsigned j = 0; j < 3; ++j) {
                    const double basis = x.barycentric[i] * y.barycentric[j];
                    const C h = green * (dot(xPanel.curls[i], yPanel.curls[j]) -
                                         k * k * basis * dot(xPanel.normal, yPanel.normal));
                    result[3 * i + j] += weight * (-basis * sourceDerivative - imaginary * h / k);
                }
            }
        }
    return result;
}
double error(const std::array<C, 12> &actual, const std::array<C, 12> &expected, unsigned channel) {
    const unsigned first = channel ? 9 : 0, last = channel ? 12 : 9;
    double difference{}, norm{};
    for (unsigned i = first; i < last; ++i) {
        difference += std::norm(actual[i] - expected[i]);
        norm += std::norm(expected[i]);
    }
    require(norm > 0, "Regular arithmetic oracle is degenerate.");
    return std::sqrt(difference / norm);
}
void diagnose() {
    constexpr double k = 2 * std::numbers::pi * 20000 / 343;
    const auto referenceRule = independentRule(16), comparisonRule = independentRule(12);
    Device<cuda_detail::Panel> panels(2);
    Device<Z> output(12);
    Device<int> status(1);
    double maximumArithmeticError{};
    for (double kr : {.25, .5, 1., 4., 8., 16., 32., 400., 4096., 1048576.})
        for (double shift : {0., 1e6}) {
            SurfaceMesh mesh;
            const double distance = kr / k;
            const std::array<Vec3, 4> tetra{
                {{0, 0, 0}, {.000021, 0, 0}, {0, .000017, 0}, {0, 0, .000013}}};
            for (unsigned body = 0; body < 2; ++body) {
                const std::uint32_t first = static_cast<std::uint32_t>(mesh.vertices.size());
                for (auto vertex : tetra)
                    mesh.vertices.push_back({vertex.x + shift + body * .6 * distance,
                                             vertex.y + shift,
                                             vertex.z + shift + body * .8 * distance});
                for (const std::array<std::uint32_t, 3> ids :
                     {std::array<std::uint32_t, 3>{0, 2, 1}, {0, 1, 3}, {0, 3, 2}, {1, 2, 3}})
                    mesh.triangles.push_back(
                        {{first + ids[0], first + ids[1], first + ids[2]}, body + 1});
            }
            const auto space = galerkin::buildSpace(mesh);
            const auto &x = space.panels[7], &y = space.panels[1];
            const auto expected = originalWeak(x, y, k, referenceRule);
            const auto refined = originalWeak(x, y, k, comparisonRule);
            require(std::max(error(refined, expected, 0), error(refined, expected, 1)) < 1e-9,
                    "Regular arithmetic oracle did not converge independently.");
            const std::array<cuda_detail::Panel, 2> host{packed(x), packed(y)};
            check(cudaMemcpy(panels.data, host.data(), sizeof(host), cudaMemcpyHostToDevice));
            for (unsigned degree : {4u, 6u, 8u}) {
                const auto rule = galerkin::triangleRule(degree);
                Device<galerkin::TriangleQuadraturePoint> deviceRule(rule.size());
                check(cudaMemcpy(deviceRule.data, rule.data(), rule.size() * sizeof(rule[0]),
                                 cudaMemcpyHostToDevice));
                for (bool contracted : {false, true}) {
                    check(cudaMemset(status.data, 0, sizeof(int)));
                    if (contracted)
                        actualContracted<<<1, 32>>>(panels.data, deviceRule.data,
                                                    static_cast<unsigned>(rule.size()), k,
                                                    output.data, status.data);
                    else
                        actualRegular<<<1, 32>>>(panels.data, deviceRule.data,
                                                 static_cast<unsigned>(rule.size()), k, output.data,
                                                 status.data);
                    check(cudaGetLastError());
                    std::array<Z, 12> deviceResult;
                    int flag{};
                    check(cudaMemcpy(deviceResult.data(), output.data, sizeof(deviceResult),
                                     cudaMemcpyDeviceToHost));
                    check(cudaMemcpy(&flag, status.data, sizeof(flag), cudaMemcpyDeviceToHost));
                    require(flag == 0, "Actual regular kernel rejected finite diagnostic data.");
                    std::array<C, 12> actual;
                    for (unsigned i = 0; i < 12; ++i)
                        actual[i] = {deviceResult[i].x, deviceResult[i].y};
                    const auto sameRule = originalWeak(x, y, k, rule);
                    const double qa = error(sameRule, expected, 0),
                                 qb = error(sameRule, expected, 1);
                    const double aa = error(actual, sameRule, 0), ab = error(actual, sameRule, 1);
                    std::cout << std::setprecision(12)
                              << (contracted ? "Tensor CUDA arithmetic kh="
                                             : "Regular CUDA arithmetic kh=")
                              << k * std::max(x.maximumEdgeLength, y.maximumEdgeLength)
                              << ", kD=" << kr << ", shift=" << shift << ", degree=" << degree
                              << ", FP64RuleA/B=" << qa << '/' << qb
                              << ", actualMixedArithmeticA/B=" << aa << '/' << ab << '\n';
                    require(qa < budget && qb < budget,
                            "Small-kh arithmetic fixture is not quadrature-resolved.");
                    require(
                        aa < budget && ab < budget,
                        "Actual regular arithmetic exceeds the independent FP64 weakform budget.");
                    maximumArithmeticError = std::max({maximumArithmeticError, aa, ab});
                }
            }
        }
    std::cout << "Actual regular arithmetic maximum A/B error=" << maximumArithmeticError << '\n';
}
} // namespace
int main() {
    int count{};
    if (cudaGetDeviceCount(&count) != cudaSuccess || !count)
        return 77;
    try {
        diagnose();
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
