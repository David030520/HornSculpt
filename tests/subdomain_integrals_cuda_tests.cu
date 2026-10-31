#include "bem_hf_test_oracle.hpp"
#include <algorithm>
#include <array>
#include <cuda/subdomain_integrals.cuh>
#include <galerkin/regular_qualification.hpp>
#include <iostream>
#include <stdexcept>

namespace {
using namespace hornsculpt;
using C = std::complex<double>;
using Rule = std::vector<galerkin::TriangleQuadraturePoint>;
using Moments = std::array<C, 22>; // pressure[9], interface flux[9], DP0 rhs[3], G
void require(bool value, const char *message) {
    if (!value)
        throw std::runtime_error(message);
}
void check(cudaError_t status) {
    if (status != cudaSuccess)
        throw std::runtime_error(cudaGetErrorString(status));
}
template <class T> struct Device {
    T *data{};
    explicit Device(std::size_t n) {
        check(cudaMalloc(reinterpret_cast<void **>(&data), n * sizeof(T)));
    }
    ~Device() {
        cudaFree(data);
    }
    Device(const Device &) = delete;
    Device &operator=(const Device &) = delete;
};
cuda_detail::Panel pack(const galerkin::Panel &p) {
    cuda_detail::Panel r{};
    for (unsigned i = 0; i < 3; ++i) {
        r.vertices[i] = p.vertices[i];
        r.nodes[i] = p.nodes[i];
        r.curls[i] = p.curls[i];
    }
    r.centroid = p.centroid;
    r.normal = p.normal;
    r.area = p.area;
    r.edge = p.maximumEdgeLength;
    return r;
}
__global__ void integrate(const cuda_detail::Panel *panel,
                          const galerkin::PairQuadraturePoint *rule, unsigned count, double k,
                          bool mixed, cuda_detail::SubdomainMoments *out, int *error) {
    cuda_detail::SubdomainMoments interface{}, physical{};
    for (unsigned i = 0; i < count; ++i) {
        cuda_detail::accumulateSubdomainMoment(interface, panel[0], panel[1], rule[i], k, false,
                                               true, mixed, error);
        cuda_detail::accumulateSubdomainMoment(physical, panel[0], panel[1], rule[i], k, false,
                                               false, mixed, error);
    }
    for (unsigned i = 0; i < 3; ++i)
        interface.rhs[i] = physical.rhs[i];
    *out = interface;
}
Rule duffy(unsigned n) {
    Rule r;
    for (const auto &[u, wu] : test::wideband_oracle::gauss(n))
        for (const auto &[v, wv] : test::wideband_oracle::gauss(n))
            r.push_back({{1 - u, u * (1 - v), u * v}, wu * wv * u});
    return r;
}
std::vector<galerkin::PairQuadraturePoint> pairs(const Rule &rule) {
    std::vector<galerkin::PairQuadraturePoint> out;
    for (const auto &x : rule)
        for (const auto &y : rule) {
            galerkin::PairQuadraturePoint p{};
            for (unsigned i = 0; i < 3; ++i) {
                p.test[i] = x.barycentric[i];
                p.trial[i] = y.barycentric[i];
            }
            p.weight = x.weight * y.weight;
            out.push_back(p);
        }
    return out;
}
Moments reference(const galerkin::Panel &x, const galerkin::Panel &y, double k, const Rule &rule) {
    using namespace test::wideband_oracle;
    constexpr C I{0, 1};
    Moments out{};
    for (const auto &a : rule)
        for (const auto &b : rule) {
            const Vec3 d = difference(x, a.barycentric, y, b.barycentric);
            const double r = std::sqrt(dot(d, d)), w = 4 * x.area * y.area * a.weight * b.weight;
            const C g = w * std::exp(-I * k * r) / (4 * std::numbers::pi * r);
            const C radial = g * (1. + I * k * r) / (r * r);
            const C A = radial * dot(d, y.normal) + I * k * dot(x.normal, y.normal) * g;
            const C B = g + I * (-radial * dot(d, x.normal)) / k;
            out[21] += g;
            for (unsigned i = 0; i < 3; ++i) {
                out[18 + i] += B * a.barycentric[i];
                for (unsigned j = 0; j < 3; ++j) {
                    const double basis = a.barycentric[i] * b.barycentric[j];
                    out[3 * i + j] += A * basis;
                    out[9 + 3 * i + j] += B * basis;
                }
            }
        }
    return out;
}
Moments gpu(const galerkin::Panel &x, const galerkin::Panel &y, double k, const Rule &rule,
            bool mixed) {
    const std::array<cuda_detail::Panel, 2> h{pack(x), pack(y)};
    const auto q = pairs(rule);
    Device<cuda_detail::Panel> p(2);
    Device<galerkin::PairQuadraturePoint> d(q.size());
    Device<cuda_detail::SubdomainMoments> o(1);
    Device<int> status(1);
    check(cudaMemcpy(p.data, h.data(), sizeof(h), cudaMemcpyHostToDevice));
    check(cudaMemcpy(d.data, q.data(), q.size() * sizeof(q[0]), cudaMemcpyHostToDevice));
    check(cudaMemset(status.data, 0, sizeof(int)));
    integrate<<<1, 1>>>(p.data, d.data, static_cast<unsigned>(q.size()), k, mixed, o.data,
                        status.data);
    check(cudaGetLastError());
    cuda_detail::SubdomainMoments z{};
    int error{};
    check(cudaMemcpy(&z, o.data, sizeof(z), cudaMemcpyDeviceToHost));
    check(cudaMemcpy(&error, status.data, sizeof(error), cudaMemcpyDeviceToHost));
    require(error == 0, "GPU rejected finite pair");
    Moments r{};
    for (unsigned i = 0; i < 9; ++i) {
        r[i] = {z.pressure[i].real, z.pressure[i].imag};
        r[9 + i] = {z.flux[i].real, z.flux[i].imag};
    }
    for (unsigned i = 0; i < 3; ++i)
        r[18 + i] = {z.rhs[i].real, z.rhs[i].imag};
    r[21] = {z.green.real, z.green.imag};
    return r;
}
double groupError(const Moments &a, const Moments &b, unsigned first, unsigned last) {
    double delta = 0, norm = 0;
    for (unsigned i = first; i < last; ++i) {
        delta += std::norm(a[i] - b[i]);
        norm += std::norm(b[i]);
    }
    require(std::isfinite(delta) && std::isfinite(norm) && norm > 0,
            "Non-finite or zero oracle moment");
    return std::sqrt(delta / norm);
}
void compare(const Moments &a, const Moments &b, double tolerance, const char *message) {
    for (auto range :
         {std::pair{0u, 9u}, std::pair{9u, 18u}, std::pair{18u, 21u}, std::pair{21u, 22u}})
        require(groupError(a, b, range.first, range.second) < tolerance, message);
}
void run() {
    for (double kr : {.25, 16., 1048576.})
        for (double scale : {1e-8, 1., 1e8, 1e16})
            for (unsigned degree : {4u, 6u, 8u}) {
                const double k = kr / scale;
                const double targetKh = degree == 4 ? .1 : degree == 6 ? .6 : 1.5;
                const double h = kr == .25 ? .01 * scale : targetKh / (1.6 * k);
                const double shift = scale * 1e3;
                SurfaceMesh mesh;
                for (unsigned body = 0; body < 2; ++body) {
                    const double x = shift + body * .6 * scale, z = shift + body * .8 * scale;
                    const auto first = static_cast<std::uint32_t>(mesh.vertices.size());
                    mesh.vertices.push_back({x, shift, z});
                    mesh.vertices.push_back({x + h, shift, z});
                    mesh.vertices.push_back(
                        {x + (body ? h * .2 : 0), shift + h, z + h * (body ? -.4 : .3)});
                    mesh.triangles.push_back({{first, first + 1, first + 2}, body + 1});
                }
                const auto space = galerkin::buildSpace(mesh);
                for (unsigned signs = 0; signs < 4; ++signs) {
                    auto x = space.panels[0], y = space.panels[1];
                    if (signs & 1) {
                        x.normal.x = -x.normal.x;
                        x.normal.y = -x.normal.y;
                        x.normal.z = -x.normal.z;
                    }
                    if (signs & 2) {
                        y.normal.x = -y.normal.x;
                        y.normal.y = -y.normal.y;
                        y.normal.z = -y.normal.z;
                    }
                    const auto refined = reference(x, y, k, duffy(20));
                    const auto convergence = reference(x, y, k, duffy(16));
                    compare(convergence, refined, 2e-8,
                            "Independent Duffy oracle did not converge");
                    require(galerkin::qualifiedRegularDegree(galerkin::makeRegularPanelBounds(x),
                                                             galerkin::makeRegularPanelBounds(y), k,
                                                             1e-4, degree) == degree,
                            "Regular fixture did not qualify for requested degree");
                    const auto rule = galerkin::triangleRule(degree);
                    const auto same = reference(x, y, k, rule);
                    const auto fp64 = gpu(x, y, k, rule, false), mixed = gpu(x, y, k, rule, true);
                    compare(fp64, same, 1e-9, "FP64 device formula differs from CPU oracle");
                    const double pressureError = groupError(mixed, same, 0, 9);
                    const double fluxError = groupError(mixed, same, 9, 18);
                    const double rhsError = groupError(mixed, same, 18, 21);
                    const double greenError = groupError(mixed, same, 21, 22);
                    if (std::max({pressureError, fluxError, rhsError, greenError}) >= 1e-6) {
                        std::cerr << "mixed_failure kr=" << kr << " scale=" << scale
                                  << " degree=" << degree << " signs=" << signs
                                  << " pressure=" << pressureError << " flux=" << fluxError
                                  << " rhs=" << rhsError << " green=" << greenError << '\n';
                        throw std::runtime_error("Mixed device arithmetic exceeds 1e-6");
                    }
                    compare(same, refined, 2e-5, "Regular quadrature exceeds 0.2 tolerance budget");
                    for (unsigned i = 0; i < 3; ++i) {
                        C sum{};
                        for (unsigned j = 0; j < 3; ++j)
                            sum += mixed[9 + 3 * i + j];
                        require(std::abs(sum - mixed[18 + i]) <=
                                    1e-12 * std::max(std::abs(mixed[18 + i]), 1e-300),
                                "Interface P1 flux does not sum to DP0 RHS");
                    }
                    std::array<C, 9> scatteredActual{}, scatteredExpected{};
                    double scatteredDifferenceNorm = 0, scatteredExpectedNorm = 0;
                    for (unsigned i = 0; i < 3; ++i)
                        for (unsigned j = 0; j < 3; ++j) {
                            const C curl =
                                C{0, -1} * test::wideband_oracle::dot(x.curls[i], y.curls[j]) / k;
                            const auto index = 3 * i + j;
                            scatteredActual[index] = mixed[index] + curl * mixed[21];
                            scatteredExpected[index] = same[index] + curl * same[21];
                            scatteredDifferenceNorm +=
                                std::norm(scatteredActual[index] - scatteredExpected[index]);
                            scatteredExpectedNorm += std::norm(scatteredExpected[index]);
                        }
                    for (unsigned i = 0; i < 3; ++i)
                        for (unsigned j = 0; j < 3; ++j) {
                            const auto index = 3 * i + j;
                            const double difference =
                                std::abs(scatteredActual[index] - scatteredExpected[index]);
                            const double expectedMagnitude = std::abs(scatteredExpected[index]);
                            if (!(difference <= 1e-6 * std::max(expectedMagnitude, 1e-300))) {
                                std::cerr
                                    << "curl_failure kr=" << kr << " scale=" << scale
                                    << " degree=" << degree << " signs=" << signs << " i=" << i
                                    << " j=" << j << " abs_error=" << difference
                                    << " abs_expected=" << expectedMagnitude
                                    << " vector_relative_error="
                                    << std::sqrt(scatteredDifferenceNorm / scatteredExpectedNorm)
                                    << '\n';
                                throw std::runtime_error(
                                    "Physical curl-scattered pressure differs from FP64");
                            }
                        }
                }
            }
}
} // namespace
int main() {
    try {
        run();
        std::cout << "Subdomain integral CUDA oracle PASS.\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "Subdomain integral CUDA oracle failed: " << e.what() << '\n';
        return 1;
    }
}
