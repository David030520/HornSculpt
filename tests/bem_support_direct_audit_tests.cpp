#include "bem_hf_test_oracle.hpp"
#include "bem_support_operator_oracle.hpp"
#include <cuda/support_direct_audit.hpp>
#include <iostream>

namespace {
using namespace hornsculpt;
using namespace hornsculpt::wideband;
using C = std::complex<double>;
using Z = cuDoubleComplex;
void check(cudaError_t status) {
    if (status != cudaSuccess)
        throw std::runtime_error(cudaGetErrorString(status));
}
void require(bool value, const char *message) {
    if (!value)
        throw std::runtime_error(message);
}
class Device {
  public:
    explicit Device(std::size_t count) : count_(count) {
        check(cudaMalloc(reinterpret_cast<void **>(&data), count * sizeof(Z)));
    }
    ~Device() {
        cudaFree(data);
    }
    Device(const Device &) = delete;
    Device &operator=(const Device &) = delete;
    void upload(std::span<const C> values, cudaStream_t stream) {
        require(values.size() == count_, "Direct audit fixture input size mismatch.");
        std::vector<Z> packed;
        for (const auto value : values)
            packed.push_back(make_cuDoubleComplex(value.real(), value.imag()));
        check(cudaMemcpyAsync(data, packed.data(), count_ * sizeof(Z), cudaMemcpyHostToDevice,
                              stream));
        check(cudaStreamSynchronize(stream));
    }
    std::vector<C> download(cudaStream_t stream) {
        std::vector<Z> packed(count_);
        check(cudaMemcpyAsync(packed.data(), data, count_ * sizeof(Z), cudaMemcpyDeviceToHost,
                              stream));
        check(cudaStreamSynchronize(stream));
        std::vector<C> result;
        for (const auto value : packed)
            result.emplace_back(value.x, value.y);
        return result;
    }
    Z *data{};

  private:
    std::size_t count_{};
};
void checkPartition(unsigned capacity) {
    constexpr double k = 20;
    auto fixture = test::makeWidebandFarBlockFixture(k);
    for (auto &panel : fixture.space.panels)
        panel.boundary = BoundaryCondition::driven;
    SupportTreeSettings settings;
    settings.leafCapacity = capacity;
    settings.nearDistanceRatio = .5;
    auto tree = makeSupportTree(fixture.space, settings);
    const std::array<std::uint32_t, 2> faces{5, 1};
    auto plan = makeSupportNearPlan(fixture.space, tree, faces, 16 * 1024 * 1024);
    SupportNearOperatorSettings quadrature;
    quadrature.quadrature.nearDistanceRatio = .5;
    auto near = std::make_shared<SupportNearOperator>(fixture.space, tree, plan, k, quadrature);
    SupportDirectAudit direct(fixture.space, tree, plan, near, k, 4);
    const auto reference =
        test::support_operator_oracle::directReference(fixture.space, k, quadrature.quadrature);
    Device p(8), q(2), output(8);
    for (unsigned rhs = 0; rhs < 2; ++rhs) {
        std::vector<C> pressure;
        for (unsigned i = 0; i < 8; ++i)
            pressure.emplace_back(std::sin(.7 * i + .3 * rhs), std::cos(.4 * i - .2 * rhs));
        const std::vector<C> flux{{.4 + .2 * rhs, -.3}, {-.7, .1 + .5 * rhs}};
        p.upload(pressure, near->stream());
        q.upload(flux, near->stream());
        for (unsigned channel = 0; channel < 2; ++channel) {
            if (channel)
                direct.applyFaceFluxDevice(q.data, output.data);
            else
                direct.applyPressureDevice(p.data, output.data);
            const auto actual = output.download(near->stream());
            const auto &input = channel ? flux : pressure,
                       &matrix = channel ? reference.b : reference.a;
            double error{}, norm{};
            for (unsigned row = 0; row < 8; ++row) {
                C expected{};
                for (unsigned column = 0; column < input.size(); ++column)
                    expected += matrix[row * input.size() + column] * input[column];
                error += std::norm(actual[row] - expected);
                norm += std::norm(expected);
            }
            const double relative = std::sqrt(error / norm);
            std::cout << "Direct weak Galerkin audit capacity=" << capacity << ", rhs=" << rhs
                      << ", A/B channel=" << channel << ", relative=" << relative << '\n';
            require(relative < 2e-7,
                    "Original direct Galerkin audit failed independent matrix action.");
        }
    }
    std::vector<C> zero(8);
    p.upload(zero, near->stream());
    direct.applyPressureDevice(p.data, output.data);
    for (const auto value : output.download(near->stream()))
        require(value == C{}, "Direct audit retained a previous pressure vector.");
}
} // namespace
int main() {
    int count{};
    if (cudaGetDeviceCount(&count) != cudaSuccess || !count)
        return 77;
    try {
        for (const unsigned capacity : {1, 2, 64})
            checkPartition(capacity);
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
