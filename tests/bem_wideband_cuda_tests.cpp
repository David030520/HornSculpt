#include "bem_diagonal_test_operator.hpp"
#include "bem_hf_test_oracle.hpp"
#include "bem_mixed_wideband_oracle.hpp"
#include "bem_support_operator_oracle.hpp"
#include "cuda/support_fmm_operator.hpp"
#include "cuda/support_m2l_graph.hpp"
#include "cuda/wideband_operator.hpp"
#include "fixture_geometry.hpp"
#include "fmm/rcr_work.hpp"
#include "fmm/translation_plan_cache.hpp"
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cuda_runtime_api.h>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numbers>
#include <set>
#include <stdexcept>
#include <stop_token>
#include <string>
#include <string_view>
#include <thread>
#include <utility>

void checkScaledRegularCuda();

namespace {
using namespace hornsculpt;
using namespace hornsculpt::wideband;
void require(bool condition, const char *message) {
    if (!condition)
        throw std::runtime_error(message);
}
double dot(Vec3 a, Vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
Complex wave(double argument) {
    return {std::cos(argument), std::sin(argument)};
}
Complex green(Vec3 delta, double k) {
    const double radius = std::sqrt(dot(delta, delta));
    return wave(-k * radius) / (4 * std::numbers::pi * radius);
}
void checkDiagonalLevel() {
    constexpr double k = 20;
    const auto rule = makeSphereRule(32);
    const std::vector<Vec3> translations{{3, 1, 0}, {-3, 1, 0}, {0, -3, 1}};
    const std::vector<Vec3> sourceOffsets{
        {.031, -.018, .009}, {-.024, .011, .017}, {.012, .021, -.013}, {-.007, -.014, -.019}};
    const std::vector<Complex> amplitudes{{.8, -.3}, {-.4, .9}, {1.1, .2}, {-.6, -.7}};
    const std::vector<Vec3> targetOffsets{{0, 0, 0}, {.017, .013, -.021}, {-.023, .009, .015}};
    DiagonalLevelPlan plan;
    plan.sourceNodes = static_cast<unsigned>(sourceOffsets.size());
    plan.targetNodes = static_cast<unsigned>(targetOffsets.size());
    plan.angularSamples = static_cast<unsigned>(rule.samples.size());
    plan.translationClasses = static_cast<unsigned>(translations.size());
    plan.targetOffsets = {0, 0, 47, 92};
    for (unsigned i = 0; i < 92; ++i)
        plan.interactions.push_back(
            {i % plan.sourceNodes, (2 * i + i / 7) % plan.translationClasses});
    for (const auto delta : translations) {
        const auto coefficients = diagonalTranslation(rule, 16, k, delta);
        plan.translationCoefficients.insert(plan.translationCoefficients.end(),
                                            coefficients.begin(), coefficients.end());
    }
    std::vector<Complex> source(static_cast<std::size_t>(plan.sourceNodes) * plan.angularSamples);
    for (unsigned node = 0; node < plan.sourceNodes; ++node)
        for (unsigned q = 0; q < plan.angularSamples; ++q)
            source[static_cast<std::size_t>(node) * plan.angularSamples + q] =
                amplitudes[node] * wave(k * dot(rule.samples[q].direction, sourceOffsets[node]));
    DiagonalM2LOperator gpu(plan);
    const auto actual = gpu.apply(source);
    std::vector<Complex> expected(actual.size());
    double coefficientError = 0, coefficientNorm = 0;
    for (unsigned target = 0; target < plan.targetNodes; ++target)
        for (unsigned q = 0; q < plan.angularSamples; ++q) {
            Complex value{};
            for (auto pair = plan.targetOffsets[target]; pair < plan.targetOffsets[target + 1];
                 ++pair) {
                const auto job = plan.interactions[pair];
                value +=
                    plan.translationCoefficients[static_cast<std::size_t>(job.translationClass) *
                                                     plan.angularSamples +
                                                 q] *
                    source[static_cast<std::size_t>(job.source) * plan.angularSamples + q];
            }
            const auto flat = static_cast<std::size_t>(target) * plan.angularSamples + q;
            expected[flat] = value;
            coefficientError += std::norm(actual[flat] - value);
            coefficientNorm += std::norm(value);
        }
    const double relativeM2L = std::sqrt(coefficientError / coefficientNorm);
    require(relativeM2L < 2e-13, "GPU diagonal class M2L disagrees with the independent CPU sum.");
    double fieldError = 0, fieldNorm = 0;
    for (unsigned target = 0; target < plan.targetNodes; ++target) {
        Complex field{}, exact{};
        for (unsigned q = 0; q < plan.angularSamples; ++q)
            field += rule.samples[q].weight *
                     actual[static_cast<std::size_t>(target) * plan.angularSamples + q] *
                     wave(-k * dot(rule.samples[q].direction, targetOffsets[target]));
        for (auto pair = plan.targetOffsets[target]; pair < plan.targetOffsets[target + 1];
             ++pair) {
            const auto job = plan.interactions[pair];
            const auto center = translations[job.translationClass], y = sourceOffsets[job.source],
                       x = targetOffsets[target];
            exact += amplitudes[job.source] *
                     green({center.x + x.x - y.x, center.y + x.y - y.y, center.z + x.z - y.z}, k);
        }
        fieldError += std::norm(field - exact);
        fieldNorm += std::norm(exact);
    }
    const double relativeField = std::sqrt(fieldError / fieldNorm);
    require(relativeField < 2e-11, "GPU diagonal M2L fails the closed-form Helmholtz field.");
    const Complex scale{-.2, .7};
    for (auto &value : source)
        value *= scale;
    const auto repeated = gpu.apply(source);
    double reuseError = 0;
    for (std::size_t j = 0; j < repeated.size(); ++j)
        reuseError += std::norm(repeated[j] - scale * expected[j]);
    require(std::sqrt(reuseError / coefficientNorm) < 2e-13,
            "GPU diagonal M2L retained stale coefficients across applications.");
    std::fill(source.begin(), source.end(), Complex{});
    const auto zero = gpu.apply(source);
    require(std::all_of(zero.begin(), zero.end(), [](auto value) { return value == Complex{}; }),
            "GPU diagonal M2L did not overwrite the previous local spectrum.");
    std::stop_source stop;
    stop.request_stop();
    bool stopped = false;
    try {
        static_cast<void>(gpu.apply(source, stop.get_token()));
    } catch (const std::runtime_error &error) {
        stopped = std::string(error.what()).find("cancelled") != std::string::npos;
    }
    require(stopped, "GPU diagonal M2L ignored a requested cancellation.");
    std::cout << "Real CUDA diagonal M2L: relative coefficients=" << relativeM2L
              << ", closed Green field=" << relativeField << ", bytes=" << gpu.deviceBytes()
              << '\n';
}
void checkAngularTransfers() {
    struct Case {
        unsigned sourceOrder, targetOrder, bandlimit;
    };
    for (const auto spec : {Case{0, 0, 0}, Case{5, 9, 5}, Case{9, 5, 5}, Case{16, 24, 12}}) {
        const auto sourceRule = makeSphereRule(spec.sourceOrder),
                   targetRule = makeSphereRule(spec.targetOrder);
        const auto plan = makeAngularTransferPlan(sourceRule, targetRule, spec.bandlimit);
        constexpr unsigned batches = 3;
        std::vector<Complex> values, expected;
        for (unsigned batch = 0; batch < batches; ++batch) {
            std::vector<Complex> single;
            for (std::size_t q = 0; q < sourceRule.samples.size(); ++q) {
                const auto d = sourceRule.samples[q].direction;
                single.emplace_back(std::sin(.17 * q + .4 * batch) + d.x * (.3 + batch),
                                    std::cos(.23 * q - .7 * batch) + d.y * (.9 - batch));
            }
            const auto reference = resampleSpectrum(single, plan);
            values.insert(values.end(), single.begin(), single.end());
            expected.insert(expected.end(), reference.begin(), reference.end());
        }
        AngularTransferOperator gpu(plan, batches);
        const auto actual = gpu.apply(values);
        double error = 0, norm = 0;
        require(actual.size() == expected.size(), "GPU angular transfer returned the wrong shape.");
        for (std::size_t q = 0; q < actual.size(); ++q) {
            error += std::norm(actual[q] - expected[q]);
            norm += std::norm(expected[q]);
        }
        const double relative = std::sqrt(error / norm);
        std::cout << "Real CUDA angular " << spec.sourceOrder << " -> " << spec.targetOrder
                  << ", bandlimit=" << spec.bandlimit << ", batches=" << batches
                  << ", relative=" << relative << ", bytes=" << gpu.deviceBytes() << '\n';
        require(relative < 2e-12, "GPU FFT/Legendre transfer fails the independent CPU transform.");
        std::fill(values.begin(), values.end(), Complex{});
        const auto zero = gpu.apply(values);
        for (const auto value : zero)
            require(value == Complex{},
                    "GPU angular transfer retained a previous batched spectrum.");
    }
}

void checkAngularCenterShifts() {
    const auto childRule = makeSphereRule(8), parentRule = makeSphereRule(12);
    const auto upPlan = makeAngularTransferPlan(childRule, parentRule, 8),
               downPlan = makeAngularTransferPlan(parentRule, childRule, 8);
    AngularTransferOperator upward(upPlan, 1), downward(downPlan, 1);
    const Vec3 shift{.011, -.008, .006};
    constexpr double k = 20;
    std::vector<Complex> outgoing, incoming;
    for (std::size_t q = 0; q < childRule.samples.size(); ++q) {
        const auto d = childRule.samples[q].direction;
        outgoing.emplace_back(std::sin(.11 * q) + d.x, std::cos(.17 * q) - d.z);
    }
    for (std::size_t q = 0; q < parentRule.samples.size(); ++q) {
        const auto d = parentRule.samples[q].direction;
        incoming.emplace_back(std::cos(.19 * q) - d.y, std::sin(.13 * q) + d.z);
    }
    const auto actualUp = upward.outgoingToParent(outgoing, k, shift),
               actualDown = downward.incomingToChild(incoming, k, shift),
               referenceUp = outgoingToParent(outgoing, childRule, parentRule, 8, k, shift),
               referenceDown = incomingToChild(incoming, parentRule, childRule, 8, k, shift);
    double error = 0, norm = 0;
    Complex left{}, right{};
    for (std::size_t q = 0; q < actualUp.size(); ++q) {
        error += std::norm(actualUp[q] - referenceUp[q]);
        norm += std::norm(referenceUp[q]);
        left += parentRule.samples[q].weight * std::conj(actualUp[q]) * incoming[q];
    }
    for (std::size_t q = 0; q < actualDown.size(); ++q) {
        error += std::norm(actualDown[q] - referenceDown[q]);
        norm += std::norm(referenceDown[q]);
        right += childRule.samples[q].weight * std::conj(outgoing[q]) * actualDown[q];
    }
    const double relative = std::sqrt(error / norm),
                 adjoint = std::abs(left - right) / (1 + std::abs(left) + std::abs(right));
    std::cout << "Real CUDA shifted angular transfers relative=" << relative
              << ", weighted Hermitian adjoint=" << adjoint << '\n';
    require(relative < 2e-12, "GPU center-shifted angular transfer fails the CPU reference.");
    require(adjoint < 2e-12, "GPU upward/downward shifts violate the weighted adjoint contract.");
}

class TestDeviceArray {
  public:
    explicit TestDeviceArray(std::size_t count) : count_(count) {
        if (cudaMalloc(reinterpret_cast<void **>(&values), count * sizeof(cuDoubleComplex)) !=
            cudaSuccess)
            throw std::runtime_error("CUDA hierarchy test allocation failed.");
    }
    ~TestDeviceArray() {
        cudaFree(values);
    }
    TestDeviceArray(const TestDeviceArray &) = delete;
    TestDeviceArray &operator=(const TestDeviceArray &) = delete;
    void upload(const std::vector<Complex> &input, cudaStream_t stream) {
        require(input.size() == count_, "CUDA hierarchy test upload shape is invalid.");
        std::vector<cuDoubleComplex> packed;
        for (const auto value : input)
            packed.push_back(make_cuDoubleComplex(value.real(), value.imag()));
        require(cudaMemcpyAsync(values, packed.data(), count_ * sizeof(cuDoubleComplex),
                                cudaMemcpyHostToDevice, stream) == cudaSuccess,
                "CUDA hierarchy upload failed.");
        require(cudaStreamSynchronize(stream) == cudaSuccess,
                "CUDA hierarchy upload completion failed.");
    }
    std::vector<Complex> download() const {
        std::vector<cuDoubleComplex> packed(count_);
        require(cudaMemcpy(packed.data(), values, count_ * sizeof(cuDoubleComplex),
                           cudaMemcpyDeviceToHost) == cudaSuccess,
                "CUDA hierarchy download failed.");
        std::vector<Complex> output;
        for (const auto value : packed)
            output.emplace_back(value.x, value.y);
        return output;
    }
    cuDoubleComplex *values{};

  private:
    std::size_t count_{};
};
class TestDeviceStatus {
  public:
    TestDeviceStatus() {
        require(cudaMalloc(reinterpret_cast<void **>(&value), sizeof(int)) == cudaSuccess,
                "CUDA test failure status allocation failed.");
    }
    ~TestDeviceStatus() {
        cudaFree(value);
    }
    TestDeviceStatus(const TestDeviceStatus &) = delete;
    TestDeviceStatus &operator=(const TestDeviceStatus &) = delete;
    void reset(cudaStream_t stream) {
        require(cudaMemsetAsync(value, 0, sizeof(int), stream) == cudaSuccess,
                "CUDA test failure status reset failed.");
    }
    int finish(cudaStream_t stream) const {
        int result = 0;
        require(cudaMemcpyAsync(&result, value, sizeof(int), cudaMemcpyDeviceToHost, stream) ==
                    cudaSuccess,
                "CUDA test failure status download failed.");
        require(cudaStreamSynchronize(stream) == cudaSuccess,
                "CUDA test failure status completion failed.");
        return result;
    }
    int *value{};
};

void checkDeviceHierarchy() {
    const auto childRule = makeSphereRule(8), middleRule = makeSphereRule(16),
               farRule = makeSphereRule(32);
    const Vec3 sourceLeafToMiddle{.01, -.02, .02}, sourceMiddleToFar{.03, -.02, .01},
        targetMiddleToFar{-.04, .03, -.01}, targetLeafToMiddle{-.02, -.01, .02},
        sourcePoint{.005, .008, -.004}, targetPoint{.003, -.004, .006}, separation{3, 1, -.5};
    constexpr double k = 20;
    const Complex amplitude{1.2, -.4};
    AngularTransferOperator upFirst(makeAngularTransferPlan(childRule, middleRule, 8), 1);
    const auto sharedStream = upFirst.stream();
    AngularTransferOperator upSecond(makeAngularTransferPlan(middleRule, farRule, 16), 1, {},
                                     sharedStream),
        downFirst(makeAngularTransferPlan(farRule, middleRule, 16), 1, {}, sharedStream),
        downSecond(makeAngularTransferPlan(middleRule, childRule, 8), 1, {}, sharedStream);
    DiagonalLevelPlan plan;
    plan.sourceNodes = plan.targetNodes = plan.translationClasses = 1;
    plan.angularSamples = static_cast<unsigned>(farRule.samples.size());
    plan.targetOffsets = {0, 1};
    plan.interactions = {{0, 0}};
    plan.translationCoefficients = diagonalTranslation(farRule, 16, k, separation);
    DiagonalM2LOperator m2l(plan, {}, sharedStream);
    TestDeviceArray source(childRule.samples.size()), sourceMiddle(middleRule.samples.size()),
        sourceFar(farRule.samples.size()), localFar(farRule.samples.size()),
        localMiddle(middleRule.samples.size()), localChild(childRule.samples.size());
    std::vector<Complex> emitted;
    for (const auto sample : childRule.samples)
        emitted.push_back(amplitude * wave(k * dot(sample.direction, sourcePoint)));
    source.upload(emitted, sharedStream);
    // Every intermediate array remains on the device. These two upload/download
    // boundaries surround the complete two-level up/M2L/two-level down chain.
    upFirst.outgoingToParentDevice(source.values, sourceMiddle.values, k, sourceLeafToMiddle);
    upSecond.outgoingToParentDevice(sourceMiddle.values, sourceFar.values, k, sourceMiddleToFar);
    m2l.applyDevice(sourceFar.values, localFar.values);
    downFirst.incomingToChildDevice(localFar.values, localMiddle.values, k, targetMiddleToFar);
    downSecond.incomingToChildDevice(localMiddle.values, localChild.values, k, targetLeafToMiddle);
    const auto local = localChild.download();
    Complex field{};
    std::array<Complex, 3> gradient{};
    for (std::size_t q = 0; q < local.size(); ++q) {
        const auto sample = childRule.samples[q];
        const auto value = local[q] * sample.weight * wave(-k * dot(sample.direction, targetPoint));
        field += value;
        gradient[0] += Complex{0, -k * sample.direction.x} * value;
        gradient[1] += Complex{0, -k * sample.direction.y} * value;
        gradient[2] += Complex{0, -k * sample.direction.z} * value;
    }
    const Vec3 delta{separation.x + targetMiddleToFar.x + targetLeafToMiddle.x + targetPoint.x -
                         sourceMiddleToFar.x - sourceLeafToMiddle.x - sourcePoint.x,
                     separation.y + targetMiddleToFar.y + targetLeafToMiddle.y + targetPoint.y -
                         sourceMiddleToFar.y - sourceLeafToMiddle.y - sourcePoint.y,
                     separation.z + targetMiddleToFar.z + targetLeafToMiddle.z + targetPoint.z -
                         sourceMiddleToFar.z - sourceLeafToMiddle.z - sourcePoint.z};
    const auto exact = amplitude * green(delta, k);
    const double r = std::sqrt(dot(delta, delta));
    const auto radial = exact * Complex{-1 / r, -k};
    double gradientError = 0, gradientNorm = 0;
    const double axes[]{delta.x, delta.y, delta.z};
    for (unsigned axis = 0; axis < 3; ++axis) {
        const auto expected = radial * (axes[axis] / r);
        gradientError += std::norm(gradient[axis] - expected);
        gradientNorm += std::norm(expected);
    }
    const double relativeField = std::abs(field - exact) / std::abs(exact),
                 relativeGradient = std::sqrt(gradientError / gradientNorm);
    std::cout << "Real CUDA resident two-level off-axis hierarchy field=" << relativeField
              << ", gradient=" << relativeGradient << '\n';
    require(relativeField < 2e-10, "CUDA resident hierarchy fails the exact off-axis Green field.");
    require(relativeGradient < 2e-10,
            "CUDA resident hierarchy fails the exact off-axis Green gradient.");
}

void checkMixedDeviceHierarchy() {
    constexpr unsigned order = 12;
    constexpr double k = 20, radius = .025;
    const auto middleRule = makeSphereRule(32), farRule = makeSphereRule(48);
    const Vec3 sourceLeafToMiddle{.01, -.02, .02}, sourceMiddleToFar{.03, -.02, .01},
        targetMiddleToFar{-.04, .03, -.01}, targetLeafToMiddle{-.02, -.01, .02},
        sourcePoint{.005, .008, -.004}, targetPoint{.003, -.004, .006}, separation{3, 1, -.5};
    const Complex amplitude{1.2, -.4};
    AngularTransferOperator sourceHandoff(makeAngularTransferPlan(middleRule, middleRule, order),
                                          1);
    const auto stream = sourceHandoff.stream();
    AngularTransferOperator up(makeAngularTransferPlan(middleRule, farRule, 32), 1, {}, stream),
        down(makeAngularTransferPlan(farRule, middleRule, 32), 1, {}, stream),
        targetHandoff(makeAngularTransferPlan(middleRule, middleRule, order), 1, {}, stream);
    DiagonalLevelPlan plan;
    plan.sourceNodes = plan.targetNodes = plan.translationClasses = 1;
    plan.angularSamples = static_cast<unsigned>(farRule.samples.size());
    plan.targetOffsets = {0, 1};
    plan.interactions = {{0, 0}};
    plan.translationCoefficients = diagonalTranslation(farRule, 24, k, separation);
    DiagonalM2LOperator m2l(plan, {}, stream);
    const auto basis = scaledRegularWaves(order, k, radius, sourcePoint);
    std::vector<Complex> moments;
    for (const auto &value : basis.waves)
        moments.push_back(amplitude * std::conj(value.value));
    TestDeviceArray leaf(moments.size()), sourceMiddle(middleRule.samples.size()),
        sourceFar(farRule.samples.size()), localFar(farRule.samples.size()),
        localMiddle(middleRule.samples.size()), localLeaf(moments.size());
    leaf.upload(moments, stream);
    sourceHandoff.scaledOutgoingToParentDevice(leaf.values, sourceMiddle.values, k, radius,
                                               sourceLeafToMiddle);
    up.outgoingToParentDevice(sourceMiddle.values, sourceFar.values, k, sourceMiddleToFar);
    m2l.applyDevice(sourceFar.values, localFar.values);
    down.incomingToChildDevice(localFar.values, localMiddle.values, k, targetMiddleToFar);
    targetHandoff.incomingToScaledChildDevice(localMiddle.values, localLeaf.values, k, radius,
                                              targetLeafToMiddle);
    const ScaledLocal local{order, radius, localLeaf.download()};
    const auto actual = evaluateScaledLocal(local, k, targetPoint);
    const Vec3 delta{separation.x + targetMiddleToFar.x + targetLeafToMiddle.x + targetPoint.x -
                         sourceMiddleToFar.x - sourceLeafToMiddle.x - sourcePoint.x,
                     separation.y + targetMiddleToFar.y + targetLeafToMiddle.y + targetPoint.y -
                         sourceMiddleToFar.y - sourceLeafToMiddle.y - sourcePoint.y,
                     separation.z + targetMiddleToFar.z + targetLeafToMiddle.z + targetPoint.z -
                         sourceMiddleToFar.z - sourceLeafToMiddle.z - sourcePoint.z};
    const auto exact = amplitude * green(delta, k);
    const double distance = std::hypot(delta.x, delta.y, delta.z);
    const auto radial = exact * Complex{-1 / distance, -k};
    const double coordinate[]{delta.x, delta.y, delta.z};
    double error = 0, norm = 0;
    for (unsigned component = 0; component < 3; ++component) {
        const auto expected = radial * (coordinate[component] / distance);
        error += std::norm(actual.gradient[component] - expected);
        norm += std::norm(expected);
    }
    const double field = std::abs(actual.value - exact) / std::abs(exact),
                 gradient = std::sqrt(error / norm);
    std::cout << "CUDA resident LF/HF mixed hierarchy field=" << field << ", gradient=" << gradient
              << '\n';
    require(field < 2e-10 && gradient < 2e-10,
            "CUDA LF/HF handoff hierarchy fails exact off-axis Green physics.");
}

void checkScaledAxial() {
    for (const auto kind : {ScaledAxialKind::multipoleToParent, ScaledAxialKind::multipoleToLocal,
                            ScaledAxialKind::localToChild}) {
        const double sourceRadius = kind == ScaledAxialKind::localToChild ? .04 : .01,
                     targetRadius = kind == ScaledAxialKind::multipoleToParent ? .04 : .01,
                     distance = kind == ScaledAxialKind::multipoleToLocal ? .2 : .02;
        const auto plan =
            makeScaledAxialPlan(kind, 6, 8, 10, sourceRadius, targetRadius, distance, 48);
        constexpr unsigned batches = 3;
        constexpr std::size_t sourceCount = 49, targetCount = 81;
        std::vector<Complex> input, expected;
        for (unsigned batch = 0; batch < batches; ++batch) {
            std::vector<Complex> one;
            for (unsigned l = 0; l <= 6; ++l)
                for (int m = -static_cast<int>(l); m <= static_cast<int>(l); ++m)
                    one.emplace_back(std::pow(.31, l) * std::cos(.3 * m + batch),
                                     std::pow(.29, l) * std::sin(.7 * m - .4 * batch));
            const auto reference = applyScaledAxialPlan(plan, one);
            input.insert(input.end(), one.begin(), one.end());
            expected.insert(expected.end(), reference.begin(), reference.end());
        }
        ScaledAxialOperator gpu(plan, batches);
        TestDeviceArray source(batches * sourceCount), target(batches * targetCount);
        source.upload(input, gpu.stream());
        gpu.applyDevice(source.values, target.values);
        const auto actual = target.download();
        double error = 0, norm = 0;
        for (std::size_t j = 0; j < actual.size(); ++j) {
            error += std::norm(actual[j] - expected[j]);
            norm += std::norm(expected[j]);
        }
        const double relative = std::sqrt(error / norm);
        std::cout << "Real CUDA scaled axial kind=" << static_cast<unsigned>(kind)
                  << ", relative=" << relative << ", bytes=" << gpu.deviceBytes() << '\n';
        require(relative < 2e-13, "CUDA scaled axial m-block application fails the CPU reference.");
        std::fill(input.begin(), input.end(), Complex{});
        source.upload(input, gpu.stream());
        gpu.applyDevice(source.values, target.values);
        for (const auto value : target.download())
            require(value == Complex{},
                    "CUDA scaled axial translation retained previous coefficients.");
    }
}

void checkCooperativeContractions() {
    constexpr unsigned capacity = 19;
    const Complex sentinel{-.41, .73};
    const double nan = std::numeric_limits<double>::quiet_NaN();
    SphericalRotationOperator streamOwner(makeSphericalRotationPlan(0, {0, 0, 1}), 1);
    require(streamOwner.stream() != nullptr,
            "Cooperative contraction gate requires a nondefault CUDA stream.");
    const auto check = [&](const std::vector<Complex> &actual,
                           const std::vector<Complex> &expected) {
        double error = 0, norm = 0;
        for (std::size_t j = 0; j < expected.size(); ++j) {
            error += std::norm(actual[j] - expected[j]);
            norm += std::norm(expected[j]);
        }
        require(norm == 0 ? error == 0 : std::sqrt(error / norm) < 2e-12,
                "Cooperative FP64 contraction differs from the independent CPU row sum.");
        for (std::size_t j = expected.size(); j < actual.size(); ++j)
            require(actual[j] == sentinel,
                    "Cooperative contraction overwrote an inactive trailing batch.");
    };
    TestDeviceStatus status;
    for (unsigned order : {0u, 1u, 8u, 17u, 33u, 65u}) {
        const auto plan = makeSphericalRotationPlan(order, {7, -5, 3});
        const auto fullState = static_cast<std::size_t>(order + 1) * (order + 1);
        SphericalRotationOperator rotation(plan, capacity, {}, streamOwner.stream());
        TestDeviceArray source(capacity * fullState), output(capacity * fullState);
        const auto identity = rotation.deviceIdentity();
        const auto bytes = rotation.deviceBytes();
        for (unsigned batches : {1u, 2u, 17u}) {
            std::vector<Complex> input(capacity * fullState, {nan, nan}), expected;
            for (unsigned batch = 0; batch < batches; ++batch) {
                std::vector<Complex> one(fullState);
                for (std::size_t j = 0; j < fullState; ++j)
                    input[batch * fullState + j] =
                        one[j] = {std::sin(.13 * j + batch), std::cos(.27 * j - .31 * batch)};
                const auto reference = rotateToAxis(plan, one);
                expected.insert(expected.end(), reference.begin(), reference.end());
            }
            source.upload(input, rotation.stream());
            output.upload(std::vector<Complex>(capacity * fullState, sentinel), rotation.stream());
            status.reset(rotation.stream());
            rotation.enqueueToAxisDeviceBatch(source.values, output.values, batches, status.value);
            require(status.finish(rotation.stream()) == 0,
                    "Cooperative forward rotation read poisoned inactive batches.");
            check(output.download(), expected);

            const auto prefixOrder = order ? order - 1 : 0;
            const auto state = static_cast<std::size_t>(prefixOrder + 1) * (prefixOrder + 1);
            input.assign(capacity * fullState, {nan, nan});
            expected.clear();
            for (unsigned batch = 0; batch < batches; ++batch) {
                std::vector<Complex> one(fullState);
                for (std::size_t j = 0; j < state; ++j)
                    input[batch * state + j] =
                        one[j] = {std::cos(.19 * j + batch), std::sin(.11 * j - .43 * batch)};
                const auto reference = rotateFromAxis(plan, one);
                expected.insert(expected.end(), reference.begin(), reference.begin() + state);
            }
            source.upload(input, rotation.stream());
            output.upload(std::vector<Complex>(capacity * fullState, sentinel), rotation.stream());
            rotation.fromAxisPrefixDeviceBatch(source.values, output.values, prefixOrder, batches);
            check(output.download(), expected);
            const auto reused = rotation.deviceIdentity();
            require(bytes == rotation.deviceBytes() && reused.offsets == identity.offsets &&
                        reused.matrix == identity.matrix && reused.phases == identity.phases,
                    "Cooperative rotation replaced immutable storage during an application.");
        }
    }
    struct AxialCase {
        unsigned sourceOrder, targetOrder;
    };
    for (const auto spec : {AxialCase{0, 0}, AxialCase{0, 17}, AxialCase{17, 0}, AxialCase{17, 33},
                            AxialCase{33, 17}, AxialCase{128, 33}}) {
        ScaledAxialPlan plan;
        plan.sourceOrder = spec.sourceOrder;
        plan.targetOrder = spec.targetOrder;
        // This gate isolates the executor from plan quadrature. Existing
        // physical M2M/M2L/L2L gates independently cover coefficient creation.
        for (unsigned m = 0; m <= std::min(spec.sourceOrder, spec.targetOrder); ++m) {
            plan.blockOffsets.push_back(plan.coefficients.size());
            for (unsigned l = m; l <= spec.targetOrder; ++l)
                for (unsigned n = m; n <= spec.sourceOrder; ++n)
                    plan.coefficients.emplace_back(std::cos(.13 * l + .21 * n + m),
                                                   std::sin(.17 * l - .31 * n + m));
        }
        plan.blockOffsets.push_back(plan.coefficients.size());
        const auto sourceState =
                       static_cast<std::size_t>(spec.sourceOrder + 1) * (spec.sourceOrder + 1),
                   targetState =
                       static_cast<std::size_t>(spec.targetOrder + 1) * (spec.targetOrder + 1);
        ScaledAxialOperator axial(plan, capacity, {}, streamOwner.stream());
        TestDeviceArray source(capacity * sourceState), output(capacity * targetState);
        for (unsigned batches : {1u, 2u, 17u}) {
            std::vector<Complex> input(capacity * sourceState, {nan, nan}), expected;
            for (unsigned batch = 0; batch < batches; ++batch) {
                std::vector<Complex> one(sourceState);
                for (std::size_t j = 0; j < sourceState; ++j)
                    input[batch * sourceState + j] =
                        one[j] = {std::sin(.13 * j + batch), std::cos(.27 * j - .31 * batch)};
                const auto reference = applyScaledAxialPlan(plan, one);
                expected.insert(expected.end(), reference.begin(), reference.end());
            }
            source.upload(input, axial.stream());
            output.upload(std::vector<Complex>(capacity * targetState, sentinel), axial.stream());
            status.reset(axial.stream());
            axial.enqueueDeviceBatch(source.values, output.values, batches, status.value);
            require(status.finish(axial.stream()) == 0,
                    "Cooperative axial translation read poisoned inactive batches.");
            check(output.download(), expected);
        }
    }
    std::cout << "Real CUDA cooperative FP64 contractions: irregular orders, 17 active batches, "
                 "compact prefixes and asymmetric axial blocks passed.\n";
}

void checkQueuedAzimuthUpdates() {
    constexpr unsigned sourceOrder = 17, nativeTargetOrder = 20, targetOrder = 11, capacity = 3;
    constexpr std::size_t sourceState = (sourceOrder + 1) * (sourceOrder + 1),
                          targetState = (targetOrder + 1) * (targetOrder + 1),
                          nativeTargetState = (nativeTargetOrder + 1) * (nativeTargetOrder + 1);
    auto sourcePlan = makeSphericalRotationPlan(sourceOrder, {7, -5, 3}),
         targetPlan = makeSphericalRotationPlan(nativeTargetOrder, {7, -5, 3});
    ScaledAxialPlan axialPlan;
    axialPlan.sourceOrder = sourceOrder;
    axialPlan.targetOrder = targetOrder;
    for (unsigned m = 0; m <= targetOrder; ++m) {
        axialPlan.blockOffsets.push_back(axialPlan.coefficients.size());
        for (unsigned l = m; l <= targetOrder; ++l)
            for (unsigned n = m; n <= sourceOrder; ++n)
                axialPlan.coefficients.emplace_back(std::cos(.13 * l + .21 * n + m) / (1 + n),
                                                    std::sin(.17 * l - .31 * n + m) / (1 + n));
    }
    axialPlan.blockOffsets.push_back(axialPlan.coefficients.size());
    SphericalRotationOperator sourceRotation(sourcePlan, capacity);
    const auto stream = sourceRotation.stream();
    SphericalRotationOperator targetRotation(targetPlan, capacity, {}, stream);
    ScaledAxialOperator axial(axialPlan, capacity, {}, stream);
    const auto sourceIdentity = sourceRotation.deviceIdentity(),
               targetIdentity = targetRotation.deviceIdentity();
    TestDeviceArray source(capacity * sourceState), invalidSource(capacity * sourceState),
        aligned(capacity * sourceState), translated(capacity * targetState),
        first(capacity * nativeTargetState), middle(capacity * nativeTargetState),
        last(capacity * nativeTargetState);
    std::vector<Complex> input(capacity * sourceState);
    for (std::size_t j = 0; j < input.size(); ++j)
        input[j] = {std::sin(.13 * j), std::cos(.27 * j)};
    source.upload(input, stream);
    auto invalid = input;
    invalid[0] = {std::numeric_limits<double>::quiet_NaN(), 0};
    invalidSource.upload(invalid, stream);
    const Complex sentinel{-.41, .73};
    TestDeviceStatus status;
    const auto enqueue = [&](const cuDoubleComplex *values, cuDoubleComplex *output,
                             unsigned batches) {
        sourceRotation.enqueueToAxisDeviceBatch(values, aligned.values, batches, status.value);
        axial.enqueueDeviceBatch(aligned.values, translated.values, batches, status.value);
        targetRotation.enqueueFromAxisPrefixDeviceBatch(translated.values, output, targetOrder,
                                                        batches, status.value);
    };
    const auto phase = [&](double phi) {
        sourceRotation.setAzimuth(phi);
        targetRotation.setAzimuth(phi);
    };
    detail::M2lChunkGraph graph;
    detail::M2lChunkGraph::Preparation timing;
    require(graph.prepare(
                stream, 128, [&] { enqueue(source.values, middle.values, 2); }, timing),
            "Queued azimuth gate could not capture its native/prefix RCR chain.");
    const auto check = [&](const TestDeviceArray &output, unsigned batches, double phi) {
        for (int m = -static_cast<int>(sourceOrder); m <= static_cast<int>(sourceOrder); ++m)
            sourcePlan.azimuthPhases[sourceOrder + m] = wave(m * phi);
        for (int m = -static_cast<int>(nativeTargetOrder); m <= static_cast<int>(nativeTargetOrder);
             ++m)
            targetPlan.azimuthPhases[nativeTargetOrder + m] = wave(m * phi);
        const auto actual = output.download();
        double error = 0, norm = 0;
        for (unsigned batch = 0; batch < batches; ++batch) {
            const std::span<const Complex> one(input.data() + batch * sourceState, sourceState);
            const auto local = applyScaledAxialPlan(axialPlan, rotateToAxis(sourcePlan, one));
            std::vector<Complex> full(nativeTargetState);
            std::copy(local.begin(), local.end(), full.begin());
            const auto reference = rotateFromAxis(targetPlan, full);
            for (std::size_t j = 0; j < targetState; ++j) {
                error += std::norm(actual[batch * targetState + j] - reference[j]);
                norm += std::norm(reference[j]);
            }
        }
        require(std::sqrt(error / norm) < 2e-12,
                "Queued phase update changed an earlier direction or a compact target prefix.");
        for (std::size_t j = batches * targetState; j < actual.size(); ++j)
            require(actual[j] == sentinel, "Queued phase/RCR chain overwrote inactive storage.");
    };
    for (const bool poison : {false, true, false}) {
        first.upload(std::vector<Complex>(capacity * nativeTargetState, sentinel), stream);
        middle.upload(std::vector<Complex>(capacity * nativeTargetState, sentinel), stream);
        last.upload(std::vector<Complex>(capacity * nativeTargetState, sentinel), stream);
        status.reset(stream);
        phase(.37);
        enqueue(poison ? invalidSource.values : source.values, first.values, 1);
        phase(-1.12);
        graph.launch(stream);
        phase(2.04);
        enqueue(source.values, last.values, 3);
        // No upload, download or synchronization separates these directions.
        require((status.finish(stream) != 0) == poison,
                "Canonical/graph direction changes lost a sticky error or caller epoch reset.");
        if (!poison)
            check(first, 1, .37);
        check(middle, 2, -1.12);
        check(last, 3, 2.04);
    }
    // A captured setter would fail immediately if it waited on its stream;
    // capturing it also checks that no temporary host phase payload is used.
    detail::M2lChunkGraph phaseGraph;
    detail::M2lChunkGraph::Preparation phaseTiming;
    require(phaseGraph.prepare(
                stream, 128,
                [&] {
                    phase(-.66);
                    enqueue(source.values, first.values, 1);
                },
                phaseTiming),
            "Device azimuth update is not capturable on its borrowed stream.");
    status.reset(stream);
    phaseGraph.launch(stream);
    require(status.finish(stream) == 0, "Captured device phase update failed.");
    check(first, 1, -.66);
    for (const double phi :
         {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity(),
          std::numeric_limits<double>::max()}) {
        bool rejected = false;
        try {
            sourceRotation.setAzimuth(phi);
        } catch (const std::invalid_argument &) {
            rejected = true;
        }
        require(rejected, "Azimuth setter enqueued a non-finite angle or scaled phase.");
    }
    const auto reusedSource = sourceRotation.deviceIdentity(),
               reusedTarget = targetRotation.deviceIdentity();
    require(reusedSource.phases == sourceIdentity.phases &&
                reusedSource.matrix == sourceIdentity.matrix &&
                reusedTarget.phases == targetIdentity.phases &&
                reusedTarget.matrix == targetIdentity.matrix,
            "Queued azimuth changes replaced captured matrix/phase identities.");
    std::cout << "CUDA queued azimuth canonical/graph/canonical order, sticky epoch, "
                 "compact prefix and phase capture passed.\n";
}

void checkRotationPrefixes() {
    constexpr unsigned capacity = 3;
    const Complex sentinel{-.27, .63};
    const double nan = std::numeric_limits<double>::quiet_NaN();
    SphericalRotationOperator streamOwner(makeSphericalRotationPlan(0, {0, 0, 1}), 1);
    for (const bool forward : {true, false})
        for (const Vec3 axis : {Vec3{7, -5, 3}, Vec3{1e-10, -2e-10, 1}, Vec3{1e-10, 2e-10, -1}}) {
            const unsigned planOrder = axis.z == 3 ? 128 : 8;
            auto plan = makeSphericalRotationPlan(planOrder, axis);
            const auto fullState = static_cast<std::size_t>(planOrder + 1) * (planOrder + 1);
            SphericalRotationOperator rotation(plan, capacity, {}, streamOwner.stream());
            const auto identity = rotation.deviceIdentity();
            const auto bytes = rotation.deviceBytes();
            TestDeviceArray source(capacity * fullState), output(capacity * fullState);
            TestDeviceStatus status;
            const auto applyPrefix = [&](unsigned order, unsigned batches) {
                if (forward)
                    rotation.toAxisPrefixDeviceBatch(source.values, output.values, order, batches);
                else
                    rotation.fromAxisPrefixDeviceBatch(source.values, output.values, order,
                                                       batches);
            };
            const auto enqueuePrefix = [&](unsigned order, unsigned batches, int *failure,
                                           std::stop_token stop = {}) {
                if (forward)
                    rotation.enqueueToAxisPrefixDeviceBatch(source.values, output.values, order,
                                                            batches, failure, stop);
                else
                    rotation.enqueueFromAxisPrefixDeviceBatch(source.values, output.values, order,
                                                              batches, failure, stop);
            };
            const std::vector<unsigned> orders =
                planOrder == 128
                    ? std::vector<unsigned>{0, 1, 31, 32, 63, 64, 95, 96, 127, 128, 1, 0}
                    : std::vector<unsigned>{0, 1, 7, 8, 1, 0};
            for (std::size_t caseId = 0; caseId < orders.size(); ++caseId) {
                const auto order = orders[caseId];
                const auto state = static_cast<std::size_t>(order + 1) * (order + 1);
                const unsigned batches = capacity - static_cast<unsigned>(caseId % capacity);
                const bool zero = caseId + 1 == orders.size(), queued = caseId % 2 == 0;
                // Changing phi must retain the native plan's phase center, including
                // both signs of m when the compact order differs from planOrder.
                const double phi = .37 - .19 * static_cast<double>(caseId);
                rotation.setAzimuth(phi);
                for (int m = -static_cast<int>(planOrder); m <= static_cast<int>(planOrder); ++m)
                    plan.azimuthPhases[planOrder + m] = wave(m * phi);
                std::vector<Complex> input(capacity * fullState, {nan, nan});
                std::vector<Complex> expected;
                for (unsigned batch = 0; batch < batches; ++batch) {
                    std::vector<Complex> fullInput(fullState);
                    for (unsigned l = 0; l <= order; ++l)
                        for (int m = -static_cast<int>(l); m <= static_cast<int>(l); ++m) {
                            const auto j = static_cast<std::size_t>(l) * l + l + m;
                            const Complex value = zero ? Complex{}
                                                       : Complex{std::cos(.3 * m + batch + caseId),
                                                                 std::sin(.7 * m - batch - .1 * l)};
                            input[batch * state + j] = value;
                            fullInput[j] = value;
                        }
                    const auto reference =
                        forward ? rotateToAxis(plan, fullInput) : rotateFromAxis(plan, fullInput);
                    expected.insert(expected.end(), reference.begin(), reference.begin() + state);
                }
                source.upload(input, rotation.stream());
                output.upload(std::vector<Complex>(capacity * fullState, sentinel),
                              rotation.stream());
                if (queued) {
                    status.reset(rotation.stream());
                    enqueuePrefix(order, batches, status.value);
                    require(status.finish(rotation.stream()) == 0,
                            "Compact rotation read poisoned trailing scratch.");
                } else
                    applyPrefix(order, batches);
                const auto actual = output.download();
                double error = 0, norm = 0;
                for (std::size_t j = 0; j < expected.size(); ++j) {
                    error += std::norm(actual[j] - expected[j]);
                    norm += std::norm(expected[j]);
                }
                require(norm == 0 ? error == 0 : std::sqrt(error / norm) < 2e-12,
                        "Compact rotation differs from the zero-extended native CPU rotation.");
                for (std::size_t j = expected.size(); j < actual.size(); ++j)
                    require(actual[j] == sentinel,
                            "Compact rotation overwrote trailing reserved storage.");
                const auto reused = rotation.deviceIdentity();
                require(
                    rotation.deviceBytes() == bytes && reused.offsets == identity.offsets &&
                        reused.matrix == identity.matrix && reused.phases == identity.phases &&
                        reused.order == planOrder && reused.capacity == identity.capacity,
                    "Rotation prefix allocated a separate plan or changed its native identity.");
            }
            // A graph reads the native phase allocation at replay time. Prefix
            // strides remain compact while the stored plan and phase center stay
            // native; later azimuth changes must not require another allocation.
            const unsigned graphOrder = std::min(17u, planOrder), graphBatches = 2;
            const auto graphState = static_cast<std::size_t>(graphOrder + 1) * (graphOrder + 1);
            std::vector<Complex> graphInput(capacity * fullState, {nan, nan});
            for (std::size_t j = 0; j < graphBatches * graphState; ++j)
                graphInput[j] = {std::sin(.13 * j + .4), std::cos(.19 * j - .7)};
            source.upload(graphInput, rotation.stream());
            output.upload(std::vector<Complex>(capacity * fullState, sentinel), rotation.stream());
            detail::M2lChunkGraph graph;
            detail::M2lChunkGraph::Preparation timing;
            require(graph.prepare(
                        rotation.stream(), 128,
                        [&] { enqueuePrefix(graphOrder, graphBatches, status.value); }, timing),
                    "Compact native rotation could not capture its actual kernel chain.");
            status.reset(rotation.stream());
            rotation.setAzimuth(.91);
            graph.launch(rotation.stream());
            rotation.setAzimuth(-.73);
            graph.launch(rotation.stream());
            require(status.finish(rotation.stream()) == 0,
                    "Compact graph rotation read poisoned trailing storage.");
            for (int m = -static_cast<int>(planOrder); m <= static_cast<int>(planOrder); ++m)
                plan.azimuthPhases[planOrder + m] = wave(-.73 * m);
            const auto graphActual = output.download();
            double graphError = 0, graphNorm = 0;
            for (unsigned batch = 0; batch < graphBatches; ++batch) {
                std::vector<Complex> fullInput(fullState);
                std::copy_n(graphInput.begin() + batch * graphState, graphState, fullInput.begin());
                const auto reference =
                    forward ? rotateToAxis(plan, fullInput) : rotateFromAxis(plan, fullInput);
                for (std::size_t j = 0; j < graphState; ++j) {
                    graphError += std::norm(graphActual[batch * graphState + j] - reference[j]);
                    graphNorm += std::norm(reference[j]);
                }
            }
            require(
                std::sqrt(graphError / graphNorm) < 2e-12,
                "Compact graph rotation differs from native zero-extension after phase update.");
            for (std::size_t j = graphBatches * graphState; j < graphActual.size(); ++j)
                require(graphActual[j] == sentinel,
                        "Compact graph rotation modified trailing reserved storage.");
            graphInput[0] = {nan, nan};
            source.upload(graphInput, rotation.stream());
            status.reset(rotation.stream());
            graph.launch(rotation.stream());
            graphInput[0] = {};
            source.upload(graphInput, rotation.stream());
            graph.launch(rotation.stream());
            require(status.finish(rotation.stream()) != 0,
                    "Compact graph replay erased an earlier non-finite prefix.");
            status.reset(rotation.stream());
            graph.launch(rotation.stream());
            require(status.finish(rotation.stream()) == 0,
                    "Compact graph replay retained an earlier explicitly reset failure.");
            const auto replayIdentity = rotation.deviceIdentity();
            require(rotation.deviceBytes() == bytes && replayIdentity.offsets == identity.offsets &&
                        replayIdentity.matrix == identity.matrix &&
                        replayIdentity.phases == identity.phases,
                    "Compact graph replay replaced its immutable native plan storage.");

            const auto reject = [&](unsigned order, unsigned batches) {
                bool rejected = false;
                try {
                    applyPrefix(order, batches);
                } catch (const std::invalid_argument &) {
                    rejected = true;
                }
                require(rejected, "Compact rotation accepted invalid order/batch dimensions.");
            };
            reject(static_cast<unsigned>(-1), 1);
            reject(planOrder + 1, 1);
            reject(0, 0);
            reject(0, capacity + 1);
            for (const auto pointers :
                 {std::pair<const cuDoubleComplex *, cuDoubleComplex *>{nullptr, output.values},
                  std::pair<const cuDoubleComplex *, cuDoubleComplex *>{source.values, nullptr},
                  std::pair<const cuDoubleComplex *, cuDoubleComplex *>{source.values,
                                                                        source.values}}) {
                bool rejected = false;
                try {
                    if (forward)
                        rotation.toAxisPrefixDeviceBatch(pointers.first, pointers.second, 0, 1);
                    else
                        rotation.fromAxisPrefixDeviceBatch(pointers.first, pointers.second, 0, 1);
                } catch (const std::invalid_argument &) {
                    rejected = true;
                }
                require(rejected, "Compact rotation accepted null or aliased device arrays.");
            }
            bool rejected = false;
            try {
                enqueuePrefix(0, 1, nullptr);
            } catch (const std::invalid_argument &) {
                rejected = true;
            }
            require(rejected, "Queued compact rotation accepted a missing failure status.");
            std::stop_source stop;
            stop.request_stop();
            bool cancelled = false;
            try {
                enqueuePrefix(0, 1, status.value, stop.get_token());
            } catch (const std::runtime_error &error) {
                cancelled =
                    std::string_view(error.what()).find("cancelled") != std::string_view::npos;
            }
            require(cancelled, "Compact rotation ignored cancellation.");
            // An early non-finite prefix must remain sticky after a valid overwrite.
            std::vector<Complex> input(capacity * fullState, {nan, nan});
            source.upload(input, rotation.stream());
            status.reset(rotation.stream());
            enqueuePrefix(0, 1, status.value);
            input[0] = {};
            source.upload(input, rotation.stream());
            enqueuePrefix(0, 1, status.value);
            require(status.finish(rotation.stream()) != 0,
                    "Compact rotation erased an earlier non-finite failure.");
            status.reset(rotation.stream());
            enqueuePrefix(0, 1, status.value);
            require(status.finish(rotation.stream()) == 0 && output.download()[0] == Complex{},
                    "Compact rotation could not reuse a reset status after cancellation.");
        }
    require(cudaStreamQuery(streamOwner.stream()) == cudaSuccess,
            "Destroying a prefix rotation invalidated its borrowed stream.");
    std::cout
        << "Real CUDA native-plan compact forward/reverse prefixes/batches/zero/cancel pass\n";
}

void checkSphericalRotations() {
    const Vec3 sourceDirection{.3, .4, std::sqrt(.75)};
    struct Case {
        unsigned order;
        Vec3 axis;
    };
    std::vector<Case> cases;
    for (const unsigned order : {0u, 8u, 19u, 20u, 21u, 24u, 40u, 64u, 128u})
        cases.push_back({order, {.37, -.51, -.78}});
    for (const unsigned order : {36u, 52u})
        for (const Vec3 axis : {Vec3{1e-10, -2e-10, 1}, Vec3{1, 0, 0}, Vec3{1e-10, 2e-10, -1},
                                Vec3{7, -5, 3}, Vec3{-2, 3, -1}})
            cases.push_back({order, axis});
    cases.push_back({128, {1, 0, 0}});
    for (const auto spec : cases) {
        const unsigned order = spec.order;
        const auto axis = spec.axis;
        const auto plan = makeSphericalRotationPlan(order, axis);
        const double length = std::hypot(axis.x, axis.y, axis.z), phi = std::atan2(axis.y, axis.x),
                     x = std::cos(phi) * sourceDirection.x + std::sin(phi) * sourceDirection.y,
                     y = -std::sin(phi) * sourceDirection.x + std::cos(phi) * sourceDirection.y;
        const double cosine = axis.z / length, sine = std::hypot(axis.x, axis.y) / length;
        const Vec3 rotated{cosine * x - sine * sourceDirection.z, y,
                           sine * x + cosine * sourceDirection.z};
        const auto original = sphericalHarmonics(order, sourceDirection),
                   expectedHarmonics = sphericalHarmonics(order, rotated);
        constexpr unsigned batches = 2;
        std::vector<Complex> input, expected;
        for (unsigned batch = 0; batch < batches; ++batch) {
            const Complex amplitude{.7 + batch, .3 - .4 * batch};
            for (std::size_t j = 0; j < original.size(); ++j) {
                input.push_back(amplitude * std::conj(original[j]));
                expected.push_back(amplitude * std::conj(expectedHarmonics[j]));
            }
        }
        SphericalRotationOperator gpu(plan, batches);
        TestDeviceArray source(input.size()), aligned(input.size()), restored(input.size());
        const auto classBytes = gpu.deviceBytes();
        source.upload(input, gpu.stream());
        gpu.toAxisDevice(source.values, aligned.values);
        gpu.fromAxisDevice(aligned.values, restored.values);
        const auto actual = aligned.download(), roundTrip = restored.download();
        double error = 0, norm = 0, returnError = 0;
        for (std::size_t j = 0; j < input.size(); ++j) {
            error += std::norm(actual[j] - expected[j]);
            norm += std::norm(expected[j]);
            returnError += std::norm(roundTrip[j] - input[j]);
        }
        const double relative = std::sqrt(error / norm), inverse = std::sqrt(returnError / norm);
        std::cout << "Real CUDA Wigner rotation p=" << order << ", axis=" << axis.x << ',' << axis.y
                  << ',' << axis.z << ", SH relative=" << relative << ", inverse=" << inverse
                  << ", bytes=" << gpu.deviceBytes() << '\n';
        require(relative < 2e-12,
                "GPU Wigner rotation fails independently evaluated spherical harmonics.");
        require(inverse < 2e-12, "GPU Wigner rotation is not unitary at its supported order.");
        const Complex reuseScale{-.31, .47};
        for (auto &value : input)
            value *= reuseScale;
        source.upload(input, gpu.stream());
        gpu.toAxisDevice(source.values, aligned.values);
        gpu.fromAxisDevice(aligned.values, restored.values);
        const auto reused = restored.download();
        double reuseError = 0, reuseNorm = 0;
        for (std::size_t j = 0; j < input.size(); ++j) {
            reuseError += std::norm(reused[j] - input[j]);
            reuseNorm += std::norm(input[j]);
        }
        require(std::sqrt(reuseError / reuseNorm) < 2e-12,
                "GPU Wigner nondefault-stream reuse retained stale input.");
        require(gpu.deviceBytes() == classBytes,
                "Device rotation retained per-class state/scratch storage.");
    }
}
void checkActiveLfBatches() {
    checkRotationPrefixes();
    constexpr unsigned order = 36, capacity = 3;
    const auto state = static_cast<std::size_t>(order + 1) * (order + 1);
    const auto rotationPlan = makeSphericalRotationPlan(order, {7, -5, 3});
    const auto axialPlan =
        makeScaledAxialPlan(ScaledAxialKind::multipoleToLocal, order, order, 10, .01, .01, .2, 160);
    SphericalRotationOperator rotation(rotationPlan, capacity);
    ScaledAxialOperator axial(axialPlan, capacity, {}, rotation.stream());
    TestDeviceArray source(capacity * state), aligned(capacity * state),
        translated(capacity * state), output(capacity * state);
    const Complex sentinel{-.27, .63};
    const std::vector<Complex> untouched(capacity * state, sentinel);
    TestDeviceStatus status;
    for (const bool queued : {false, true})
        for (unsigned active : {3u, 1u, 2u, 3u, 1u}) {
            std::vector<Complex> input(capacity * state);
            for (unsigned batch = 0; batch < capacity; ++batch)
                for (unsigned l = 0; l <= order; ++l)
                    for (int m = -static_cast<int>(l); m <= static_cast<int>(l); ++m)
                        input[batch * state + static_cast<std::size_t>(l) * l + l + m] =
                            active == 1
                                ? Complex{}
                                : std::pow(.1, l) * Complex{std::cos(.3 * m + batch + active),
                                                            std::sin(.7 * m - batch)};
            source.upload(input, rotation.stream());
            aligned.upload(untouched, rotation.stream());
            translated.upload(untouched, rotation.stream());
            output.upload(untouched, rotation.stream());
            if (queued) {
                status.reset(rotation.stream());
                rotation.enqueueToAxisDeviceBatch(source.values, aligned.values, active,
                                                  status.value);
                axial.enqueueDeviceBatch(aligned.values, translated.values, active, status.value);
                rotation.enqueueFromAxisDeviceBatch(translated.values, output.values, active,
                                                    status.value);
                require(status.finish(rotation.stream()) == 0,
                        "Queued finite RCR coefficients reported an error.");
            } else {
                rotation.toAxisDeviceBatch(source.values, aligned.values, active);
                axial.applyDeviceBatch(aligned.values, translated.values, active);
                rotation.fromAxisDeviceBatch(translated.values, output.values, active);
            }
            const auto a = aligned.download(), t = translated.download(),
                       actual = output.download();
            double error = 0, norm = 0;
            for (unsigned batch = 0; batch < active; ++batch) {
                const std::span<const Complex> one(input.data() + batch * state, state);
                const auto reference = rotateFromAxis(
                    rotationPlan, applyScaledAxialPlan(axialPlan, rotateToAxis(rotationPlan, one)));
                for (std::size_t j = 0; j < state; ++j) {
                    error += std::norm(actual[batch * state + j] - reference[j]);
                    norm += std::norm(reference[j]);
                }
            }
            require(norm == 0 ? error == 0 : std::sqrt(error / norm) < 2e-12,
                    "Active LF batch RCR fails its independent CPU composition.");
            for (std::size_t j = active * state; j < capacity * state; ++j)
                require(a[j] == sentinel && t[j] == sentinel && actual[j] == sentinel,
                        "Active LF batch modified trailing reserved storage.");
        }
    // An invalid intermediate must survive later finite overwrites. This
    // tests accumulation rather than merely rechecking the final output.
    std::vector<Complex> invalid(capacity * state);
    invalid[0] = {std::numeric_limits<double>::quiet_NaN(), 0};
    TestDeviceArray zeros(capacity * state);
    zeros.upload(std::vector<Complex>(capacity * state), rotation.stream());
    source.upload(invalid, rotation.stream());
    status.reset(rotation.stream());
    // No intervening upload or synchronization: a first invalid pipeline is
    // overwritten by 31 valid pipelines before the single completion check.
    for (unsigned batch = 0; batch < 32; ++batch) {
        rotation.enqueueToAxisDeviceBatch(batch ? zeros.values : source.values, aligned.values, 1,
                                          status.value);
        axial.enqueueDeviceBatch(aligned.values, translated.values, 1, status.value);
        rotation.enqueueFromAxisDeviceBatch(translated.values, output.values, 1, status.value);
    }
    require(status.finish(rotation.stream()) != 0,
            "Queued RCR erased a non-finite intermediate failure.");
    status.reset(rotation.stream());
    for (unsigned batch = 0; batch < 32; ++batch) {
        rotation.enqueueToAxisDeviceBatch(zeros.values, aligned.values, 1, status.value);
        axial.enqueueDeviceBatch(aligned.values, translated.values, 1, status.value);
        rotation.enqueueFromAxisDeviceBatch(translated.values, output.values, 1, status.value);
    }
    require(status.finish(rotation.stream()) == 0,
            "Queued RCR retained a failure after explicit status reset and valid reuse.");
    const auto reused = output.download();
    for (std::size_t j = 0; j < state; ++j)
        require(reused[j] == Complex{}, "Queued RCR zero reuse retained an earlier pipeline.");
    bool rejected = false;
    try {
        rotation.toAxisDeviceBatch(source.values, aligned.values, capacity + 1);
    } catch (const std::invalid_argument &) {
        rejected = true;
    }
    require(rejected, "Rotation accepted an active batch larger than capacity.");
    rejected = false;
    try {
        rotation.toAxisDeviceBatch(source.values, source.values, 1);
    } catch (const std::invalid_argument &) {
        rejected = true;
    }
    require(rejected, "Rotation accepted an output aliasing its source.");
    rejected = false;
    try {
        axial.applyDeviceBatch(source.values, translated.values, capacity + 1);
    } catch (const std::invalid_argument &) {
        rejected = true;
    }
    require(rejected, "Axial translation accepted an active batch larger than capacity.");
    std::cout << "Real CUDA active LF batches full/tail/zero/reuse pass\n";
}

void checkGraphStickyReplay() {
    constexpr unsigned order = 12;
    const auto state = std::size_t(order + 1) * (order + 1);
    SphericalRotationOperator rotation(makeSphericalRotationPlan(order, {7, -5, 3}), 1);
    ScaledAxialOperator axial(
        makeScaledAxialPlan(ScaledAxialKind::multipoleToLocal, order, order, 10, .01, .01, .2, 96),
        1, {}, rotation.stream());
    TestDeviceArray source(state), zero(state), aligned(state), translated(state), output(state);
    TestDeviceStatus status;
    source.upload(std::vector<Complex>(state), rotation.stream());
    zero.upload(std::vector<Complex>(state), rotation.stream());
    const auto chain = [&] {
        for (unsigned i = 0; i < 32; ++i) {
            rotation.enqueueToAxisDeviceBatch(i ? zero.values : source.values, aligned.values, 1,
                                              status.value);
            axial.enqueueDeviceBatch(aligned.values, translated.values, 1, status.value);
            rotation.enqueueFromAxisDeviceBatch(translated.values, output.values, 1, status.value);
        }
    };
    detail::M2lChunkGraph graph;
    detail::M2lChunkGraph::Preparation timing;
    require(graph.prepare(rotation.stream(), 4096, chain, timing) && timing.nodes > 1,
            "Actual RCR graph did not instantiate its real captured nodes.");
    status.reset(rotation.stream());
    graph.launch(rotation.stream());
    require(status.finish(rotation.stream()) == 0, "Finite RCR graph replay failed.");
    auto invalid = std::vector<Complex>(state);
    invalid[0] = {std::numeric_limits<double>::quiet_NaN(), 0};
    source.upload(invalid, rotation.stream());
    status.reset(rotation.stream());
    graph.launch(rotation.stream());
    require(status.finish(rotation.stream()) != 0,
            "Graph replay erased an early NaN after 31 finite overwrites.");
    for (auto value : output.download())
        require(value == Complex{}, "Late graph zero pipeline did not overwrite its output.");
    source.upload(std::vector<Complex>(state), rotation.stream());
    status.reset(rotation.stream());
    graph.launch(rotation.stream());
    require(status.finish(rotation.stream()) == 0,
            "Caller epoch reset did not permit valid graph reuse.");
    detail::M2lChunkGraph denied;
    detail::M2lChunkGraph::Preparation deniedTiming;
    require(!denied.prepare(rotation.stream(), 1, chain, deniedTiming) && !denied &&
                deniedTiming.nodes == timing.nodes,
            "Graph node cap did not reject the actual recorded node count.");
    status.reset(rotation.stream());
    chain();
    require(status.finish(rotation.stream()) == 0,
            "Canonical execution after graph node-cap rejection failed.");
    for (auto value : output.download())
        require(value == Complex{}, "Canonical graph-budget path changed numerical work.");
    detail::M2lChunkGraph cancelledGraph;
    detail::M2lChunkGraph::Preparation cancelledTiming;
    bool stopped = false;
    try {
        cancelledGraph.prepare(
            rotation.stream(), 4096,
            [&] {
                chain();
                throw std::runtime_error("Graph capture cancellation probe");
            },
            cancelledTiming);
    } catch (const std::runtime_error &error) {
        stopped =
            std::string_view(error.what()).find("cancellation probe") != std::string_view::npos;
    }
    cudaStreamCaptureStatus captureStatus{};
    require(stopped && cudaStreamIsCapturing(rotation.stream(), &captureStatus) == cudaSuccess &&
                captureStatus == cudaStreamCaptureStatusNone,
            "Capture exception left the borrowed stream capturing.");
    status.reset(rotation.stream());
    graph.launch(rotation.stream());
    require(status.finish(rotation.stream()) == 0,
            "Capture exception corrupted an earlier valid graph executable.");
    std::cout << "CUDA graph sticky/reset/node-budget/capture-cleanup replay pass, nodes="
              << timing.nodes << '\n';
}

void checkM2lGraphLifecycle() {
    auto fixture = test::makeMixedWidebandFixture(20, 4);
    const auto size = fixture.space.vertices.size();
    for (auto &panel : fixture.space.panels)
        panel.boundary = BoundaryCondition::driven;
    SupportTreeSettings treeSettings;
    treeSettings.leafCapacity = 2;
    treeSettings.nearDistanceRatio = .5;
    auto tree = std::make_shared<SupportTree>(makeSupportTree(fixture.space, treeSettings));
    const std::array<std::uint32_t, 2> faces{5, 1};
    auto plan = std::make_shared<SupportNearPlan>(
        makeSupportNearPlan(fixture.space, *tree, faces, 16 * 1024 * 1024));
    SupportNearOperatorSettings nearSettings;
    nearSettings.quadrature.nearDistanceRatio = .5;
    auto near =
        std::make_shared<SupportNearOperator>(fixture.space, *tree, *plan, 20, nearSettings);
    auto pc = std::make_shared<SupportBlockJacobi>(near->pressureMatrix(), tree->pressureOrder, 32,
                                                   std::stop_token{}, near->stream());
    SupportFmmOperatorSettings settings;
    settings.role = SupportFmmOperatorRole::distancePrimary;
    settings.panelRule = fixture.triangleRule;
    auto canonicalSettings = settings;
    canonicalSettings.maximumM2lGraphEntries = 0;
    SupportFmmOperator canonical(fixture.space, tree, plan, near, pc, 20, canonicalSettings);
    SupportFmmOperator graph(fixture.space, tree, plan, near, pc, 20, settings);
    TestDeviceArray pressure(size), flux(2), direct(size), captured(size);
    for (unsigned repeat = 0; repeat < 4; ++repeat) {
        std::vector<Complex> p;
        for (std::size_t i = 0; i < size; ++i)
            p.emplace_back(std::cos(.7 * i + repeat), std::sin(.9 * i - .3 * repeat));
        pressure.upload(p, graph.stream());
        flux.upload({{.3 + repeat, -.7}, {-.4, .2 + .3 * repeat}}, graph.stream());
        for (unsigned role = 0; role < 2; ++role) {
            if (role) {
                canonical.applyFaceFluxDevice(flux.values, direct.values);
                graph.applyFaceFluxDevice(flux.values, captured.values);
            } else {
                canonical.applyPressureDevice(pressure.values, direct.values);
                graph.applyPressureDevice(pressure.values, captured.values);
            }
            const auto expected = direct.download(), actual = captured.download();
            require(expected == actual,
                    "Graph/canonical M2L paths changed the same ordered FP64 kernel chain.");
        }
        if (repeat == 0)
            require(graph.workStatistics().m2lGraphCaptureAttempts == 0,
                    "Single-contact Ap/Bq paid graph recording costs or shared a role key.");
    }
    const auto work = graph.workStatistics();
    require(work.m2lGraphLaunches && work.m2lGraphReuses &&
                work.m2lGraphInstantiations == work.m2lGraphCachedExecutables &&
                work.m2lGraphCachedNodes <= settings.maximumM2lGraphNodes &&
                work.m2lGraphOwnedHostPayloadBytes <= settings.maximumM2lGraphHostBytes,
            "Repeated Ap/Bq did not use bounded, genuinely instantiated graph executables.");
    require(graph.hasDistanceTargetGraphs(),
            "Graph lifecycle gate did not capture/execute a genuinely shortened endpoint class.");
    require(canonical.workStatistics().m2lGraphLaunches == 0 &&
                canonical.workStatistics().m2lGraphDirectChunks &&
                canonical.workStatistics().m2lGraphAdmissionClosures == 1,
            "Zero metadata budget did not retain canonical execution and close admission once.");
    const auto launches = work.m2lGraphLaunches;
    graph.releaseDeviceState();
    require(!graph.hasDistanceTargetGraphs(),
            "Hierarchy release retained an endpoint-prefix executable.");
    require(graph.workStatistics().m2lGraphCachedExecutables == 0 &&
                graph.workStatistics().m2lGraphOwnedHostPayloadBytes == 0 &&
                graph.workStatistics().m2lGraphInvalidatedExecutables >=
                    work.m2lGraphCachedExecutables,
            "State release retained captured matrix/state references.");
    pressure.upload(std::vector<Complex>(size), graph.stream());
    graph.applyPressureDevice(pressure.values, captured.values);
    require(graph.role() == SupportFmmOperatorRole::distancePrimary,
            "Rebuilt graph hierarchy lost its primary class-order policy.");
    require(graph.workStatistics().m2lGraphLaunches == launches,
            "Rebuilt hierarchy reused a graph from an earlier state generation.");
    for (const auto value : captured.download())
        require(value == Complex{}, "Graph hierarchy rebuild retained an earlier RHS.");
    std::cout << "CUDA M2L graph Ap/Bq exact parity/reuse/budget/release pass, launches="
              << launches << ", nodes=" << work.m2lGraphCachedNodes << '\n';
}

void checkRcrQueuePolicy() {
    require(makeRcrQueuePolicy(0, 0, 32).workPerLink == 5,
            "RCR work bound omits an endpoint phase for constant modes.");
    for (const unsigned source : {0u, 1u, 19u, 24u, 40u, 52u, 80u, 128u})
        for (const unsigned target : {0u, 19u, 52u, 128u})
            for (const unsigned reserved : {1u, 3u, 32u, 256u}) {
                const auto policy = makeRcrQueuePolicy(source, target, reserved);
                require(policy.batchCapacity && policy.batchCapacity <= reserved &&
                            policy.maximumChunkBatches && policy.maximumChunkBatches <= 32,
                        "RCR queue policy exceeds reserved batches or chunk count.");
                require(policy.workPerLink * policy.batchCapacity * policy.maximumChunkBatches <=
                            RcrQueuePolicy::maximumChunkWork,
                        "RCR queue policy exceeds its conservative work cap.");
            }
    require(makeRcrQueuePolicy(19, 19, 32).maximumChunkBatches == 32 &&
                makeRcrQueuePolicy(80, 80, 32).batchCapacity < 32 &&
                makeRcrQueuePolicy(128, 128, 32).maximumChunkBatches == 1,
            "RCR queue does not shorten high-order chunks.");
    bool rejected = false;
    try {
        (void)makeRcrQueuePolicy(129, 20, 32);
    } catch (const std::invalid_argument &) {
        rejected = true;
    }
    require(rejected, "RCR queue accepted an unsupported LF order.");
}

void checkResidentRcr() {
    constexpr unsigned order = 10;
    constexpr double k = 10, sourceRadius = .01, targetRadius = .015, distance = .2;
    const Vec3 axis{.3, -.4, .5};
    const double axisNorm = std::hypot(axis.x, axis.y, axis.z);
    const Vec3 separation{distance * axis.x / axisNorm, distance * axis.y / axisNorm,
                          distance * axis.z / axisNorm},
        sourcePoint{.003, .001, -.002}, targetPoint{-.002, .004, .001};
    const Complex amplitude{1.1, -.6};
    const auto sourceBasis = scaledRegularWaves(order, k, sourceRadius, sourcePoint);
    std::vector<Complex> sourceCoefficients;
    for (const auto &wave : sourceBasis.waves)
        sourceCoefficients.push_back(amplitude * std::conj(wave.value));
    const auto rotationPlan = makeSphericalRotationPlan(order, separation);
    const auto axialPlan = makeScaledAxialPlan(ScaledAxialKind::multipoleToLocal, order, order, k,
                                               sourceRadius, targetRadius, distance, 48);
    SphericalRotationOperator rotation(rotationPlan, 1);
    ScaledAxialOperator axial(axialPlan, 1, {}, rotation.stream());
    const auto count = sourceCoefficients.size();
    TestDeviceArray original(count), aligned(count), localAligned(count), localGlobal(count);
    original.upload(sourceCoefficients, rotation.stream());
    rotation.toAxisDevice(original.values, aligned.values);
    axial.applyDevice(aligned.values, localAligned.values);
    rotation.fromAxisDevice(localAligned.values, localGlobal.values);
    const ScaledLocal local{order, targetRadius, localGlobal.download()};
    const auto actual = evaluateScaledLocal(local, k, targetPoint);
    const Vec3 delta{separation.x + targetPoint.x - sourcePoint.x,
                     separation.y + targetPoint.y - sourcePoint.y,
                     separation.z + targetPoint.z - sourcePoint.z};
    const auto exact = amplitude * green(delta, k);
    const double r = std::sqrt(dot(delta, delta));
    const auto radial = exact * Complex{-1 / r, -k};
    const double axes[]{delta.x, delta.y, delta.z};
    double error = 0, norm = 0;
    for (unsigned i = 0; i < 3; ++i) {
        const auto expected = radial * (axes[i] / r);
        error += std::norm(actual.gradient[i] - expected);
        norm += std::norm(expected);
    }
    const double field = std::abs(actual.value - exact) / std::abs(exact),
                 gradient = std::sqrt(error / norm);
    std::cout << "Real CUDA resident RCR field=" << field << ", gradient=" << gradient << '\n';
    require(field < 2e-9, "CUDA scaled RCR fails the independent off-axis Green field.");
    require(gradient < 2e-9, "CUDA scaled RCR fails the independent off-axis Green gradient.");
    // A high source multipole contributes to low local degrees. Retaining a
    // local prefix must still read the entire source column after rotation.
    for (const unsigned targetOrder : {0u, 1u, 4u, order - 1}) {
        const auto rectangular =
            makeScaledAxialPlan(ScaledAxialKind::multipoleToLocal, order, targetOrder, k,
                                sourceRadius, targetRadius, distance, axialPlan.polarCount);
        ScaledAxialOperator reduced(rectangular, 1, {}, rotation.stream());
        const auto prefix = std::size_t(targetOrder + 1) * (targetOrder + 1);
        TestDeviceArray compactAligned(prefix), compactGlobal(prefix);
        for (const bool zero : {false, true}) {
            std::vector<Complex> high(count);
            if (!zero) {
                high[std::size_t(order) * order] = {.7, -.4};
                high[std::size_t(order) * order + 2 * order] = {-.7, .4};
                high[std::size_t(order) * order + order] = {.2, .9};
            }
            original.upload(high, rotation.stream());
            rotation.toAxisDevice(original.values, aligned.values);
            axial.applyDevice(aligned.values, localAligned.values);
            rotation.fromAxisDevice(localAligned.values, localGlobal.values);
            reduced.applyDevice(aligned.values, compactAligned.values);
            rotation.fromAxisPrefixDeviceBatch(compactAligned.values, compactGlobal.values,
                                               targetOrder, 1);
            const auto full = localGlobal.download(), actualPrefix = compactGlobal.download();
            double difference = 0, magnitude = 0;
            for (std::size_t i = 0; i < prefix; ++i) {
                difference += std::norm(actualPrefix[i] - full[i]);
                magnitude += std::norm(full[i]);
            }
            require(zero ? difference == 0 && magnitude == 0
                         : magnitude > 0 && std::sqrt(difference / magnitude) < 3e-12,
                    "Rectangular RCR lost high source modes in retained local degrees.");
        }
    }
}

void checkSupportLeafPipeline() {
    for (const auto representation :
         {LevelRepresentation::regularSpherical, LevelRepresentation::planeWave}) {
        for (const auto cacheBudget :
             {std::size_t{0}, std::size_t{18000}, std::size_t{256 * 1024}}) {
            if (representation == LevelRepresentation::planeWave && cacheBudget)
                continue;
            const double k =
                representation == LevelRepresentation::regularSpherical ? 200 : 366.366;
            const auto fixture = test::makeWidebandFarBlockFixture(k);
            SupportTree tree;
            tree.nodes.resize(2);
            tree.nodes[0].center = fixture.sourceCenter;
            tree.nodes[1].center = fixture.targetCenter;
            for (unsigned i = 0; i < 8; ++i) {
                tree.pressureOrder.push_back(i);
                tree.fluxOrder.push_back(i);
                tree.pressureLeaf.push_back(i / 4);
                tree.fluxLeaf.push_back(i / 4);
            }
            for (unsigned i = 0; i < 2; ++i) {
                tree.nodes[i].firstPressure = tree.nodes[i].firstFlux = 4 * i;
                tree.nodes[i].pressureCount = tree.nodes[i].fluxCount = 4;
            }
            SupportLeafSettings sourceSettings;
            sourceSettings.nodes = {0};
            sourceSettings.maximumMomentCacheBytes = cacheBudget;
            sourceSettings.representation = representation;
            sourceSettings.order = 24;
            sourceSettings.wavenumber = k;
            sourceSettings.normalizationRadius = .03;
            sourceSettings.activeFaceCount = 8;
            sourceSettings.angularRule = makeSphereRule(64);
            sourceSettings.panelRule = fixture.triangleRule;
            for (unsigned i = 0; i < 8; ++i)
                sourceSettings.panelToActiveFaceColumn.push_back(i);
            auto targetSettings = sourceSettings;
            targetSettings.nodes = {1};
            cudaStream_t stream{};
            require(cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking) == cudaSuccess,
                    "Leaf test stream creation failed.");
            struct StreamGuard {
                cudaStream_t stream;
                ~StreamGuard() {
                    cudaStreamDestroy(stream);
                }
            } guard{stream};
            auto invalidRule = sourceSettings;
            for (auto &q : invalidRule.panelRule)
                q.weight = 0;
            bool rejected = false;
            try {
                SupportLeafOperator invalid(fixture.space, tree, invalidRule, {}, stream);
            } catch (const std::invalid_argument &) {
                rejected = true;
            }
            require(rejected, "Support leaf accepted a zero-measure triangle rule.");
            SupportLeafOperator sourceLeaf(fixture.space, tree, sourceSettings, {}, stream),
                targetLeaf(fixture.space, tree, targetSettings, {}, stream);
            const auto stateCount = sourceLeaf.stateSize();
            TestDeviceArray input(8), output(8), moment(stateCount), aligned(stateCount),
                localAligned(stateCount), local(stateCount);
            const Vec3 separation{fixture.targetCenter.x - fixture.sourceCenter.x,
                                  fixture.targetCenter.y - fixture.sourceCenter.y,
                                  fixture.targetCenter.z - fixture.sourceCenter.z};
            std::unique_ptr<SphericalRotationOperator> rotation;
            std::unique_ptr<ScaledAxialOperator> axial;
            std::unique_ptr<DiagonalM2LOperator> diagonal;
            if (representation == LevelRepresentation::regularSpherical) {
                rotation = std::make_unique<SphericalRotationOperator>(
                    makeSphericalRotationPlan(sourceSettings.order, separation), 1,
                    std::stop_token{}, stream);
                axial = std::make_unique<ScaledAxialOperator>(
                    makeScaledAxialPlan(ScaledAxialKind::multipoleToLocal, sourceSettings.order,
                                        sourceSettings.order, k, .03, .03,
                                        std::hypot(separation.x, separation.y, separation.z), 96),
                    1, std::stop_token{}, stream);
            } else {
                DiagonalLevelPlan plan;
                plan.sourceNodes = plan.targetNodes = plan.translationClasses = 1;
                plan.angularSamples = static_cast<unsigned>(stateCount);
                plan.targetOffsets = {0, 1};
                plan.interactions = {{0, 0}};
                plan.translationCoefficients =
                    diagonalTranslation(sourceSettings.angularRule, 32, k, separation);
                diagonal = std::make_unique<DiagonalM2LOperator>(plan, std::stop_token{}, stream);
            }
            for (unsigned rhs = 0; rhs < 2; ++rhs)
                for (unsigned channel = 0; channel < 2; ++channel) {
                    input.upload(channel ? fixture.flux[rhs] : fixture.pressure[rhs], stream);
                    if (channel)
                        sourceLeaf.emitFaceFluxDevice(input.values, moment.values);
                    else
                        sourceLeaf.emitPressureDevice(input.values, moment.values);
                    if (axial) {
                        rotation->toAxisDevice(moment.values, aligned.values);
                        axial->applyDevice(aligned.values, localAligned.values);
                        rotation->fromAxisDevice(localAligned.values, local.values);
                    } else
                        diagonal->applyDevice(moment.values, local.values);
                    const auto &expected =
                        channel ? fixture.expectedFlux[rhs] : fixture.expectedPressure[rhs];
                    const Complex nearValue = expected[4] * Complex{.17, -.09};
                    output.upload(std::vector<Complex>(8, nearValue), stream);
                    targetLeaf.testDevice(local.values, output.values);
                    const auto actual = output.download();
                    double error = 0, norm = 0;
                    for (unsigned i = 0; i < 8; ++i) {
                        error += std::norm(actual[i] - nearValue - expected[i]);
                        norm += std::norm(expected[i]);
                    }
                    const double relative = std::sqrt(error / norm);
                    std::cout << "CUDA complete-support leaf pipeline representation="
                              << static_cast<unsigned>(representation)
                              << ", cacheBudget=" << cacheBudget << ", channel=" << channel
                              << ", rhs=" << rhs << ", relative=" << relative << '\n';
                    require(
                        relative < 2e-9,
                        "CUDA complete-support leaf emission/RCR/testing failed independent BM.");
                }
            input.upload(std::vector<Complex>(8), stream);
            for (unsigned channel = 0; channel < 2; ++channel) {
                if (channel)
                    sourceLeaf.emitFaceFluxDevice(input.values, moment.values);
                else
                    sourceLeaf.emitPressureDevice(input.values, moment.values);
                if (axial) {
                    rotation->toAxisDevice(moment.values, aligned.values);
                    axial->applyDevice(aligned.values, localAligned.values);
                    rotation->fromAxisDevice(localAligned.values, local.values);
                } else
                    diagonal->applyDevice(moment.values, local.values);
                output.upload(std::vector<Complex>(8), stream);
                targetLeaf.testDevice(local.values, output.values);
                for (const auto value : output.download())
                    require(value == Complex{},
                            "Complete-support leaf retained a previous excitation.");
            }
        }
    }
}

void compareIntegratedLeafValues(const std::vector<Complex> &actual,
                                 const std::vector<Complex> &expected) {
    double error = 0, norm = 0;
    require(actual.size() == expected.size(), "Integrated leaf result size changed.");
    for (std::size_t i = 0; i < actual.size(); ++i) {
        require(std::isfinite(actual[i].real()) && std::isfinite(actual[i].imag()),
                "Integrated leaf produced a nonfinite value.");
        error += std::norm(actual[i] - expected[i]);
        norm += std::norm(expected[i]);
    }
    require(norm ? std::sqrt(error / norm) < 3e-12 : error == 0,
            "Integrated leaf differs from the canonical quadrature integrand.");
}

struct LeafBatchFixture {
    galerkin::Space space;
    SupportTree tree;
    SupportLeafSettings settings;
    std::vector<std::size_t> panelCounts;
    std::vector<unsigned> selected;
    unsigned cancellationNode{};
};
LeafBatchFixture makeLeafBatchFixture() {
    SurfaceMesh mesh;
    // Closed, oriented bipyramids provide genuine P1 stars of degrees 4/5/12
    // and a small, long-running 96-panel star for delayed cancellation.
    for (unsigned body = 0; body < 93; ++body) {
        const unsigned sides = body == 90 ? 5 : body == 91 ? 12 : body == 92 ? 96 : 4;
        const auto first = static_cast<std::uint32_t>(mesh.vertices.size());
        const Vec3 center{.04 * (body % 8), .04 * ((body / 8) % 4), .04 * (body / 32)};
        mesh.vertices.push_back({center.x, center.y, center.z + .003});
        mesh.vertices.push_back({center.x, center.y, center.z - .003});
        for (unsigned i = 0; i < sides; ++i) {
            const double phi = .17 + 2 * std::numbers::pi * i / sides;
            mesh.vertices.push_back(
                {center.x + .003 * std::cos(phi), center.y + .003 * std::sin(phi), center.z});
        }
        for (unsigned i = 0; i < sides; ++i) {
            const auto a = first + 2 + i, b = first + 2 + (i + 1) % sides;
            mesh.triangles.push_back({{first, a, b},
                                      static_cast<std::uint32_t>(mesh.triangles.size() + 1),
                                      BoundaryCondition::driven});
            mesh.triangles.push_back({{first + 1, b, a},
                                      static_cast<std::uint32_t>(mesh.triangles.size() + 1),
                                      BoundaryCondition::driven});
        }
    }
    LeafBatchFixture fixture;
    fixture.space = galerkin::buildSpace(mesh);
    SupportTreeSettings treeSettings;
    treeSettings.leafCapacity = 1;
    fixture.tree = makeSupportTree(fixture.space, treeSettings);
    auto &settings = fixture.settings;
    std::vector<bool> inactive(fixture.space.panels.size());
    bool empty = false;
    for (unsigned id = 0; id < fixture.tree.nodes.size(); ++id) {
        const auto &node = fixture.tree.nodes[id];
        if (!node.leaf())
            continue;
        require(node.pressureCount <= 1 && node.fluxCount <= 1,
                "Leaf batch fixture reached a depth limit instead of splitting actual supports.");
        std::set<unsigned> star;
        if (node.pressureCount) {
            const auto dof = fixture.tree.pressureOrder[node.firstPressure];
            for (auto j = fixture.space.incidenceOffsets[dof];
                 j < fixture.space.incidenceOffsets[dof + 1]; ++j)
                star.insert(fixture.space.incidences[j].panel);
        }
        // Each panel has its own physical face tag. Wall faces may be absent
        // from DP0 while their geometry remains in every incident P1 star.
        for (unsigned j = 0; j < node.fluxCount; ++j) {
            const auto p = fixture.tree.fluxOrder[node.firstFlux + j];
            if ((node.pressureCount && !star.contains(p)) || (!node.pressureCount && !empty)) {
                inactive[p] = true;
                empty |= !node.pressureCount;
            }
        }
    }
    settings.panelToActiveFaceColumn.resize(inactive.size(), noSupportNode);
    for (unsigned p = 0; p < inactive.size(); ++p) {
        if (inactive[p])
            fixture.space.panels[p].boundary = BoundaryCondition::wall;
        else
            settings.panelToActiveFaceColumn[p] = settings.activeFaceCount++;
    }
    std::map<std::size_t, unsigned> representatives;
    fixture.panelCounts.resize(fixture.tree.nodes.size());
    for (unsigned id = 0; id < fixture.tree.nodes.size(); ++id) {
        const auto &node = fixture.tree.nodes[id];
        if (!node.leaf())
            continue;
        std::set<unsigned> panels;
        for (unsigned d = 0; d < node.pressureCount; ++d) {
            const auto dof = fixture.tree.pressureOrder[node.firstPressure + d];
            for (auto j = fixture.space.incidenceOffsets[dof];
                 j < fixture.space.incidenceOffsets[dof + 1]; ++j)
                panels.insert(fixture.space.incidences[j].panel);
        }
        for (unsigned j = 0; j < node.fluxCount; ++j) {
            const auto p = fixture.tree.fluxOrder[node.firstFlux + j];
            if (!inactive[p])
                panels.insert(p);
        }
        fixture.panelCounts[id] = panels.size();
        representatives.try_emplace(panels.size(), id);
    }
    for (const unsigned count : {0u, 1u, 4u, 5u, 12u, 96u}) {
        const auto found = representatives.find(count);
        require(found != representatives.end(),
                "Real leaf fixture lacks a required empty/C4/panel-slice boundary.");
        settings.nodes.push_back(found->second);
    }
    fixture.selected = settings.nodes;
    fixture.cancellationNode = representatives.at(96);
    // 521 is prime and exceeds the hard 512-block wave cap. For any wave
    // capacity 2..512 it also forces a final partial wave.
    std::set<unsigned> chosen(settings.nodes.begin(), settings.nodes.end());
    for (unsigned id = 0; id < fixture.tree.nodes.size() && settings.nodes.size() < 521; ++id)
        if (fixture.tree.nodes[id].leaf() && chosen.insert(id).second)
            settings.nodes.push_back(id);
    require(settings.nodes.size() == 521, "Leaf fixture did not cross a real occupancy wave.");
    settings.order = 1;
    settings.wavenumber = 200;
    settings.normalizationRadius = .03;
    return fixture;
}

void checkIntegratedLeafBatches(cudaStream_t stream) {
    auto fixture = makeLeafBatchFixture();
    const auto exercise = [&](SupportLeafSettings settings, bool wave) {
        constexpr std::size_t fullBudget = 16 * 1024 * 1024;
        std::size_t partialBudget = 0;
        for (unsigned pass = 0; pass < (wave ? 3u : 2u); ++pass) {
            const auto budget = pass == 0 ? fullBudget : pass == 1 ? partialBudget : 0;
            auto family = std::make_shared<SupportLeafMomentFamily>(budget);
            for (const unsigned degree : {4u, 12u}) {
                settings.panelRule = galerkin::triangleRule(degree);
                settings.maximumMomentCacheBytes = 0;
                settings.momentFamily.reset();
                SupportLeafOperator canonical(fixture.space, fixture.tree, settings, {}, stream);
                settings.maximumMomentCacheBytes = budget;
                settings.momentFamily = family;
                SupportLeafOperator integrated(fixture.space, fixture.tree, settings, {}, stream);
                const auto stats = integrated.momentStatistics();
                if (budget) {
                    require(stats.bytes && stats.bytes <= budget &&
                                stats.topologyReused == (degree == 12) &&
                                stats.constructionMilliseconds > 0,
                            "Sliced family rule change did not retain topology and reintegrate.");
                    if (pass == 0) {
                        require(stats.cachedPressureDofs == stats.pressureDofs &&
                                    stats.cachedFaceSlots == stats.faceSlots,
                                "Sliced full family did not cache all P1/DP0 rows.");
                        const auto p = std::max(settings.order, settings.maximumMomentOrder);
                        const auto positive = std::size_t(p + 1) * (p + 2) / 2;
                        const auto metadata = stats.bytes -
                                              stats.cachedPressureDofs * positive * 32 -
                                              stats.cachedFaceSlots * positive * 16;
                        partialBudget = metadata + (stats.pressureDofs / 2) * positive * 32;
                    } else
                        require(stats.cachedPressureDofs > 0 &&
                                    stats.cachedPressureDofs < stats.pressureDofs,
                                "Sliced partial family failed to retain a canonical P1 remainder.");
                } else
                    require(!stats.bytes, "Zero cache budget retained sliced moments.");
                const auto states = settings.nodes.size() * canonical.stateSize();
                TestDeviceArray p(fixture.space.vertices.size()), q(settings.activeFaceCount),
                    actual(states), expected(states), locals(states),
                    tested(fixture.space.vertices.size()), reference(fixture.space.vertices.size());
                for (unsigned rhs = 0; rhs < 3; ++rhs) {
                    std::vector<Complex> pressure(fixture.space.vertices.size()),
                        flux(settings.activeFaceCount), local(states);
                    if (rhs != 2) {
                        for (std::size_t j = 0; j < pressure.size(); ++j)
                            pressure[j] = {std::sin(.31 * j + rhs), std::cos(.47 * j - rhs)};
                        for (std::size_t j = 1; j < flux.size(); ++j)
                            flux[j] = {.7 * std::cos(.23 * j + rhs), -.4 * std::sin(.17 * j - rhs)};
                        for (std::size_t j = 0; j < local.size(); ++j)
                            local[j] = {.1 * std::sin(.17 * j), .2 * std::cos(.23 * j)};
                    }
                    p.upload(pressure, stream);
                    q.upload(flux, stream);
                    locals.upload(local, stream);
                    integrated.emitPressureDevice(p.values, actual.values);
                    canonical.emitPressureDevice(p.values, expected.values);
                    compareIntegratedLeafValues(actual.download(), expected.download());
                    integrated.emitFaceFluxDevice(q.values, actual.values);
                    canonical.emitFaceFluxDevice(q.values, expected.values);
                    compareIntegratedLeafValues(actual.download(), expected.download());
                    tested.upload(std::vector<Complex>(pressure.size(), {.3, -.2}), stream);
                    reference.upload(std::vector<Complex>(pressure.size(), {.3, -.2}), stream);
                    integrated.testDevice(locals.values, tested.values);
                    canonical.testDevice(locals.values, reference.values);
                    compareIntegratedLeafValues(tested.download(), reference.download());
                }
            }
        }
    };
    exercise(fixture.settings, true);
    auto small = fixture.settings;
    small.nodes.assign(fixture.selected.begin(), fixture.selected.begin() + 5);
    for (const unsigned order : {31u, 33u, 36u, 128u}) {
        small.order = order;
        small.maximumMomentOrder = order;
        exercise(small, false);
    }
    auto interrupted = fixture.settings;
    interrupted.nodes = {fixture.cancellationNode};
    interrupted.order = interrupted.maximumMomentOrder = 128;
    interrupted.panelRule = galerkin::triangleRule(12);
    interrupted.maximumMomentCacheBytes = 16 * 1024 * 1024;
    interrupted.momentFamily =
        std::make_shared<SupportLeafMomentFamily>(interrupted.maximumMomentCacheBytes);
    const auto warmStart = std::chrono::steady_clock::now();
    double momentMilliseconds = 0;
    {
        SupportLeafOperator warm(fixture.space, fixture.tree, interrupted, {}, stream);
        momentMilliseconds = warm.momentStatistics().constructionMilliseconds;
        require(momentMilliseconds > 0,
                "Cancellation control did not build actual sliced moments.");
    }
    const double warmMilliseconds =
        std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - warmStart)
            .count();
    interrupted.momentFamily =
        std::make_shared<SupportLeafMomentFamily>(interrupted.maximumMomentCacheBytes);
    std::stop_source stopping;
    const auto delay = std::chrono::duration<double, std::milli>(
        std::max(.5, warmMilliseconds - .5 * momentMilliseconds));
    bool stopped = false;
    const auto cancelStart = std::chrono::steady_clock::now();
    {
        // The warm control anchors a genuinely delayed constructor cancellation.
        // It cannot identify the exact slice without an independent GPU trace.
        std::jthread requester([&] {
            std::this_thread::sleep_for(delay);
            stopping.request_stop();
        });
        try {
            SupportLeafOperator invalid(fixture.space, fixture.tree, interrupted,
                                        stopping.get_token(), stream);
        } catch (const std::runtime_error &error) {
            stopped = std::string_view(error.what()).find("cancelled") != std::string_view::npos;
        }
    }
    require(stopped && interrupted.momentFamily->deviceBytes() == 0 &&
                cudaStreamQuery(stream) == cudaSuccess,
            "Delayed sliced construction published partial moments or failed to drain.");
    require(std::chrono::steady_clock::now() - cancelStart < std::chrono::seconds(15),
            "Sliced moment cancellation exceeded the existing latency budget.");
    interrupted.order = interrupted.maximumMomentOrder = 1;
    SupportLeafOperator rebuilt(fixture.space, fixture.tree, interrupted, {}, stream);
    require(rebuilt.momentStatistics().bytes &&
                rebuilt.momentStatistics().constructionMilliseconds > 0,
            "Delayed cancelled family could not rebuild complete sliced moments.");
    auto plain = interrupted;
    plain.maximumMomentCacheBytes = 0;
    plain.momentFamily.reset();
    SupportLeafOperator canonical(fixture.space, fixture.tree, plain, {}, stream);
    TestDeviceArray pressure(fixture.space.vertices.size()), actual(rebuilt.stateSize()),
        expected(canonical.stateSize());
    std::vector<Complex> pressureValues(fixture.space.vertices.size());
    for (std::size_t j = 0; j < pressureValues.size(); ++j)
        pressureValues[j] = {std::sin(.31 * j), std::cos(.47 * j)};
    pressure.upload(pressureValues, stream);
    rebuilt.emitPressureDevice(pressure.values, actual.values);
    canonical.emitPressureDevice(pressure.values, expected.values);
    compareIntegratedLeafValues(actual.download(), expected.download());
    std::cout << "CUDA integrated 521 actual leaves, panels0/1/4/5/12/96, q4/q12, "
                 "storage31/33/36/128, full/partial/none and delayed cancel/rebuild PASS; "
              << "warm/delay_ms=" << warmMilliseconds << '/' << delay.count() << '\n';
}

void checkIntegratedLeafMoments() {
    const auto fixture = test::makeWidebandFarBlockFixture(200, 4);
    SupportTree tree;
    tree.nodes.resize(2);
    tree.nodes[0].center = fixture.sourceCenter;
    tree.nodes[1].center = fixture.targetCenter;
    for (unsigned i = 0; i < 8; ++i) {
        tree.pressureOrder.push_back(i);
        tree.fluxOrder.push_back(i);
        tree.pressureLeaf.push_back(i / 4);
        tree.fluxLeaf.push_back(i / 4);
    }
    for (unsigned i = 0; i < 2; ++i) {
        tree.nodes[i].firstPressure = tree.nodes[i].firstFlux = 4 * i;
        tree.nodes[i].pressureCount = tree.nodes[i].fluxCount = 4;
    }
    cudaStream_t stream{};
    require(cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking) == cudaSuccess,
            "Integrated leaf test stream creation failed.");
    struct Guard {
        cudaStream_t stream;
        ~Guard() {
            cudaStreamDestroy(stream);
        }
    } guard{stream};
    const auto compare = compareIntegratedLeafValues;
    checkIntegratedLeafBatches(stream);
    for (const unsigned order : {0u, 1u, 20u, 40u, 80u, 126u, 128u}) {
        SupportLeafSettings settings;
        settings.nodes = {0};
        settings.order = order;
        settings.wavenumber = 200;
        settings.normalizationRadius = .03;
        settings.panelRule = fixture.triangleRule;
        settings.activeFaceCount = 2;
        // Several original DP0 panels share a face; the other body's panels
        // cannot leak into this leaf's aggregated face signature.
        settings.panelToActiveFaceColumn = {0, 0, 1, 1, 0, 0, 1, 1};
        SupportLeafOperator canonical(fixture.space, tree, settings, {}, stream);
        const auto state = canonical.stateSize();
        TestDeviceArray input(8), flux(2), original(state), actual(state), locals(state),
            expected(8), result(8);
        std::vector<Complex> localValues(state);
        for (std::size_t i = 0; i < state; ++i)
            localValues[i] = {.1 * std::sin(.17 * i), .2 * std::cos(.23 * i)};
        locals.upload(localValues, stream);
        const auto positive = std::size_t(order + 1) * (order + 2) / 2;
        for (const auto budget :
             {std::size_t{1}, 2 * positive * 32 + 512, std::size_t{16 * 1024 * 1024}}) {
            auto cachedSettings = settings;
            cachedSettings.maximumMomentCacheBytes = budget;
            cachedSettings.maximumMomentOrder = std::min(128u, order + 2);
            cachedSettings.momentFamily = std::make_shared<SupportLeafMomentFamily>(budget);
            SupportLeafOperator cached(fixture.space, tree, cachedSettings, {}, stream);
            const auto stats = cached.momentStatistics();
            require(stats.bytes <= budget,
                    "Integrated signatures exceeded their exact payload cap.");
            require(cachedSettings.momentFamily->deviceBytes() == stats.bytes &&
                        cached.deviceBytes() == cached.ownDeviceBytes() + stats.bytes,
                    "Family signatures were duplicated or omitted from device accounting.");
            {
                auto copiedRule = cachedSettings;
                require(copiedRule.panelRule.data() != cachedSettings.panelRule.data(),
                        "Family identity gate did not allocate a separate equal-valued rule.");
                SupportLeafOperator reused(fixture.space, tree, copiedRule, {}, stream);
                require(reused.momentStatistics().constructionMilliseconds == 0 &&
                            reused.momentStatistics().bytes == stats.bytes &&
                            cachedSettings.momentFamily->deviceBytes() == stats.bytes,
                        "Identical family degrees were reintegrated instead of borrowed.");
            }
            if (stats.bytes) {
                auto lowered = cachedSettings;
                lowered.maximumMomentCacheBytes = 1;
                SupportLeafOperator uncached(fixture.space, tree, lowered, {}, stream);
                require(uncached.momentStatistics().bytes == 0 &&
                            cachedSettings.momentFamily->deviceBytes() == stats.bytes,
                        "A reused family entry ignored a lowered per-leaf cache cap.");
                input.upload(fixture.pressure[1], stream);
                uncached.emitPressureDevice(input.values, actual.values);
                canonical.emitPressureDevice(input.values, original.values);
                compare(actual.download(), original.download());
                auto changed = cachedSettings;
                changed.panelRule = galerkin::triangleRule(6);
                bool leased = false;
                try {
                    SupportLeafOperator invalid(fixture.space, tree, changed, {}, stream);
                } catch (const std::runtime_error &) {
                    leased = true;
                }
                require(leased && cachedSettings.momentFamily->deviceBytes() == stats.bytes,
                        "An active family lease was mutated by a different integration rule.");
                for (unsigned key = 0; key < 3; ++key) {
                    auto changedKey = cachedSettings;
                    if (key == 0)
                        changedKey.wavenumber = 201;
                    if (key == 1)
                        changedKey.normalizationRadius = .031;
                    if (key == 2)
                        changedKey.panelToActiveFaceColumn[0] = 1;
                    leased = false;
                    try {
                        SupportLeafOperator invalid(fixture.space, tree, changedKey, {}, stream);
                    } catch (const std::runtime_error &) {
                        leased = true;
                    }
                    require(leased && cachedSettings.momentFamily->deviceBytes() == stats.bytes,
                            "A family lease accepted a changed k/radius/physical face assignment.");
                }
            }
            if (budget == 1)
                require(!stats.bytes,
                        "Insufficient optional cache budget did not keep the canonical integrand.");
            if (budget == 16 * 1024 * 1024)
                require(stats.cachedPressureDofs == 4 && stats.cachedFaceSlots == 2,
                        "Full integrated cache did not cover complete P1 and leaf/face sources.");
            if (order >= 20 && budget == 2 * positive * 32 + 512)
                require(stats.cachedPressureDofs > 0 && stats.cachedPressureDofs < 4,
                        "Partial integrated cache test did not execute an uncached remainder.");
            for (unsigned rhs = 0; rhs < 2; ++rhs) {
                input.upload(fixture.pressure[rhs], stream);
                flux.upload({{.7 + .1 * rhs, -.2}, {-.4, .8 - .2 * rhs}}, stream);
                canonical.emitPressureDevice(input.values, original.values);
                cached.emitPressureDevice(input.values, actual.values);
                compare(actual.download(), original.download());
                canonical.emitFaceFluxDevice(flux.values, original.values);
                cached.emitFaceFluxDevice(flux.values, actual.values);
                compare(actual.download(), original.download());
                expected.upload(std::vector<Complex>(8, {.3, -.2}), stream);
                result.upload(std::vector<Complex>(8, {.3, -.2}), stream);
                canonical.testDevice(locals.values, expected.values);
                cached.testDevice(locals.values, result.values);
                compare(result.download(), expected.download());
            }
            input.upload(std::vector<Complex>(8), stream);
            flux.upload(std::vector<Complex>(2), stream);
            cached.emitPressureDevice(input.values, actual.values);
            compare(actual.download(), std::vector<Complex>(state));
            cached.emitFaceFluxDevice(flux.values, actual.values);
            compare(actual.download(), std::vector<Complex>(state));
            locals.upload(std::vector<Complex>(state), stream);
            result.upload(std::vector<Complex>(8), stream);
            cached.testDevice(locals.values, result.values);
            compare(result.download(), std::vector<Complex>(8));
            locals.upload(localValues, stream);
            std::stop_source stopped;
            stopped.request_stop();
            bool rejected = false;
            try {
                cached.emitPressureDevice(input.values, actual.values, stopped.get_token());
            } catch (const std::runtime_error &) {
                rejected = true;
            }
            require(rejected, "Integrated leaf ignored cancellation.");
            auto insufficient = cachedSettings;
            insufficient.maximumDeviceBytes = canonical.deviceBytes() - 1;
            rejected = false;
            try {
                SupportLeafOperator invalid(fixture.space, tree, insufficient, {}, stream);
            } catch (const std::runtime_error &) {
                rejected = true;
            }
            require(rejected, "Leaf geometry/device budget failed to reject before construction.");
        }
        std::cout << "CUDA integrated leaf moments "
                     "prefix/full/partial/uncached/face/zero/reuse/budget/cancel order="
                  << order << '\n';
    }
    SupportLeafSettings miss;
    miss.nodes = {0};
    miss.order = 2;
    miss.wavenumber = 200;
    miss.normalizationRadius = .03;
    miss.panelRule = fixture.triangleRule;
    miss.activeFaceCount = 2;
    miss.panelToActiveFaceColumn = {0, 0, 1, 1, 0, 0, 1, 1};
    miss.maximumMomentCacheBytes = 16 * 1024 * 1024;
    miss.momentFamily = std::make_shared<SupportLeafMomentFamily>(miss.maximumMomentCacheBytes);
    std::size_t previousBytes = 0;
    for (unsigned pass = 0; pass < 5; ++pass) {
        if (pass == 1)
            miss.maximumMomentOrder = 6;
        if (pass == 2)
            miss.panelRule = galerkin::triangleRule(6);
        if (pass == 3) {
            miss.panelRule.front().barycentric[0] += 1e-7;
            miss.panelRule.front().barycentric[1] -= 1e-7;
        }
        if (pass == 4)
            miss.panelRule.front().weight = std::nextafter(miss.panelRule.front().weight,
                                                           std::numeric_limits<double>::infinity());
        SupportLeafOperator integrated(fixture.space, tree, miss, {}, stream);
        const auto stats = integrated.momentStatistics();
        require(stats.bytes && stats.constructionMilliseconds > 0 &&
                    stats.bytes == miss.momentFamily->deviceBytes(),
                "A genuine degree/rule miss did not integrate a new bounded family bank.");
        if (pass == 1)
            require(stats.bytes > previousBytes, "A higher family degree was only relabelled.");
        previousBytes = stats.bytes;
        auto plainSettings = miss;
        plainSettings.maximumMomentCacheBytes = 0;
        plainSettings.momentFamily.reset();
        SupportLeafOperator plain(fixture.space, tree, plainSettings, {}, stream);
        TestDeviceArray p(8), q(2), actual(integrated.stateSize()), expected(plain.stateSize());
        for (unsigned rhs = 0; rhs < 2; ++rhs) {
            p.upload(fixture.pressure[rhs], stream);
            q.upload({{.7, -.2}, {-.4, .8}}, stream);
            integrated.emitPressureDevice(p.values, actual.values);
            plain.emitPressureDevice(p.values, expected.values);
            compare(actual.download(), expected.download());
            integrated.emitFaceFluxDevice(q.values, actual.values);
            plain.emitFaceFluxDevice(q.values, expected.values);
            compare(actual.download(), expected.download());
        }
        auto grow = miss;
        grow.maximumMomentOrder = 7;
        bool leased = false;
        try {
            SupportLeafOperator invalid(fixture.space, tree, grow, {}, stream);
        } catch (const std::runtime_error &) {
            leased = true;
        }
        require(leased && miss.momentFamily->deviceBytes() == stats.bytes,
                "Growing a leased family silently changed its table stride.");
        std::stop_source stopped;
        stopped.request_stop();
        bool cancelled = false;
        try {
            SupportLeafOperator invalid(fixture.space, tree, grow, stopped.get_token(), stream);
        } catch (const std::runtime_error &) {
            cancelled = true;
        }
        require(cancelled && miss.momentFamily->deviceBytes() == stats.bytes,
                "Cancelled family preparation destroyed a valid reusable bank.");
    }
    // Split a complete tetrahedral P1 star between two leaves, with distinct
    // expansion origins. Each leaf still contains all incident physical panels;
    // original DP0 panels belong to exactly one of the two leaves.
    auto splitTree = tree;
    splitTree.nodes.resize(3);
    splitTree.pressureLeaf = {0, 0, 1, 1, 2, 2, 2, 2};
    splitTree.fluxLeaf = {0, 1, 0, 1, 2, 2, 2, 2};
    splitTree.fluxOrder = {0, 2, 1, 3, 4, 5, 6, 7};
    splitTree.nodes[0].center = {fixture.sourceCenter.x + .001, fixture.sourceCenter.y,
                                 fixture.sourceCenter.z};
    splitTree.nodes[1].center = {fixture.sourceCenter.x, fixture.sourceCenter.y - .001,
                                 fixture.sourceCenter.z};
    splitTree.nodes[2].center = fixture.targetCenter;
    for (unsigned node = 0; node < 3; ++node) {
        splitTree.nodes[node].firstPressure = splitTree.nodes[node].firstFlux = 2 * node;
        splitTree.nodes[node].pressureCount = splitTree.nodes[node].fluxCount = node < 2 ? 2 : 4;
    }
    const unsigned splitOrder = 20;
    const auto splitPositive = std::size_t(splitOrder + 1) * (splitOrder + 2) / 2;
    for (const auto budget :
         {std::size_t{0}, std::size_t{1}, 3 * splitPositive * 32 + 512, std::size_t{1024 * 1024}}) {
        auto sharedTopology = std::make_shared<SupportLeafMomentFamily>(budget);
        std::size_t topologyBytes = 0;
        for (const unsigned degree : {6u, 8u}) {
            SupportLeafSettings split;
            split.nodes = {0, 1};
            split.order = splitOrder;
            split.wavenumber = 200;
            split.normalizationRadius = .04;
            split.panelRule = galerkin::triangleRule(degree);
            split.activeFaceCount = 2;
            split.panelToActiveFaceColumn = {0, 1, 1, 0, 0, 1, 1, 0};
            SupportLeafOperator plain(fixture.space, splitTree, split, {}, stream);
            split.maximumMomentCacheBytes = budget;
            split.momentFamily = sharedTopology;
            SupportLeafOperator integrated(fixture.space, splitTree, split, {}, stream);
            const auto stats = integrated.momentStatistics();
            require(stats.bytes <= budget && sharedTopology->deviceBytes() == stats.bytes,
                    "Split-star moments/topology exceeded or duplicated their shared budget.");
            if (budget > 1) {
                require(stats.cachedPressureDofs == (budget < 1024 * 1024 ? 3u : 4u) &&
                            stats.pressureIncidences == 12 && stats.panelReferences == 8 &&
                            stats.fluxPanelReferences == 4 && stats.topologyBytes == 8 * 16 &&
                            stats.topologyHostBytes >= stats.topologyBytes,
                        "Split-star topology lost an incidence, panel copy or original DP0 panel.");
                require(stats.constructionMilliseconds > 0 &&
                            stats.constructionPanelReferences > 0 &&
                            stats.constructionPanelReferences <= 8,
                        "Independent panel rule did not genuinely integrate the active topology.");
                require(stats.topologyReused == (degree == 8) &&
                            (degree != 8 || topologyBytes == stats.topologyBytes),
                        "A rule change rebuilt topology or reused dependent moment integrals.");
                if (budget == 1024 * 1024)
                    require(stats.cachedFaceSlots == 4 && stats.constructionPanelReferences == 8,
                            "Full split-star cache did not fuse P1 and DP0 panel work.");
                topologyBytes = stats.topologyBytes;
            } else
                require(!stats.bytes && !stats.topologyBytes,
                        "An insufficient optional budget retained an unusable topology.");
            const auto state = 2 * plain.stateSize();
            TestDeviceArray p(8), q(2), actual(state), expected(state), locals(state), tested(8),
                reference(8);
            p.upload(fixture.pressure[1], stream);
            std::vector<Complex> localValues(state);
            for (std::size_t j = 0; j < state; ++j)
                localValues[j] = {.1 * std::sin(.17 * j), .2 * std::cos(.23 * j)};
            locals.upload(localValues, stream);
            integrated.emitPressureDevice(p.values, actual.values);
            plain.emitPressureDevice(p.values, expected.values);
            compare(actual.download(), expected.download());
            // A valid active face with zero velocity must survive construction
            // and later accept a nonzero complex excitation without rebuilding.
            for (const auto &velocities :
                 {std::vector<Complex>{{}, {-.4, .8}}, std::vector<Complex>{{.7, -.2}, {}},
                  std::vector<Complex>{{.7, -.2}, {-.4, .8}}}) {
                q.upload(velocities, stream);
                integrated.emitFaceFluxDevice(q.values, actual.values);
                plain.emitFaceFluxDevice(q.values, expected.values);
                compare(actual.download(), expected.download());
            }
            tested.upload(std::vector<Complex>(8, {.3, -.2}), stream);
            reference.upload(std::vector<Complex>(8, {.3, -.2}), stream);
            integrated.testDevice(locals.values, tested.values);
            plain.testDevice(locals.values, reference.values);
            compare(tested.download(), reference.download());
        }
    }
    std::cout << "CUDA integrated split-star panel topology q6/q8, full/partial/none, "
                 "complex/zero active faces PASS\n";
    auto interrupted = miss;
    interrupted.maximumMomentOrder = 128;
    interrupted.panelRule = galerkin::triangleRule(14);
    interrupted.momentFamily = std::make_shared<SupportLeafMomentFamily>(16 * 1024 * 1024);
    std::stop_source preparingStop;
    bool rejected = false;
    {
        std::jthread requester([&] {
            std::this_thread::sleep_for(std::chrono::milliseconds(1));
            preparingStop.request_stop();
        });
        try {
            SupportLeafOperator invalid(fixture.space, tree, interrupted, preparingStop.get_token(),
                                        stream);
        } catch (const std::runtime_error &) {
            rejected = true;
        }
    }
    require(rejected && interrupted.momentFamily->deviceBytes() == 0 &&
                cudaStreamQuery(stream) == cudaSuccess,
            "Interrupted moment integration published partial tables or retained inflight work.");
    interrupted.maximumMomentOrder = 6;
    interrupted.panelRule = fixture.triangleRule;
    SupportLeafOperator afterCancel(fixture.space, tree, interrupted, {}, stream);
    require(afterCancel.momentStatistics().bytes &&
                afterCancel.momentStatistics().constructionMilliseconds > 0,
            "A cancelled family could not build complete replacement moments.");
    std::cout << "CUDA integrated leaf family degree/rule miss, lease and cancellation PASS\n";
}

void checkCompleteSupportOperator() {
    auto fixture = test::makeWidebandFarBlockFixture(20);
    for (auto &panel : fixture.space.panels)
        panel.boundary = BoundaryCondition::driven;
    SupportTreeSettings treeSettings;
    treeSettings.leafCapacity = 2;
    treeSettings.nearDistanceRatio = .5;
    auto tree = std::make_shared<SupportTree>(makeSupportTree(fixture.space, treeSettings));
    require(tree->pressurePressure.representedFarPairs > 0 &&
                tree->pressureFlux.representedFarPairs > 0,
            "Complete operator fixture has no real far partition.");
    const std::array<std::uint32_t, 2> faces{5, 1};
    auto nearPlan = std::make_shared<SupportNearPlan>(
        makeSupportNearPlan(fixture.space, *tree, faces, 16 * 1024 * 1024));
    SupportNearOperatorSettings nearSettings;
    nearSettings.quadrature.nearDistanceRatio = .5;
    auto near =
        std::make_shared<SupportNearOperator>(fixture.space, *tree, *nearPlan, 20, nearSettings);
    auto pc = std::make_shared<SupportBlockJacobi>(near->pressureMatrix(), tree->pressureOrder, 32,
                                                   std::stop_token{}, near->stream());
    SupportFmmOperatorSettings settings;
    settings.role = SupportFmmOperatorRole::distancePrimary;
    settings.panelRule = fixture.triangleRule;
    settings.maximumClassCacheBytes = 64 * 1024 * 1024;
    SupportFmmOperator gpu(fixture.space, tree, nearPlan, near, pc, 20, settings);
    const auto reference =
        test::support_operator_oracle::directReference(fixture.space, 20, nearSettings.quadrature);
    TestDeviceArray pressure(8), flux(2), output(8);
    std::array<std::array<std::vector<Complex>, 2>, 3> primaryActions;
    for (unsigned rhs = 0; rhs < 3; ++rhs) {
        std::vector<Complex> p, q;
        for (unsigned i = 0; i < 8; ++i)
            p.emplace_back(std::cos(.7 * i + .2 * rhs), std::sin(.9 * i - .3 * rhs));
        q = rhs == 0   ? std::vector<Complex>{{.3, -.7}, {}}
            : rhs == 1 ? std::vector<Complex>{{.5, -.7}, {-.5, .7}}
                       : std::vector<Complex>(2);
        for (unsigned channel = 0; channel < 2; ++channel) {
            if (channel) {
                flux.upload(q, gpu.stream());
                gpu.applyFaceFluxDevice(flux.values, output.values);
            } else {
                pressure.upload(p, gpu.stream());
                gpu.applyPressureDevice(pressure.values, output.values);
            }
            const auto actual = output.download();
            primaryActions[rhs][channel] = actual;
            double error = 0, norm = 0;
            for (unsigned i = 0; i < 8; ++i) {
                Complex expected{};
                const auto &input = channel ? q : p;
                const auto &matrix = channel ? reference.b : reference.a;
                for (unsigned j = 0; j < input.size(); ++j)
                    expected += matrix[i * input.size() + j] * input[j];
                error += std::norm(actual[i] - expected);
                norm += std::norm(expected);
            }
            const double relative = norm ? std::sqrt(error / norm) : std::sqrt(error);
            std::cout << "CUDA complete near+wideband operator channel=" << channel
                      << ", rhs=" << rhs << ", relative=" << relative
                      << ", live=" << gpu.deviceBytes() << ", own=" << gpu.ownDeviceBytes() << '\n';
            require(relative < 2e-6,
                    "Complete CUDA near+wideband operator fails direct Galerkin action.");
        }
    }
    gpu.releaseDeviceState();
    {
        SupportFmmOperator native(gpu, settings.accuracySeed);
        require(native.role() == SupportFmmOperatorRole::nativeAudit &&
                    std::equal(native.levelEvidence().begin(), native.levelEvidence().end(),
                               gpu.levelEvidence().begin(), gpu.levelEvidence().end()),
                "Same-seed native reference changed the full level plan.");
        for (const auto *partition : {&tree->pressurePressure, &tree->pressureFlux})
            for (const auto pair : partition->farInteractions) {
                const auto &source = tree->nodes[pair.source], &target = tree->nodes[pair.target];
                const double distance =
                    std::hypot(target.center.x - source.center.x, target.center.y - source.center.y,
                               target.center.z - source.center.z);
                const auto found =
                    std::find_if(native.levelEvidence().begin(), native.levelEvidence().end(),
                                 [&](const auto &level) { return level.level == source.level; });
                if (found == native.levelEvidence().end())
                    continue;
                require(native.m2lSourceOrder(source.level, distance) == found->momentBandlimit &&
                            native.m2lTargetOrder(source.level, distance) == found->momentBandlimit,
                        "Same-seed native reference shortened an M2L endpoint.");
            }
        for (unsigned rhs = 0; rhs < 3; ++rhs) {
            std::vector<Complex> p;
            for (unsigned i = 0; i < 8; ++i)
                p.emplace_back(std::cos(.7 * i + .2 * rhs), std::sin(.9 * i - .3 * rhs));
            const auto q = rhs == 0   ? std::vector<Complex>{{.3, -.7}, {}}
                           : rhs == 1 ? std::vector<Complex>{{.5, -.7}, {-.5, .7}}
                                      : std::vector<Complex>(2);
            for (unsigned channel = 0; channel < 2; ++channel) {
                if (channel) {
                    flux.upload(q, native.stream());
                    native.applyFaceFluxDevice(flux.values, output.values);
                } else {
                    pressure.upload(p, native.stream());
                    native.applyPressureDevice(pressure.values, output.values);
                }
                const auto expected = output.download();
                double difference = 0, magnitude = 0;
                for (unsigned i = 0; i < 8; ++i) {
                    difference += std::norm(primaryActions[rhs][channel][i] - expected[i]);
                    magnitude += std::norm(expected[i]);
                }
                require(magnitude ? std::sqrt(difference / magnitude) < 2e-6 : difference == 0,
                        "Two-endpoint compact BM action fails the uncut same-seed operator.");
            }
        }
        native.releaseDeviceState();
    }
    std::vector<Complex> zero(8);
    pressure.upload(zero, gpu.stream());
    gpu.applyPressureDevice(pressure.values, output.values);
    for (const auto value : output.download())
        require(value == Complex{}, "Complete support operator retained a previous RHS.");
    std::stop_source cancelled;
    cancelled.request_stop();
    bool stopped = false;
    try {
        gpu.applyPressureDevice(pressure.values, output.values, cancelled.get_token());
    } catch (const std::exception &) {
        stopped = true;
    }
    require(stopped, "Complete support operator ignored cancellation.");
    gpu.releaseDeviceState();
    require(gpu.ownDeviceBytes() == 0, "Released support hierarchy retains device state.");
    SupportFmmOperator audit(gpu, settings.accuracySeed / 100);
    require(audit.role() == SupportFmmOperatorRole::nativeAudit,
            "Seed audit inherited the primary distance cutoff.");
    for (const auto pair : tree->pressurePressure.farInteractions) {
        const auto &source = tree->nodes[pair.source], &target = tree->nodes[pair.target];
        const auto distance =
            std::hypot(target.center.x - source.center.x, target.center.y - source.center.y,
                       target.center.z - source.center.z);
        const auto found =
            std::find_if(audit.levelEvidence().begin(), audit.levelEvidence().end(),
                         [&](const auto &level) { return level.level == source.level; });
        require(found != audit.levelEvidence().end() &&
                    audit.m2lSourceOrder(source.level, distance) == found->momentBandlimit &&
                    audit.m2lTargetOrder(source.level, distance) == found->momentBandlimit,
                "Seed audit used a shortened M2L endpoint.");
    }
    require(audit.sharedDeviceBytes() == gpu.sharedDeviceBytes(),
            "Audit duplicated shared near/preconditioner/family accounting.");
    require(audit.deviceBytes() == audit.ownDeviceBytes() + audit.sharedDeviceBytes() &&
                audit.stationaryDeviceBytes() >= audit.deviceBytes() &&
                audit.peakDeviceBytes() >= audit.stationaryDeviceBytes(),
            "Family residency was omitted or combined from unrelated historical phases.");
    require(audit.workStatistics().leafMomentConstructionMilliseconds == 0 &&
                audit.workStatistics().leafMomentCacheBytes ==
                    gpu.workStatistics().leafMomentCacheBytes,
            "Same-rule seed audit reintegrated an already prepared maximum-degree family.");
    const auto baseEvidence = gpu.levelEvidence(), auditEvidence = audit.levelEvidence();
    require(!std::equal(baseEvidence.begin(), baseEvidence.end(), auditEvidence.begin(),
                        auditEvidence.end()),
            "Stricter support audit did not change its actual level plan.");
    audit.releaseDeviceState();
    audit.applyPressureDevice(pressure.values, output.values);
    require(audit.role() == SupportFmmOperatorRole::nativeAudit,
            "Released seed audit rebuilt with the primary role.");
    audit.releaseDeviceState();
    {
        SupportFmmOperator quadrature(gpu, settings.accuracySeed / 100, galerkin::triangleRule(6),
                                      SupportFmmOperatorRole::nativeAudit);
        quadrature.releaseDeviceState();
        quadrature.applyPressureDevice(pressure.values, output.values);
        require(quadrature.role() == SupportFmmOperatorRole::nativeAudit,
                "Rebuilt higher-quadrature audit inherited a primary cutoff.");
        for (const auto pair : tree->pressurePressure.farInteractions) {
            const auto &source = tree->nodes[pair.source], &target = tree->nodes[pair.target];
            const auto found =
                std::find_if(quadrature.levelEvidence().begin(), quadrature.levelEvidence().end(),
                             [&](const auto &level) { return level.level == source.level; });
            const double distance =
                std::hypot(target.center.x - source.center.x, target.center.y - source.center.y,
                           target.center.z - source.center.z);
            require(found != quadrature.levelEvidence().end() &&
                        quadrature.m2lSourceOrder(source.level, distance) ==
                            found->momentBandlimit &&
                        quadrature.m2lTargetOrder(source.level, distance) == found->momentBandlimit,
                    "Higher-quadrature audit shortened a native M2L endpoint.");
        }
        quadrature.releaseDeviceState();
    }
    {
        SupportFmmOperator promoted(gpu, settings.accuracySeed, galerkin::triangleRule(6),
                                    gpu.role());
        promoted.releaseDeviceState();
        promoted.applyPressureDevice(pressure.values, output.values);
        require(promoted.role() == SupportFmmOperatorRole::distancePrimary,
                "Rule promotion/rebuild lost the primary distance-order role.");
        promoted.releaseDeviceState();
    }
    gpu.applyPressureDevice(pressure.values, output.values);
    require(gpu.role() == SupportFmmOperatorRole::distancePrimary,
            "Released primary hierarchy rebuilt with an audit role.");
    require(gpu.ownDeviceBytes() > 0, "Support hierarchy did not rebuild after release.");
    for (const auto value : output.download())
        require(value == Complex{}, "Rebuilt support hierarchy retained an earlier RHS.");
}

void checkRcrCancellation() {
    const auto space = galerkin::buildSpace(testing::sphere(3));
    const double k = 2 * std::numbers::pi * 1000 / 343;
    SupportTreeSettings treeSettings;
    treeSettings.leafCapacity = 4;
    treeSettings.nearDistanceRatio = 1;
    auto tree = std::make_shared<SupportTree>(makeSupportTree(space, treeSettings));
    const std::array<std::uint32_t, 1> faces{1};
    auto nearPlan = std::make_shared<SupportNearPlan>(
        makeSupportNearPlan(space, *tree, faces, 256 * 1024 * 1024));
    SupportNearOperatorSettings nearSettings;
    nearSettings.quadrature.nearDistanceRatio = 1;
    auto near = std::make_shared<SupportNearOperator>(space, *tree, *nearPlan, k, nearSettings);
    auto pc = std::make_shared<SupportBlockJacobi>(near->pressureMatrix(), tree->pressureOrder, 32,
                                                   std::stop_token{}, near->stream());
    SupportFmmOperatorSettings settings;
    settings.role = SupportFmmOperatorRole::distancePrimary;
    settings.panelRule = galerkin::triangleRule(4);
    SupportFmmOperator gpu(space, tree, nearPlan, near, pc, k, settings);
    TestDeviceArray pressure(space.vertices.size()), output(space.vertices.size());
    std::vector<Complex> input;
    for (unsigned i = 0; i < space.vertices.size(); ++i)
        input.emplace_back(std::sin(.31 * i), std::cos(.47 * i));
    pressure.upload(input, gpu.stream());
    gpu.applyPressureDevice(pressure.values, output.values);
    const auto reference = output.download();
    gpu.applyPressureDevice(pressure.values, output.values);
    require(gpu.workStatistics().m2lGraphLaunches > 0,
            "Cancellation probe did not warm an actual M2L graph executable.");
    const auto &work = gpu.workStatistics();
    require(work.rcrBatchesQueued && work.rcrCompletionChecks <= work.rcrBatchesQueued &&
                work.maximumRcrChunkBatches <= 32 &&
                work.maximumRcrChunkWork <= RcrQueuePolicy::maximumChunkWork,
            "Complete RCR application exceeds its bounded queue contract.");
    const double warmMilliseconds = work.emissionMilliseconds + work.upwardMilliseconds +
                                    work.translationMilliseconds + work.downwardMilliseconds +
                                    work.testingMilliseconds;
    std::stop_source stop;
    // Warm plans and inputs already exist. Request during actual application,
    // rather than only checking a token that was stopped before any launch.
    const auto delay =
        std::chrono::duration<double, std::milli>(std::clamp(warmMilliseconds / 2, 1.0, 50.0));
    std::jthread requester([&] {
        std::this_thread::sleep_for(delay);
        stop.request_stop();
    });
    bool stopped = false;
    const auto start = std::chrono::steady_clock::now();
    try {
        gpu.applyPressureDevice(pressure.values, output.values, stop.get_token());
    } catch (const std::runtime_error &error) {
        stopped = std::string_view(error.what()).find("cancelled") != std::string_view::npos;
    }
    requester.join();
    require(stopped, "In-flight RCR application did not propagate cancellation.");
    require(cudaStreamQuery(gpu.stream()) == cudaSuccess,
            "Cancelled support application returned with in-flight work on its borrowed stream.");
    require(std::chrono::steady_clock::now() - start < std::chrono::seconds(15),
            "Bounded RCR cancellation exceeded the existing cancellation latency budget.");
    gpu.applyPressureDevice(pressure.values, output.values);
    const auto actual = output.download();
    double error = 0, norm = 0;
    for (std::size_t i = 0; i < actual.size(); ++i) {
        error += std::norm(actual[i] - reference[i]);
        norm += std::norm(reference[i]);
    }
    require(std::sqrt(error / norm) < 2e-12,
            "Cancelled RCR application corrupted subsequent valid reuse.");
    std::cout << "CUDA bounded RCR in-flight cancel/drain/reuse pass, batches/checks="
              << gpu.workStatistics().rcrBatchesQueued << '/'
              << gpu.workStatistics().rcrCompletionChecks << '\n';
}

void checkMixedCompleteOperator() {
    // k=55 was selected using CPU level plans only, before evaluating any
    // physical errors: fixed mesh/partition and seeds 5e-6, 5e-7, 5e-8.
    constexpr double k = 55;
    auto geometry = test::makeMixedWidebandFixture(k, 12);
    SupportTreeSettings treeSettings;
    treeSettings.leafCapacity = 2;
    treeSettings.nearDistanceRatio = .5;
    auto tree = std::make_shared<SupportTree>(makeSupportTree(geometry.space, treeSettings));
    std::vector<unsigned char> pressureMask(32 * 32), fluxMask(32 * 32);
    for (const auto pair : tree->pressurePressure.farInteractions) {
        const auto &target = tree->nodes[pair.target], &source = tree->nodes[pair.source];
        for (unsigned i = 0; i < target.pressureCount; ++i)
            for (unsigned j = 0; j < source.pressureCount; ++j)
                pressureMask[32 * tree->pressureOrder[target.firstPressure + i] +
                             tree->pressureOrder[source.firstPressure + j]] = 1;
    }
    for (const auto pair : tree->pressureFlux.farInteractions) {
        const auto &target = tree->nodes[pair.target], &source = tree->nodes[pair.source];
        for (unsigned i = 0; i < target.pressureCount; ++i)
            for (unsigned j = 0; j < source.fluxCount; ++j)
                fluxMask[32 * tree->pressureOrder[target.firstPressure + i] +
                         tree->fluxOrder[source.firstFlux + j]] = 1;
    }
    auto coarse = test::makeMixedWidebandFixture(k, 12, pressureMask, fluxMask);
    auto fixture = test::makeMixedWidebandFixture(k, 16, pressureMask, fluxMask);
    auto relativeTargets = [&](const auto &actual, const auto &expected) {
        double error = 0, norm = 0;
        for (unsigned i = 16; i < 32; ++i) {
            error += std::norm(actual[i] - expected[i]);
            norm += std::norm(expected[i]);
        }
        require(norm > 0, "Mixed physical gate has no nonzero isolated far contribution.");
        return std::sqrt(error / norm);
    };
    for (unsigned rhs = 0; rhs < 2; ++rhs) {
        require(relativeTargets(coarse.expectedPressure[rhs], fixture.expectedPressure[rhs]) < 2e-9,
                "Mixed weakform pressure oracle is not quadrature-converged.");
        require(relativeTargets(coarse.expectedFlux[rhs], fixture.expectedFlux[rhs]) < 2e-9,
                "Mixed weakform flux oracle is not quadrature-converged.");
    }
    const std::array<std::uint32_t, 2> faces{1, 5};
    auto nearPlan = std::make_shared<SupportNearPlan>(
        makeSupportNearPlan(fixture.space, *tree, faces, 16 * 1024 * 1024));
    SupportNearOperatorSettings nearSettings;
    nearSettings.quadrature.nearDistanceRatio = .5;
    auto near =
        std::make_shared<SupportNearOperator>(fixture.space, *tree, *nearPlan, k, nearSettings);
    auto pc = std::make_shared<SupportBlockJacobi>(near->pressureMatrix(), tree->pressureOrder, 32,
                                                   std::stop_token{}, near->stream());
    SupportFmmOperatorSettings settings;
    settings.role = SupportFmmOperatorRole::distancePrimary;
    settings.panelRule = fixture.triangleRule;
    constexpr std::array<double, 3> seeds{5e-6, 5e-7, 5e-8};
    settings.accuracySeed = seeds[0];
    SupportFmmOperator base(fixture.space, tree, nearPlan, near, pc, k, settings);
    std::vector<LevelRepresentation> previous;
    bool changedRepresentation = false;
    bool shortenedBothEndpoints = false;
    for (unsigned seedIndex = 0; seedIndex < seeds.size(); ++seedIndex) {
        std::unique_ptr<SupportFmmOperator> audit;
        if (seedIndex)
            audit = std::make_unique<SupportFmmOperator>(base, seeds[seedIndex]);
        auto &gpu = audit ? *audit : base;
        require(gpu.role() == (seedIndex ? SupportFmmOperatorRole::nativeAudit
                                         : SupportFmmOperatorRole::distancePrimary),
                "Mixed hierarchy has the wrong primary/audit role.");
        for (const auto *partition : {&tree->pressurePressure, &tree->pressureFlux})
            for (const auto pair : partition->farInteractions) {
                const auto &source = tree->nodes[pair.source], &target = tree->nodes[pair.target];
                const double distance =
                    std::hypot(target.center.x - source.center.x, target.center.y - source.center.y,
                               target.center.z - source.center.z);
                const auto found =
                    std::find_if(gpu.levelEvidence().begin(), gpu.levelEvidence().end(),
                                 [&](const auto &level) { return level.level == source.level; });
                if (found == gpu.levelEvidence().end())
                    continue;
                const auto sourceOrder = gpu.m2lSourceOrder(source.level, distance),
                           targetOrder = gpu.m2lTargetOrder(source.level, distance);
                if (seedIndex)
                    require(sourceOrder == found->momentBandlimit &&
                                targetOrder == found->momentBandlimit,
                            "Mixed independent audit shortened a native M2L endpoint.");
                else
                    shortenedBothEndpoints |= sourceOrder < found->momentBandlimit &&
                                              targetOrder < found->momentBandlimit;
            }
        std::vector<LevelRepresentation> representation(tree->levels.size(),
                                                        LevelRepresentation::inactive);
        bool hf = false, lf = false, handoff = false;
        for (const auto &level : gpu.levelEvidence()) {
            const bool plane = level.representation == BemFastMultipoleRepresentation::planeWave;
            representation[level.level] =
                plane ? LevelRepresentation::planeWave : LevelRepresentation::regularSpherical;
            hf |= plane;
            lf |= !plane;
        }
        for (const auto &node : tree->nodes)
            if (node.parent != noSupportNode && node.activeChannels)
                handoff |= representation[node.level] == LevelRepresentation::regularSpherical &&
                           representation[tree->nodes[node.parent].level] ==
                               LevelRepresentation::planeWave;
        require(hf && lf && handoff,
                "Complete mixed gate did not execute actual LF/HF levels and their handoff.");
        if (!previous.empty())
            for (unsigned level = 0; level < representation.size(); ++level)
                changedRepresentation |=
                    previous[level] == LevelRepresentation::planeWave &&
                    representation[level] == LevelRepresentation::regularSpherical;
        previous = representation;
        TestDeviceArray pressure(32), flux(2), output(32), nearOutput(32);
        std::array<std::array<std::vector<Complex>, 2>, 2> primaryActions;
        for (unsigned rhs = 0; rhs < 2; ++rhs)
            for (unsigned channel = 0; channel < 2; ++channel) {
                if (channel) {
                    flux.upload({fixture.flux[rhs][0], fixture.flux[rhs][16]}, gpu.stream());
                    gpu.applyFaceFluxDevice(flux.values, output.values);
                } else {
                    pressure.upload(fixture.pressure[rhs], gpu.stream());
                    gpu.applyPressureDevice(pressure.values, output.values);
                }
                auto actual = output.download();
                if (channel)
                    near->applyFaceFluxDevice(flux.values, nearOutput.values);
                else
                    near->applyPressureDevice(pressure.values, nearOutput.values);
                const auto exactNear = nearOutput.download();
                for (unsigned i = 0; i < 32; ++i)
                    actual[i] -= exactNear[i];
                if (!seedIndex)
                    primaryActions[rhs][channel] = actual;
                const double relative = relativeTargets(
                    actual, channel ? fixture.expectedFlux[rhs] : fixture.expectedPressure[rhs]);
                std::cout << "CUDA complete mixed LF/HF seed=" << seeds[seedIndex]
                          << ", far-only channel=" << channel << ", rhs=" << rhs
                          << ", relative=" << relative << ", live=" << gpu.deviceBytes() << '\n';
                require(
                    relative < 2e-6,
                    "Complete mixed LF/HF hierarchy fails the isolated independent weak-BM field.");
            }
        gpu.releaseDeviceState();
        require(gpu.ownDeviceBytes() == 0,
                "Mixed seed gate retained a hierarchy before constructing the next seed.");
        if (!seedIndex) {
            require(shortenedBothEndpoints,
                    "Same-seed mixed comparison has no actual two-endpoint cutoff.");
            SupportFmmOperator native(base, seeds[0]);
            require(native.role() == SupportFmmOperatorRole::nativeAudit &&
                        std::equal(native.levelEvidence().begin(), native.levelEvidence().end(),
                                   base.levelEvidence().begin(), base.levelEvidence().end()),
                    "Mixed same-seed native reference changed its level plan.");
            for (const auto *partition : {&tree->pressurePressure, &tree->pressureFlux})
                for (const auto pair : partition->farInteractions) {
                    const auto &source = tree->nodes[pair.source],
                               &target = tree->nodes[pair.target];
                    const double distance = std::hypot(target.center.x - source.center.x,
                                                       target.center.y - source.center.y,
                                                       target.center.z - source.center.z);
                    const auto found = std::find_if(
                        native.levelEvidence().begin(), native.levelEvidence().end(),
                        [&](const auto &level) { return level.level == source.level; });
                    if (found == native.levelEvidence().end())
                        continue;
                    require(
                        native.m2lSourceOrder(source.level, distance) == found->momentBandlimit &&
                            native.m2lTargetOrder(source.level, distance) == found->momentBandlimit,
                        "Mixed same-seed native reference shortened an M2L endpoint.");
                }
            for (unsigned rhs = 0; rhs < 2; ++rhs)
                for (unsigned channel = 0; channel < 2; ++channel) {
                    if (channel) {
                        flux.upload({fixture.flux[rhs][0], fixture.flux[rhs][16]}, native.stream());
                        native.applyFaceFluxDevice(flux.values, output.values);
                        near->applyFaceFluxDevice(flux.values, nearOutput.values);
                    } else {
                        pressure.upload(fixture.pressure[rhs], native.stream());
                        native.applyPressureDevice(pressure.values, output.values);
                        near->applyPressureDevice(pressure.values, nearOutput.values);
                    }
                    auto actual = output.download();
                    const auto exactNear = nearOutput.download();
                    for (unsigned i = 0; i < 32; ++i)
                        actual[i] -= exactNear[i];
                    require(relativeTargets(actual, primaryActions[rhs][channel]) < 2e-6,
                            "Shortened mixed BM action fails the uncut same-seed reference.");
                }
            native.releaseDeviceState();
        }
    }
    require(changedRepresentation,
            "Mixed fixed-mesh seed gate did not execute an actual HF-to-LF level change.");
    require(shortenedBothEndpoints,
            "Mixed BM accuracy gate did not exercise a real two-endpoint distance class.");
}

void checkBoundedCpuPlans() {
    TranslationPlanCache cache(8 * 1024 * 1024, 4);
    std::vector<AxialPlanRequest> requests;
    for (const double distance : {.2, .3, .4, .5})
        requests.push_back(
            {ScaledAxialKind::multipoleToLocal, 6, 8, 40, 10, .01, .02, distance, 0});
    cache.prepare(requests);
    auto first = cache.axial(requests.front()), second = cache.axial(requests.front());
    require(first == second, "Identical immutable axial CPU plans were not shared.");
    require(cache.liveBytes() > 0 && cache.peakBytes() <= 8 * 1024 * 1024,
            "Parallel axial CPU preparation exceeded its bounded reservation.");
    std::stop_source cancelled;
    cancelled.request_stop();
    bool stopped = false;
    try {
        cache.prepare(requests, cancelled.get_token());
    } catch (const std::exception &) {
        stopped = true;
    }
    require(stopped, "CPU axial worker pool ignored cancellation.");
    TranslationPlanCache small(1, 1);
    bool budget = false;
    try {
        (void)small.axial(requests.front());
    } catch (const std::exception &) {
        budget = true;
    }
    require(budget, "CPU axial builder did not reject insufficient budget before allocation.");
}

void printWork(const SupportFmmOperator &gpu, const char *label) {
    const auto &work = gpu.workStatistics();
    const auto host = gpu.hostPlanStatistics();
    std::cout << label << " phases_ms near/P2M/up/M2L/down/L2P=" << work.nearMilliseconds << '/'
              << work.emissionMilliseconds << '/' << work.upwardMilliseconds << '/'
              << work.translationMilliseconds << '/' << work.downwardMilliseconds << '/'
              << work.testingMilliseconds
              << ", cpuPlanPrepareWall_ms=" << work.cpuPlanPreparationMilliseconds
              << ", classBuild_ms=" << work.classPreparationMilliseconds
              << ", azimuthUploads=" << work.azimuthUploads
              << ", azimuth_ms=" << work.azimuthUploadMilliseconds
              << ", rcr_batches/checks=" << work.rcrBatchesQueued << '/' << work.rcrCompletionChecks
              << ", rcr_max_chunk_batches/work=" << work.maximumRcrChunkBatches << '/'
              << work.maximumRcrChunkWork
              << ", hostPlanHits/builds/builderElapsedSum_ms=" << host.hits << '/'
              << host.plansBuilt << '/' << host.builderMilliseconds;
    for (unsigned kind = 0; kind < 4; ++kind)
        std::cout << ", cacheKind" << kind << " hits/misses=" << work.classCacheHits[kind] << '/'
                  << work.classCacheMisses[kind];
    std::cout << ", live/peak=" << gpu.deviceBytes() << '/' << gpu.peakDeviceBytes() << std::endl;
}
void profileSphereOperator() {
    using Clock = std::chrono::steady_clock;
    auto start = Clock::now();
    const auto elapsed = [&] {
        return std::chrono::duration<double, std::milli>(Clock::now() - start).count();
    };
    const auto space = galerkin::buildSpace(testing::sphere(3));
    const double k = 2 * std::numbers::pi * 1000 / 343;
    SupportTreeSettings treeSettings;
    treeSettings.leafCapacity = 4;
    treeSettings.nearDistanceRatio = 1;
    auto tree = std::make_shared<SupportTree>(makeSupportTree(space, treeSettings));
    const std::array<std::uint32_t, 1> faces{1};
    auto nearPlan = std::make_shared<SupportNearPlan>(
        makeSupportNearPlan(space, *tree, faces, 256 * 1024 * 1024));
    std::cout << "PROFILE sphere512T tree+CSR_ms=" << elapsed() << std::endl;
    start = Clock::now();
    SupportNearOperatorSettings nearSettings;
    nearSettings.quadrature.nearDistanceRatio = 1;
    auto near = std::make_shared<SupportNearOperator>(space, *tree, *nearPlan, k, nearSettings);
    std::cout << "PROFILE sphere512T near_ms=" << elapsed() << std::endl;
    start = Clock::now();
    auto pc = std::make_shared<SupportBlockJacobi>(near->pressureMatrix(), tree->pressureOrder, 32,
                                                   std::stop_token{}, near->stream());
    std::cout << "PROFILE sphere512T PC_ms=" << elapsed() << std::endl;
    start = Clock::now();
    SupportFmmOperatorSettings settings;
    settings.panelRule = galerkin::triangleRule(4);
    SupportFmmOperator gpu(space, tree, nearPlan, near, pc, k, settings);
    std::cout << "PROFILE sphere512T FMMprepare_ms=" << elapsed()
              << ", classBatch=" << settings.classBatchNodes << std::endl;
    for (const auto &level : gpu.levelEvidence())
        std::cout << "PROFILE level=" << level.level
                  << ", representation=" << static_cast<unsigned>(level.representation)
                  << ", p=" << level.momentBandlimit << std::endl;
    TestDeviceArray pressure(space.vertices.size()), flux(1), output(space.vertices.size());
    std::vector<Complex> p;
    for (unsigned i = 0; i < space.vertices.size(); ++i)
        p.emplace_back(std::sin(.31 * i), std::cos(.47 * i));
    pressure.upload(p, gpu.stream());
    flux.upload({{.7, -.2}}, gpu.stream());
    for (unsigned call = 0; call < 3; ++call) {
        start = Clock::now();
        gpu.applyPressureDevice(pressure.values, output.values);
        std::cout << "PROFILE sphere512T pressure call=" << call << ", wall_ms=" << elapsed()
                  << '\n';
        printWork(gpu, "PROFILE");
    }
    start = Clock::now();
    gpu.applyFaceFluxDevice(flux.values, output.values);
    std::cout << "PROFILE sphere512T faceFlux wall_ms=" << elapsed() << '\n';
    printWork(gpu, "PROFILE");
}

} // namespace

int main(int argc, char **argv) {
    int devices = 0;
    if (cudaGetDeviceCount(&devices) != cudaSuccess || devices == 0) {
        std::cerr << "CUDA hardware/runtime unavailable; wideband GPU test did not execute.\n";
        return 77;
    }
    try {
        if (argc == 2 && std::string_view(argv[1]) == "--profile-sphere") {
            profileSphereOperator();
            return 0;
        }
        if (argc == 2 && std::string_view(argv[1]) == "--profile-complete") {
            checkCompleteSupportOperator();
            return 0;
        }
        if (argc == 2 && std::string_view(argv[1]) == "--profile-mixed") {
            checkMixedCompleteOperator();
            return 0;
        }
        if (argc == 2 && std::string_view(argv[1]) == "--graph") {
            checkQueuedAzimuthUpdates();
            checkGraphStickyReplay();
            checkM2lGraphLifecycle();
            checkRcrCancellation();
            return 0;
        }
        if (argc == 2 && std::string_view(argv[1]) == "--contractions") {
            checkQueuedAzimuthUpdates();
            checkCooperativeContractions();
            checkScaledAxial();
            checkSphericalRotations();
            checkActiveLfBatches();
            return 0;
        }
        checkDiagonalLevel();
        checkAngularTransfers();
        checkAngularCenterShifts();
        checkDeviceHierarchy();
        checkMixedDeviceHierarchy();
        checkScaledAxial();
        checkCooperativeContractions();
        checkQueuedAzimuthUpdates();
        checkSphericalRotations();
        checkResidentRcr();
        checkRcrQueuePolicy();
        checkActiveLfBatches();
        checkGraphStickyReplay();
        checkM2lGraphLifecycle();
        checkScaledRegularCuda();
        checkIntegratedLeafMoments();
        checkSupportLeafPipeline();
        checkCompleteSupportOperator();
        checkRcrCancellation();
        checkMixedCompleteOperator();
        checkBoundedCpuPlans();
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
