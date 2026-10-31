#include "axisymmetric_assembly.hpp"
#include <algorithm>
#include <cmath>
#include <cuda_runtime.h>
#include <sstream>
#include <stdexcept>

namespace hornsculpt::cuda_detail {
namespace {
using namespace galerkin;
using Complex = cuDoubleComplex;
void check(cudaError_t result, const char *operation) {
    if (result != cudaSuccess)
        throw std::runtime_error(std::string(operation) + ": " + cudaGetErrorString(result));
}
void cancelled(std::stop_token stop) {
    if (stop.stop_requested())
        throw std::runtime_error("Meridian CUDA calculation cancelled.");
}
template <class T> class Buffer {
  public:
    explicit Buffer(std::size_t count) : count_(count) {
        if (count)
            check(cudaMalloc(reinterpret_cast<void **>(&pointer_), count * sizeof(T)),
                  "Ring allocation");
    }
    ~Buffer() {
        if (pointer_)
            cudaFree(pointer_);
    }
    Buffer(const Buffer &) = delete;
    Buffer &operator=(const Buffer &) = delete;
    T *data() const {
        return pointer_;
    }
    std::size_t bytes() const {
        return count_ * sizeof(T);
    }
    void upload(const std::vector<T> &values) {
        if (values.size() != count_)
            throw std::logic_error("Ring upload size mismatch.");
        if (count_)
            check(cudaMemcpy(pointer_, values.data(), bytes(), cudaMemcpyHostToDevice),
                  "Ring upload");
    }
    std::vector<T> download() const {
        std::vector<T> values(count_);
        if (count_)
            check(cudaMemcpy(values.data(), pointer_, bytes(), cudaMemcpyDeviceToHost),
                  "Ring download");
        return values;
    }
    void clear() {
        if (count_)
            check(cudaMemset(pointer_, 0, bytes()), "Ring clear");
    }

  private:
    std::size_t count_;
    T *pointer_{};
};
__device__ RingComplex multiply(RingComplex a, RingComplex b) {
    return {a.real * b.real - a.imag * b.imag, a.real * b.imag + a.imag * b.real};
}
__device__ RingComplex promote(Complex z) {
    return {z.x, z.y};
}
__device__ void atomicValue(Complex *target, RingComplex value, int *error) {
    if (!isfinite(value.real) || !isfinite(value.imag)) {
        atomicCAS(error, 0, 2);
        return;
    }
    atomicAdd(&target->x, value.real);
    atomicAdd(&target->y, value.imag);
}
// Shared storage is scalar POD: no shared-memory constructors or point clouds.
struct PairWorkspace {
    double partial[2][48], evaluation[48], low[48], reference[48], integral[48];
    double u0[52], u1[52], v0[52], v1[52];
    unsigned mappings[52], depths[52];
    double currentU0, currentU1, currentV0, currentV1;
    unsigned currentMapping, currentDepth, pending, work, failure, done;
    double referenceNorm, parameterArea, angularTolerance, localFactor, error, angularError,
        meridianError;
};
__device__ void store(double *out, const RingPairMoments &value) {
    for (unsigned i = 0; i < 16; ++i) {
        out[3 * i] = value.value[i].real;
        out[3 * i + 1] = value.value[i].imag;
        out[3 * i + 2] = value.angularError[i];
    }
}
__device__ RingPairMoments load(const double *in) {
    RingPairMoments value;
    for (unsigned i = 0; i < 16; ++i) {
        value.value[i] = {in[3 * i], in[3 * i + 1]};
        value.angularError[i] = in[3 * i + 2];
    }
    return value;
}
__device__ void add(double *target, const double *value) {
    for (unsigned i = 0; i < 48; ++i)
        target[i] += value[i];
}
__device__ void evaluate(PairWorkspace &w, const RingSegment &a, const RingSegment &b,
                         RingPairPatch patch, const RingRulePoint *rule, unsigned order, double k,
                         double angularTolerance) {
    RingPairMoments moment{};
    for (unsigned index = threadIdx.x; index < order * order; index += blockDim.x) {
        const auto point = mapRingPairPoint(a, b, patch, rule[index / order], rule[index % order]);
        const unsigned status = accumulateRingPairPoint(moment, point, k, angularTolerance);
        if (status)
            atomicCAS(&w.failure, 0u, status);
    }
    for (unsigned offset = 16; offset; offset /= 2)
        for (unsigned i = 0; i < 16; ++i) {
            moment.value[i].real += __shfl_down_sync(0xffffffff, moment.value[i].real, offset);
            moment.value[i].imag += __shfl_down_sync(0xffffffff, moment.value[i].imag, offset);
            moment.angularError[i] += __shfl_down_sync(0xffffffff, moment.angularError[i], offset);
        }
    if ((threadIdx.x & 31) == 0)
        store(w.partial[threadIdx.x / 32], moment);
    __syncthreads();
    if (threadIdx.x == 0) {
        for (unsigned i = 0; i < 48; ++i)
            w.evaluation[i] = w.partial[0][i] + w.partial[1][i];
    }
    __syncthreads();
}
__device__ void push(PairWorkspace &w, RingPairPatch p) {
    const unsigned i = w.pending++;
    w.u0[i] = p.firstU;
    w.u1[i] = p.lastU;
    w.v0[i] = p.firstV;
    w.v1[i] = p.lastV;
    w.mappings[i] = p.mapping;
    w.depths[i] = p.depth;
}
__device__ RingPairPatch current(const PairWorkspace &w) {
    return {w.currentU0, w.currentU1, w.currentV0, w.currentV1, w.currentMapping, w.currentDepth};
}
__device__ void integrate(PairWorkspace &w, const RingSegment &a, const RingSegment &b,
                          const RingPairPatch *roots, RingPairRange range,
                          const RingRulePoint *regular, unsigned regularOrder,
                          const RingRulePoint *coarse, unsigned coarseOrder,
                          const RingRulePoint *fine, unsigned fineOrder, double k, double tolerance,
                          unsigned maximumDepth) {
    if (threadIdx.x == 0) {
        w.failure = w.done = 0;
        w.parameterArea = 0;
        w.localFactor = 1;
        w.angularTolerance = tolerance * .005 < 1e-8 ? tolerance * .005 : 1e-8;
        for (unsigned i = 0; i < 48; ++i)
            w.reference[i] = w.integral[i] = 0;
    }
    __syncthreads();
    if (!range.count) {
        evaluate(w, a, b, {0, 1, 0, 1, 0, 0}, regular, regularOrder, k, w.angularTolerance);
        if (threadIdx.x == 0)
            for (unsigned i = 0; i < 48; ++i)
                w.integral[i] = w.evaluation[i];
        __syncthreads();
        return;
    }
    for (unsigned root = 0; root < range.count; ++root) {
        const auto patch = roots[range.first + root];
        evaluate(w, a, b, patch, fine, fineOrder, k, w.angularTolerance);
        if (threadIdx.x == 0) {
            add(w.reference, w.evaluation);
            w.parameterArea += (patch.lastU - patch.firstU) * (patch.lastV - patch.firstV);
        }
        __syncthreads();
    }
    if (threadIdx.x == 0) {
        w.referenceNorm = ringOperatorNorm(load(w.reference), a, b, k);
        if (!(w.referenceNorm > 0))
            w.failure = 2;
    }
    __syncthreads();
    for (unsigned attempt = 0; attempt < 3 && !w.failure; ++attempt) {
        if (threadIdx.x == 0) {
            w.pending = w.work = 0;
            w.error = w.angularError = w.meridianError = 0;
            for (unsigned i = 0; i < 48; ++i)
                w.integral[i] = 0;
            for (unsigned root = 0; root < range.count; ++root)
                push(w, roots[range.first + root]);
        }
        __syncthreads();
        while (w.pending && !w.failure) {
            if (threadIdx.x == 0) {
                const unsigned i = --w.pending;
                w.currentU0 = w.u0[i];
                w.currentU1 = w.u1[i];
                w.currentV0 = w.v0[i];
                w.currentV1 = w.v1[i];
                w.currentMapping = w.mappings[i];
                w.currentDepth = w.depths[i];
            }
            __syncthreads();
            const auto patch = current(w);
            evaluate(w, a, b, patch, coarse, coarseOrder, k, w.angularTolerance);
            if (threadIdx.x == 0)
                for (unsigned i = 0; i < 48; ++i)
                    w.low[i] = w.evaluation[i];
            __syncthreads();
            evaluate(w, a, b, patch, fine, fineOrder, k, w.angularTolerance);
            if (threadIdx.x == 0 && !w.failure) {
                const auto difference = ringSubtract(load(w.evaluation), load(w.low));
                const double spatial = ringOperatorNorm(difference, a, b, k);
                const double angular = ringOperatorNorm(difference, a, b, k, true);
                const double fraction =
                    (patch.lastU - patch.firstU) * (patch.lastV - patch.firstV) / w.parameterArea;
                const double budget = tolerance * .1 * w.referenceNorm * fraction * w.localFactor;
                if (spatial <= budget || patch.depth == maximumDepth) {
                    add(w.integral, w.evaluation);
                    w.error += spatial + angular;
                    w.angularError += angular;
                    w.meridianError += spatial;
                } else if (++w.work > 50000 || w.pending + 4 > 52)
                    w.failure = 4;
                else {
                    const double u = .5 * (patch.firstU + patch.lastU),
                                 v = .5 * (patch.firstV + patch.lastV);
                    for (unsigned i = 0; i < 4; ++i)
                        push(w, {i & 1 ? u : patch.firstU, i & 1 ? patch.lastU : u,
                                 i & 2 ? v : patch.firstV, i & 2 ? patch.lastV : v, patch.mapping,
                                 patch.depth + 1});
                }
            }
            __syncthreads();
        }
        if (threadIdx.x == 0 && !w.failure) {
            const double norm = ringOperatorNorm(load(w.integral), a, b, k);
            if (norm > 0 && w.error <= tolerance * .2 * norm)
                w.done = 1;
            else if (attempt == 2)
                w.failure = 5;
            else if (w.angularError > w.meridianError)
                w.angularTolerance *= .1;
            else
                w.localFactor *= .25;
        }
        __syncthreads();
        if (w.done)
            break;
    }
}
__global__ void assemblePairs(const RingSegment *segments, unsigned segmentCount,
                              const RingPairRange *ranges, const RingPairPatch *roots,
                              const RingRulePoint *regular, unsigned regularOrder,
                              const RingRulePoint *coarse, unsigned coarseOrder,
                              const RingRulePoint *fine, unsigned fineOrder, const int *apertures,
                              const unsigned *representatives, unsigned pressureCount,
                              unsigned dimension, unsigned panelCount, unsigned rhsCount, double k,
                              double tolerance, unsigned maximumDepth, const double *scales,
                              const Complex *flux, Complex *matrix, Complex *rhs, int *error,
                              unsigned *failedPair, unsigned first) {
    const unsigned pair = first + blockIdx.x, ti = pair / segmentCount, si = pair % segmentCount;
    const auto a = segments[ti], b = segments[si];
    __shared__ PairWorkspace w;
    integrate(w, a, b, roots, ranges[pair], regular, regularOrder, coarse, coarseOrder, fine,
              fineOrder, k, tolerance, maximumDepth);
    if (threadIdx.x || w.failure) {
        if (threadIdx.x == 0 && w.failure && atomicCAS(error, 0, static_cast<int>(w.failure)) == 0)
            *failedPair = pair;
        return;
    }
    const auto m = load(w.integral);
    const unsigned rows[2]{a.first, a.second}, columns[2]{b.first, b.second};
    RingComplex cosine{}, green{};
    for (unsigned q = 0; q < 4; ++q) {
        cosine = ringAdd(cosine, m.value[4 + q]);
        green = ringAdd(green, m.value[q]);
    }
    for (unsigned i = 0; i < 2; ++i) {
        RingComplex single{}, adjoint{};
        for (unsigned j = 0; j < 2; ++j) {
            const unsigned q = 2 * i + j;
            single = ringAdd(single, m.value[q]);
            adjoint = ringAdd(adjoint, m.value[8 + q]);
            const double ci = (i ? -1.0 : 1.0) / a.length, cj = (j ? -1.0 : 1.0) / b.length;
            const auto h = ringAdd(ringScale(cosine, ci * cj),
                                   ringAdd(ringScale(m.value[4 + q], -k * k * a.nr * b.nr),
                                           ringScale(m.value[q], -k * k * a.nz * b.nz)));
            const auto value = ringAdd(m.value[12 + q], {h.imag / k, -h.real / k});
            atomicValue(matrix + rows[i] + static_cast<std::size_t>(columns[j]) * dimension,
                        ringScale(value, scales[rows[i]]), error);
        }
        const RingComplex coefficient{single.real - adjoint.imag / k,
                                      single.imag + adjoint.real / k};
        if (apertures[si] >= 0)
            atomicValue(matrix + rows[i] +
                            static_cast<std::size_t>(pressureCount + apertures[si]) * dimension,
                        ringScale(coefficient, -k * scales[rows[i]]), error);
        else
            for (unsigned r = 0; r < rhsCount; ++r)
                atomicValue(
                    rhs + rows[i] + static_cast<std::size_t>(r) * dimension,
                    ringScale(multiply(coefficient,
                                       promote(flux[representatives[si] +
                                                    static_cast<std::size_t>(r) * panelCount])),
                              scales[rows[i]]),
                    error);
    }
    if (apertures[ti] >= 0 && apertures[si] >= 0) {
        const unsigned row = pressureCount + apertures[ti], column = pressureCount + apertures[si];
        atomicValue(matrix + row + static_cast<std::size_t>(column) * dimension,
                    ringScale(green, 2 * k * scales[row]), error);
    }
}
__global__ void addMass(const RingSegment *segments, unsigned count, const int *apertures,
                        const unsigned *representatives, unsigned pressureCount, unsigned dimension,
                        unsigned panelCount, unsigned rhsCount, double k, const double *scales,
                        const Complex *flux, Complex *matrix, Complex *rhs, int *error) {
    const unsigned t = blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= count)
        return;
    const auto e = segments[t];
    const double r0 = e.radius, r1 = e.radius + e.dr;
    const double f = 2 * ring_detail::pi * e.length;
    const double mass[4]{f * (3 * r0 + r1) / 12, f * (r0 + r1) / 12, f * (r0 + r1) / 12,
                         f * (r0 + 3 * r1) / 12};
    const double mixed[2]{f * (2 * r0 + r1) / 6, f * (r0 + 2 * r1) / 6};
    const unsigned nodes[2]{e.first, e.second};
    for (unsigned i = 0; i < 2; ++i) {
        const unsigned row = nodes[i];
        for (unsigned j = 0; j < 2; ++j)
            atomicValue(matrix + row + static_cast<std::size_t>(nodes[j]) * dimension,
                        {.5 * mass[2 * i + j] * scales[row], 0}, error);
        if (apertures[t] >= 0)
            atomicValue(matrix + row +
                            static_cast<std::size_t>(pressureCount + apertures[t]) * dimension,
                        {0, .5 * mixed[i] * scales[row]}, error);
        else
            for (unsigned r = 0; r < rhsCount; ++r)
                atomicValue(
                    rhs + row + static_cast<std::size_t>(r) * dimension,
                    multiply(
                        {0, -.5 * mixed[i] * scales[row] / k},
                        promote(
                            flux[representatives[t] + static_cast<std::size_t>(r) * panelCount])),
                    error);
    }
    if (apertures[t] >= 0) {
        const unsigned row = pressureCount + apertures[t];
        for (unsigned i = 0; i < 2; ++i)
            atomicValue(matrix + row + static_cast<std::size_t>(nodes[i]) * dimension,
                        {mixed[i] * scales[row], 0}, error);
    }
}
struct FieldIntegral {
    RingComplex value{};
    double magnitude{}, angularError{};
    unsigned status{};
};
__device__ FieldIntegral fieldInterval(const RingSegment &segment, const RingObservation &observer,
                                       double first, double last, const RingRulePoint *rule,
                                       unsigned order, double k, double angularTolerance,
                                       bool sinhMap) {
    FieldIntegral out;
    for (unsigned i = 0; i < order; ++i) {
        const double u = first + (last - first) * rule[i].position;
        auto mapped = ring_pair_detail::SinhPoint{u, 1, 0};
        if (sinhMap)
            mapped = ring_pair_detail::projectedSinh(observer.radius - segment.radius,
                                                     observer.axial - segment.axial, segment.dr,
                                                     segment.dz, segment.length, u);
        if (mapped.status) {
            out.status = mapped.status;
            return out;
        }
        const double r = segment.radius + mapped.position * segment.dr;
        RingGeometry g{observer.radius,
                       r,
                       observer.radius - r,
                       observer.axial - segment.axial - mapped.position * segment.dz,
                       0,
                       0,
                       0,
                       0,
                       0,
                       0};
        const auto kernel = ringKernels(g, k, angularTolerance);
        if (kernel.status) {
            out.status = kernel.status;
            return out;
        }
        const double weight =
            (last - first) * rule[i].weight * mapped.jacobian * segment.length * r;
        out.value = ringAdd(out.value, ringScale(kernel.value[0], weight));
        out.magnitude += weight * ringMagnitude(kernel.value[0]);
        out.angularError += weight * kernel.angularAbsoluteError[0];
    }
    return out;
}
__device__ FieldIntegral finiteField(const RingSegment &segment, const RingObservation &observer,
                                     const RingRulePoint *coarse, unsigned coarseCount,
                                     const RingRulePoint *fine, unsigned fineCount, double k,
                                     double tolerance, unsigned maxDepth) {
    const double center =
        ring_pair_detail::clamped(((observer.radius - segment.radius) * segment.dr +
                                   (observer.axial - segment.axial) * segment.dz) /
                                  (segment.length * segment.length));
    const double distance = hypot(observer.radius - segment.radius - center * segment.dr,
                                  observer.axial - segment.axial - center * segment.dz);
    const bool sinhMap = distance < segment.length;
    double angularTolerance = tolerance * .005 < 1e-8 ? tolerance * .005 : 1e-8;
    const auto reference =
        fieldInterval(segment, observer, 0, 1, fine, fineCount, k, angularTolerance, sinhMap);
    if (reference.status)
        return reference;
    double localFactor = 1;
    for (unsigned attempt = 0; attempt < 3; ++attempt) {
        struct Patch {
            double first, last;
            unsigned depth;
        };
        Patch stack[18];
        unsigned pending = 1, work = 0;
        stack[0] = {0, 1, 0};
        FieldIntegral result;
        double spatialError = 0;
        while (pending) {
            const auto p = stack[--pending];
            const auto low = fieldInterval(segment, observer, p.first, p.last, coarse, coarseCount,
                                           k, angularTolerance, sinhMap);
            const auto high = fieldInterval(segment, observer, p.first, p.last, fine, fineCount, k,
                                            angularTolerance, sinhMap);
            if (low.status || high.status) {
                result.status = low.status ? low.status : high.status;
                return result;
            }
            const double error =
                ringMagnitude({high.value.real - low.value.real, high.value.imag - low.value.imag});
            if (error <= tolerance * .1 * reference.magnitude * (p.last - p.first) * localFactor ||
                p.depth == maxDepth) {
                result.value = ringAdd(result.value, high.value);
                result.magnitude += high.magnitude;
                result.angularError += high.angularError + low.angularError;
                spatialError += error;
            } else {
                if (pending + 2 > 18 || ++work > 50000) {
                    result.status = 4;
                    return result;
                }
                const double mid = .5 * (p.first + p.last);
                stack[pending++] = {mid, p.last, p.depth + 1};
                stack[pending++] = {p.first, mid, p.depth + 1};
            }
        }
        if (spatialError + result.angularError <= tolerance * .2 * result.magnitude)
            return result;
        if (result.angularError > spatialError)
            angularTolerance *= .1;
        else
            localFactor *= .25;
    }
    FieldIntegral failed;
    failed.status = 5;
    return failed;
}
__global__ void fields(const RingSegment *segments, unsigned count, const int *apertures,
                       const RingObservation *observers, unsigned observerCount, bool far,
                       double referenceDistance, const RingRulePoint *coarse, unsigned coarseCount,
                       const RingRulePoint *fine, unsigned fineCount, unsigned pressureCount,
                       unsigned dimension, unsigned rhsCount, double k, double tolerance,
                       unsigned maxDepth, const Complex *solution, Complex *output, int *error) {
    const unsigned o = blockIdx.x, rhs = blockIdx.y;
    const auto observer = observers[o];
    RingComplex value{};
    if (!observer.blocked)
        for (unsigned t = threadIdx.x; t < count; t += blockDim.x) {
            if (apertures[t] < 0)
                continue;
            const auto segment = segments[t];
            const auto q = ringScale(promote(solution[pressureCount + apertures[t] +
                                                      static_cast<std::size_t>(rhs) * dimension]),
                                     k);
            if (q.real == 0 && q.imag == 0)
                continue;
            RingComplex integral{};
            if (far) {
                for (unsigned j = 0; j < fineCount; ++j) {
                    const double u = fine[j].position, radial = segment.radius + u * segment.dr;
                    const double phase =
                        k * (observer.originPhase +
                             observer.directionZ * (segment.axial + u * segment.dz) -
                             referenceDistance);
                    double sine, cosine;
                    sincos(phase, &sine, &cosine);
                    const double factor = -fine[j].weight * segment.length * radial *
                                          j0(k * radial * observer.directionR) / referenceDistance;
                    integral = ringAdd(integral, {factor * cosine, factor * sine});
                }
            } else {
                const auto integrated = finiteField(segment, observer, coarse, coarseCount, fine,
                                                    fineCount, k, tolerance, maxDepth);
                if (integrated.status) {
                    atomicCAS(error, 0, static_cast<int>(integrated.status));
                    continue;
                }
                integral = ringScale(integrated.value, -2);
            }
            value = ringAdd(value, multiply(integral, q));
        }
    for (unsigned offset = 16; offset; offset /= 2) {
        value.real += __shfl_down_sync(0xffffffff, value.real, offset);
        value.imag += __shfl_down_sync(0xffffffff, value.imag, offset);
    }
    __shared__ double partial[2][2];
    if ((threadIdx.x & 31) == 0) {
        partial[threadIdx.x / 32][0] = value.real;
        partial[threadIdx.x / 32][1] = value.imag;
    }
    __syncthreads();
    if (threadIdx.x == 0)
        output[o + static_cast<std::size_t>(rhs) * observerCount] = {partial[0][0] + partial[1][0],
                                                                     partial[0][1] + partial[1][1]};
}
} // namespace

class AxisymmetricAssembly::Implementation {
  public:
    Implementation(const galerkin::RingInteractionPlan &plan, const std::vector<int> &aperture,
                   const std::vector<unsigned> &representatives, const std::vector<unsigned> &faces)
        : segments(plan.segments.size()), ranges(plan.pairs.size()), roots(plan.roots.size()),
          regular(plan.regular.size()), coarse(plan.coarse.size()), fine(plan.fine.size()),
          apertures(aperture.size()), panels(representatives.size()), error(1), failedPair(1),
          count(static_cast<unsigned>(plan.segments.size())), regularOrder(plan.regularOrder),
          coarseOrder(plan.coarseOrder), fineOrder(plan.fineOrder), faceIds(faces) {
        if (aperture.size() != count || representatives.size() != count || faces.size() != count)
            throw std::invalid_argument("Continuous ring metadata size mismatch.");
        segments.upload(plan.segments);
        ranges.upload(plan.pairs);
        roots.upload(plan.roots);
        regular.upload(plan.regular);
        coarse.upload(plan.coarse);
        fine.upload(plan.fine);
        apertures.upload(aperture);
        panels.upload(representatives);
    }
    void validate(const char *context) {
        const int status = error.download()[0];
        if (!status)
            return;
        const auto pair = failedPair.download()[0], a = pair / count, b = pair % count;
        std::ostringstream text;
        text << context << ": status=" << status << ", intervals " << a << '/' << b;
        if (a < faceIds.size() && b < faceIds.size())
            text << " (faces " << faceIds[a] << '/' << faceIds[b] << ')';
        throw std::runtime_error(text.str());
    }
    std::size_t bytes() const {
        return segments.bytes() + ranges.bytes() + roots.bytes() + regular.bytes() +
               coarse.bytes() + fine.bytes() + apertures.bytes() + panels.bytes() + error.bytes() +
               failedPair.bytes();
    }
    Buffer<galerkin::RingSegment> segments;
    Buffer<galerkin::RingPairRange> ranges;
    Buffer<galerkin::RingPairPatch> roots;
    Buffer<galerkin::RingRulePoint> regular, coarse, fine;
    Buffer<int> apertures;
    Buffer<unsigned> panels;
    Buffer<int> error;
    Buffer<unsigned> failedPair;
    unsigned count, regularOrder, coarseOrder, fineOrder;
    std::vector<unsigned> faceIds;
};
AxisymmetricAssembly::AxisymmetricAssembly(const galerkin::RingInteractionPlan &plan,
                                           const std::vector<int> &apertures,
                                           const std::vector<unsigned> &representatives,
                                           const std::vector<unsigned> &faces)
    : implementation_(std::make_unique<Implementation>(plan, apertures, representatives, faces)) {}
AxisymmetricAssembly::~AxisymmetricAssembly() = default;
std::size_t AxisymmetricAssembly::bytes() const {
    return implementation_->bytes();
}
void AxisymmetricAssembly::assemble(unsigned pressureCount, unsigned dimension, unsigned panelCount,
                                    unsigned rhsCount, double k, double tolerance,
                                    unsigned maximumDepth, const double *rowScales,
                                    const cuDoubleComplex *knownFlux, cuDoubleComplex *matrix,
                                    cuDoubleComplex *rhs, std::stop_token stop) {
    auto &s = *implementation_;
    s.error.clear();
    s.failedPair.clear();
    const std::size_t count = static_cast<std::size_t>(s.count) * s.count;
    for (std::size_t first = 0; first < count; first += 4096) {
        cancelled(stop);
        const unsigned batch = static_cast<unsigned>(std::min<std::size_t>(4096, count - first));
        assemblePairs<<<batch, 64>>>(
            s.segments.data(), s.count, s.ranges.data(), s.roots.data(), s.regular.data(),
            s.regularOrder, s.coarse.data(), s.coarseOrder, s.fine.data(), s.fineOrder,
            s.apertures.data(), s.panels.data(), pressureCount, dimension, panelCount, rhsCount, k,
            tolerance, maximumDepth, rowScales, knownFlux, matrix, rhs, s.error.data(),
            s.failedPair.data(), static_cast<unsigned>(first));
        check(cudaDeviceSynchronize(), "Continuous meridian pair assembly");
        s.validate("Meridian quadrature failed");
    }
    addMass<<<(s.count + 127) / 128, 128>>>(
        s.segments.data(), s.count, s.apertures.data(), s.panels.data(), pressureCount, dimension,
        panelCount, rhsCount, k, rowScales, knownFlux, matrix, rhs, s.error.data());
    check(cudaDeviceSynchronize(), "Continuous meridian mass assembly");
    s.validate("Meridian mass failed");
}
void AxisymmetricAssembly::field(const std::vector<RingObservation> &observations, bool far,
                                 double referenceDistance, unsigned pressureCount,
                                 unsigned dimension, unsigned rhsCount, double k, double tolerance,
                                 unsigned maxDepth, const cuDoubleComplex *solution,
                                 cuDoubleComplex *output, std::stop_token stop) {
    auto &s = *implementation_;
    cancelled(stop);
    s.error.clear();
    s.failedPair.clear();
    Buffer<RingObservation> observers(observations.size());
    observers.upload(observations);
    fields<<<dim3(static_cast<unsigned>(observations.size()), rhsCount), 64>>>(
        s.segments.data(), s.count, s.apertures.data(), observers.data(),
        static_cast<unsigned>(observations.size()), far, referenceDistance, s.coarse.data(),
        s.coarseOrder, s.fine.data(), s.fineOrder, pressureCount, dimension, rhsCount, k, tolerance,
        maxDepth, solution, output, s.error.data());
    check(cudaDeviceSynchronize(), "Continuous meridian Rayleigh field");
    s.validate("Meridian field quadrature failed");
}
} // namespace hornsculpt::cuda_detail
