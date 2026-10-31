#include "../axisymmetric_reduction.hpp"
#include "../baffle_exterior.hpp"
#include "../galerkin/interaction_plan.hpp"
#include "../galerkin/space.hpp"
#include "../lem_coupling.hpp"
#include "../planning_frequencies.hpp"
#include "axisymmetric_assembly.hpp"
#include "baffled_solver.hpp"
#include "panel_integrals.cuh"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cublas_v2.h>
#include <cusolverDn.h>
#include <hornsculpt/meshing.hpp>
#include <hornsculpt/field_evaluation.hpp>
#include <hornsculpt/surface_measures.hpp>
#include <limits>
#include <map>
#include <numbers>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>

namespace hornsculpt::cuda_detail {
namespace {
using Complex = cuDoubleComplex;
using Clock = std::chrono::steady_clock;
void check(cudaError_t status, const char *operation) {
    if (status != cudaSuccess)
        throw std::runtime_error(std::string(operation) + ": " + cudaGetErrorString(status));
}
void check(cublasStatus_t status, const char *operation) {
    if (status != CUBLAS_STATUS_SUCCESS)
        throw std::runtime_error(std::string(operation) + ": cuBLAS " + std::to_string(status));
}
void check(cusolverStatus_t status, const char *operation) {
    if (status != CUSOLVER_STATUS_SUCCESS)
        throw std::runtime_error(std::string(operation) + ": cuSOLVER " + std::to_string(status));
}
void cancelled(std::stop_token stop) {
    if (stop.stop_requested())
        throw std::runtime_error("Baffled BEM calculation cancelled.");
}
double elapsed(Clock::time_point start) {
    return std::chrono::duration<double, std::milli>(Clock::now() - start).count();
}
bool finite(Vec3 p) {
    return std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z);
}
template <class T> class Device {
  public:
    explicit Device(std::size_t count) : count_(count) {
        if (count > std::numeric_limits<std::size_t>::max() / sizeof(T))
            throw std::invalid_argument("Baffled BEM allocation overflow.");
        if (count)
            check(cudaMalloc(reinterpret_cast<void **>(&pointer_), count * sizeof(T)),
                  "Baffled allocation");
    }
    ~Device() {
        if (pointer_)
            cudaFree(pointer_);
    }
    Device(const Device &) = delete;
    Device &operator=(const Device &) = delete;
    T *data() const {
        return pointer_;
    }
    std::size_t bytes() const {
        return count_ * sizeof(T);
    }
    void upload(const std::vector<T> &value) {
        if (value.size() != count_)
            throw std::logic_error("Baffled upload size mismatch.");
        if (count_)
            check(cudaMemcpy(pointer_, value.data(), bytes(), cudaMemcpyHostToDevice),
                  "Baffled upload");
    }
    std::vector<T> download() const {
        std::vector<T> value(count_);
        if (count_)
            check(cudaMemcpy(value.data(), pointer_, bytes(), cudaMemcpyDeviceToHost),
                  "Baffled download");
        return value;
    }
    void clear() {
        if (count_)
            check(cudaMemset(pointer_, 0, bytes()), "Baffled clear");
    }

  private:
    std::size_t count_;
    T *pointer_{};
};
class Handles {
  public:
    Handles() {
        check(cublasCreate(&blas), "Baffled cuBLAS create");
        const auto status = cusolverDnCreate(&solver);
        if (status != CUSOLVER_STATUS_SUCCESS) {
            cublasDestroy(blas);
            check(status, "Baffled cuSOLVER create");
        }
    }
    ~Handles() {
        cusolverDnDestroy(solver);
        cublasDestroy(blas);
    }
    cublasHandle_t blas{};
    cusolverDnHandle_t solver{};
};
__device__ void atomicValue(Complex *output, DoubleComplex value, int *error) {
    if (!isfinite(value.real) || !isfinite(value.imag)) {
        atomicExch(error, 1);
        return;
    }
    atomicAdd(&output->x, value.real);
    atomicAdd(&output->y, value.imag);
}
// Same regularized H and Sauter--Schwab rules as the exterior solver; only
// interior trace signs differ. nu points OUT OF AIR, alpha=-i/k.
__device__ void interiorSample(PairMoments &result, const Panel &test, const Panel &trial,
                               galerkin::PairQuadraturePoint sample, Vec3 delta, double k,
                               bool coincident, int *error) {
    const double r = length(delta);
    if (!(r > 0) || !isfinite(r)) {
        atomicExch(error, 1);
        return;
    }
    double sine, cosine;
    sincos(k * r, &sine, &cosine);
    const double scale = sample.weight * 4 * test.area * trial.area / (12.56637061435917295385 * r);
    const DoubleComplex green{cosine * scale, -sine * scale};
    result.green = plus(result.green, green);
    const DoubleComplex radial{green.real - green.imag * k * r, green.imag + green.real * k * r};
    const auto d = times(radial, coincident ? 0 : dot(delta, trial.normal) / (r * r));
    const auto adjoint = times(radial, coincident ? 0 : -dot(delta, test.normal) / (r * r));
    const double nn = dot(test.normal, trial.normal);
    const DoubleComplex a{d.real - k * nn * green.imag, d.imag + k * nn * green.real};
    const DoubleComplex b{green.real - adjoint.imag / k, green.imag + adjoint.real / k};
    for (unsigned i = 0; i < 3; ++i) {
        result.rhs[i] = plus(result.rhs[i], times(b, sample.test[i]));
        for (unsigned j = 0; j < 3; ++j)
            result.matrix[3 * i + j] =
                plus(result.matrix[3 * i + j], times(a, sample.test[i] * sample.trial[j]));
    }
}
__global__ void
assemblePairs(const Panel *panels, unsigned panelCount, unsigned pressureCount, unsigned dimension,
              const int *apertures, const double *rowScale, const Complex *flux, unsigned rhsCount,
              double k, const galerkin::TriangleQuadraturePoint *regular, unsigned regularCount,
              const galerkin::GaussPoint *gauss, unsigned singularOrder,
              const galerkin::NearInteraction *near, const unsigned *nearOffsets,
              const unsigned *nearTrials, const galerkin::PairQuadraturePatch *patches,
              const galerkin::TriangleQuadraturePoint *nearRules, unsigned baseCount,
              unsigned highCount, Complex *matrix, Complex *rhs, int *error, unsigned firstTrial) {
    const unsigned ti = blockIdx.x, si = firstTrial + blockIdx.y;
    if (ti >= panelCount || si >= panelCount)
        return;
    const auto &test = panels[ti], &trial = panels[si];
    const auto permutation = galerkin::classifyPairIds(test.nodes, trial.nodes);
    PairMoments moments{};
    if (permutation.topology != galerkin::PairTopology::disjoint) {
        const unsigned perSector = singularOrder * singularOrder * singularOrder * singularOrder;
        const unsigned count = galerkin::pairSectorCount(permutation.topology) * perSector;
        for (unsigned index = threadIdx.x; index < count; index += blockDim.x) {
            unsigned tuple = index % perSector;
            const auto w = tuple % singularOrder;
            tuple /= singularOrder;
            const auto v = tuple % singularOrder;
            tuple /= singularOrder;
            const auto u = tuple % singularOrder;
            tuple /= singularOrder;
            const auto r = tuple;
            auto point = galerkin::canonicalPairTransform(permutation.topology, index / perSector,
                                                          gauss[r].position, gauss[u].position,
                                                          gauss[v].position, gauss[w].position);
            point.weight *= gauss[r].weight * gauss[u].weight * gauss[v].weight * gauss[w].weight;
            point = galerkin::permutePairPoint(point, permutation);
            interiorSample(moments, test, trial, point,
                           pairDelta(test, trial, point.test, point.trial), k,
                           permutation.topology == galerkin::PairTopology::coincident, error);
        }
    } else {
        unsigned lo = nearOffsets[ti], hi = nearOffsets[ti + 1];
        while (lo < hi) {
            const unsigned mid = lo + (hi - lo) / 2;
            if (nearTrials[mid] < si)
                lo = mid + 1;
            else
                hi = mid;
        }
        if (lo < nearOffsets[ti + 1] && nearTrials[lo] == si) {
            const auto &job = near[lo];
            for (unsigned p = 0; p < job.patchCount; ++p) {
                const auto &patch = patches[job.firstPatch + p];
                const bool high = patch.degree == 14;
                const auto *rule = nearRules + (high ? baseCount : 0);
                const unsigned count = high ? highCount : baseCount;
                const unsigned points = galerkin::nearPatchPointCount(patch, count);
                for (unsigned index = threadIdx.x; index < points; index += blockDim.x) {
                    auto point = galerkin::mapNearPatchPoint(rule, count, patch, index,
                                                             job.swapBarycentric ? trial : test,
                                                             job.swapBarycentric ? test : trial);
                    if (job.swapBarycentric)
                        for (unsigned i = 0; i < 3; ++i) {
                            const double swap = point.test[i];
                            point.test[i] = point.trial[i];
                            point.trial[i] = swap;
                        }
                    interiorSample(moments, test, trial, point,
                                   pairDelta(test, trial, point.test, point.trial), k, false,
                                   error);
                }
            }
        } else {
            for (unsigned index = threadIdx.x; index < regularCount * regularCount;
                 index += blockDim.x) {
                const auto &x = regular[index / regularCount], &y = regular[index % regularCount];
                galerkin::PairQuadraturePoint point{};
                for (unsigned i = 0; i < 3; ++i) {
                    point.test[i] = x.barycentric[i];
                    point.trial[i] = y.barycentric[i];
                }
                point.weight = x.weight * y.weight;
                interiorSample(moments, test, trial, point,
                               pairDelta(test, trial, point.test, point.trial), k, false, error);
            }
        }
    }
    reduceWarp(moments);
    __shared__ PairMoments sums[2];
    if ((threadIdx.x & 31) == 0)
        sums[threadIdx.x / 32] = moments;
    __syncthreads();
    if (threadIdx.x != 0)
        return;
    const auto &other = sums[1];
    for (unsigned i = 0; i < 9; ++i)
        moments.matrix[i] = plus(moments.matrix[i], other.matrix[i]);
    for (unsigned i = 0; i < 3; ++i)
        moments.rhs[i] = plus(moments.rhs[i], other.rhs[i]);
    moments.green = plus(moments.green, other.green);
    for (unsigned i = 0; i < 3; ++i) {
        const unsigned row = test.nodes[i];
        const double scale = rowScale[row];
        for (unsigned j = 0; j < 3; ++j) {
            const double curl = dot(test.curls[i], trial.curls[j]) / k;
            const auto value = plus(moments.matrix[3 * i + j],
                                    {moments.green.imag * curl, -moments.green.real * curl});
            atomicValue(matrix + row + static_cast<std::size_t>(trial.nodes[j]) * dimension,
                        times(value, scale), error);
        }
        if (apertures[si] >= 0) {
            const unsigned column = pressureCount + apertures[si];
            // The unknown is q_A/k, balancing pressure and flux columns in SI.
            atomicValue(matrix + row + static_cast<std::size_t>(column) * dimension,
                        times(moments.rhs[i], -k * scale), error);
        } else {
            for (unsigned r = 0; r < rhsCount; ++r)
                atomicValue(
                    rhs + row + static_cast<std::size_t>(r) * dimension,
                    times(product(moments.rhs[i],
                                  promoted(flux[si + static_cast<std::size_t>(r) * panelCount])),
                          scale),
                    error);
        }
    }
    if (apertures[ti] >= 0 && apertures[si] >= 0) {
        const unsigned row = pressureCount + apertures[ti], column = pressureCount + apertures[si];
        atomicValue(matrix + row + static_cast<std::size_t>(column) * dimension,
                    times(moments.green, 2 * k * rowScale[row]), error);
    }
}
__global__ void addInteriorMass(const Panel *panels, unsigned panelCount, unsigned pressureCount,
                                unsigned dimension, const int *apertures, const double *rowScale,
                                const Complex *flux, unsigned rhsCount, double k, Complex *matrix,
                                Complex *rhs, int *error) {
    const unsigned t = blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= panelCount)
        return;
    const auto &panel = panels[t];
    for (unsigned i = 0; i < 3; ++i) {
        const unsigned row = panel.nodes[i];
        const double scale = rowScale[row];
        for (unsigned j = 0; j < 3; ++j)
            atomicValue(matrix + row + static_cast<std::size_t>(panel.nodes[j]) * dimension,
                        {scale * panel.area * (i == j ? 2.0 : 1.0) / 24.0, 0}, error);
        if (apertures[t] >= 0) {
            const unsigned column = pressureCount + apertures[t];
            atomicValue(matrix + row + static_cast<std::size_t>(column) * dimension,
                        {0, scale * panel.area / 6}, error);
        } else {
            const DoubleComplex coefficient{0, -scale * panel.area / (6 * k)};
            for (unsigned r = 0; r < rhsCount; ++r)
                atomicValue(rhs + row + static_cast<std::size_t>(r) * dimension,
                            product(coefficient,
                                    promoted(flux[t + static_cast<std::size_t>(r) * panelCount])),
                            error);
        }
    }
    if (apertures[t] >= 0) {
        const unsigned row = pressureCount + apertures[t];
        for (unsigned j = 0; j < 3; ++j)
            atomicValue(matrix + row + static_cast<std::size_t>(panel.nodes[j]) * dimension,
                        {rowScale[row] * panel.area / 3, 0}, error);
    }
}
__device__ FieldIntegral layerPatch(Vec3 observer, const Panel &panel, const FieldPatch &patch,
                                    double k, const galerkin::TriangleQuadraturePoint *rule,
                                    unsigned ruleCount, int *error) {
    FieldIntegral result{};
    for (unsigned s = 0; s < ruleCount; ++s) {
        const auto &point = rule[s];
        const auto bary = add(
            add(mul(patch.bary[0], point.barycentric[0]), mul(patch.bary[1], point.barycentric[1])),
            mul(patch.bary[2], point.barycentric[2]));
        const double basis[3]{bary.x, bary.y, bary.z};
        const double radius =
            length(sub(sub(observer, panel.vertices[0]), localPoint(panel, basis)));
        if (!(radius > 0)) {
            atomicExch(error, 1);
            return {};
        }
        double sine, cosine;
        sincos(k * radius, &sine, &cosine);
        const double scale =
            panel.area * patch.fraction * point.weight / (6.28318530717958647693 * radius);
        result.value = plus(result.value, {scale * cosine, -scale * sine});
        result.absoluteScale += scale;
    }
    return result;
}
__device__ DoubleComplex adaptiveLayer(Vec3 observer, const Panel &panel, double k,
                                       const galerkin::TriangleQuadraturePoint *rule,
                                       unsigned ruleCount, double tolerance, unsigned maximumDepth,
                                       int *error) {
    // Depth-first traversal needs at most 1+3*depth live entries.
    FieldPatch stack[52];
    stack[0] = {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}, 1, 0};
    unsigned pending = 1;
    DoubleComplex value{};
    while (pending) {
        const auto patch = stack[--pending];
        const auto coarse = layerPatch(observer, panel, patch, k, rule, ruleCount, error);
        FieldIntegral fine{};
        for (unsigned c = 0; c < 4; ++c) {
            const auto child =
                layerPatch(observer, panel, childPatch(patch, c), k, rule, ruleCount, error);
            fine.value = plus(fine.value, child.value);
            fine.absoluteScale += child.absoluteScale;
        }
        const double estimate =
            magnitude({fine.value.real - coarse.value.real, fine.value.imag - coarse.value.imag});
        if (estimate <= tolerance * fmax(fine.absoluteScale, 1e-300))
            value = plus(value, fine.value);
        else {
            if (patch.depth >= maximumDepth || pending + 4 > 52) {
                atomicExch(error, 3);
                return {};
            }
            for (unsigned c = 0; c < 4; ++c)
                stack[pending++] = childPatch(patch, c);
        }
    }
    return value;
}
__global__ void rayleighField(const Panel *panels, unsigned panelCount, unsigned pressureCount,
                              unsigned dimension, const int *apertures, const Complex *solution,
                              unsigned rhsCount, const Vec3 *observers, unsigned observerCount,
                              bool farField, Vec3 farOrigin, Vec3 baffleNormal, double distance,
                              double k, const galerkin::TriangleQuadraturePoint *rule,
                              unsigned ruleCount, double fieldTolerance, unsigned fieldDepth,
                              Complex *output, int *error) {
    const unsigned observer = blockIdx.x * blockDim.x + threadIdx.x, r = blockIdx.y;
    if (observer >= observerCount || r >= rhsCount)
        return;
    Vec3 position = observers[observer];
    if (farField) {
        position = sub(position, farOrigin);
        position = mul(position, 1 / length(position));
        if (dot(position, baffleNormal) < 0) {
            output[observer + r * observerCount] = {};
            return;
        }
    }
    DoubleComplex value{};
    for (unsigned t = 0; t < panelCount; ++t) {
        if (apertures[t] < 0)
            continue;
        const auto &panel = panels[t];
        const auto q = promoted(
            solution[pressureCount + apertures[t] + static_cast<std::size_t>(r) * dimension]);
        if (!farField) {
            const auto integral = adaptiveLayer(position, panel, k, rule, ruleCount, fieldTolerance,
                                                fieldDepth, error);
            value = plus(value, times(product(integral, q), -2 * k));
            continue;
        }
        for (unsigned s = 0; s < ruleCount; ++s) {
            const auto &point = rule[s];
            const auto y = add(panel.vertices[0], localPoint(panel, point.barycentric));
            double radius = 1, phase;
            if (farField)
                phase = -k * dot(position, sub(y, farOrigin));
            else {
                radius = length(sub(position, y));
                phase = k * radius;
            }
            if (!(radius > 0)) {
                atomicExch(error, 1);
                return;
            }
            double sine, cosine;
            sincos(phase, &sine, &cosine);
            const double scale = -k * panel.area * point.weight / (std::numbers::pi * radius);
            value = plus(value, product({scale * cosine, -scale * sine}, q));
        }
    }
    if (farField) {
        double sine, cosine;
        sincos(k * distance, &sine, &cosine);
        value = product(value, {cosine / distance, -sine / distance});
    }
    if (!isfinite(value.real) || !isfinite(value.imag))
        atomicExch(error, 1);
    output[observer + static_cast<std::size_t>(r) * observerCount] = {value.real, value.imag};
}
std::complex<double> host(Complex z) {
    return {z.x, z.y};
}
std::size_t productSize(std::size_t a, std::size_t b) {
    if (a && b > std::numeric_limits<std::size_t>::max() / a)
        throw std::invalid_argument("Baffled matrix size overflow.");
    return a * b;
}
} // namespace

BemResult solveBaffledCavity(const SurfaceMesh &mesh, Vec3 worldOrigin, const BemRequest &request,
                             BaffledPlanCache &planCache, std::stop_token stop,
                             std::function<void(const BemProgress &)> progress) {
    cancelled(stop);
    if (!request.infiniteBaffle)
        throw std::logic_error("Baffled request has no plane.");
    if (request.coupling && request.excitations.size() > 1)
        throw std::invalid_argument(
            "Baffled interface coupling permits one combined prescribed excitation.");
    if (request.numerics.acceleration != BemAcceleration::dense)
        throw std::invalid_argument(
            "Unendliche Schallwand: dichten GPU-Pfad wählen; FMM wird noch nicht unterstützt.");
    const auto &plane = *request.infiniteBaffle;
    if (!finite(plane.origin) || !finite(plane.normal) ||
        std::abs(length(plane.normal) - 1) > 1e-12 || plane.apertureFaceIds.empty())
        throw std::invalid_argument(
            "Unendliche Schallwand: Ursprung, Einheitsnormale und Aperturflächen prüfen.");
    if (!std::isfinite(request.medium.density) || request.medium.density <= 0 ||
        !std::isfinite(request.medium.soundSpeed) || request.medium.soundSpeed <= 0 ||
        request.frequencies.empty() || (!request.coupling && request.excitations.empty()))
        throw std::invalid_argument(
            "Baffled BEM requires medium, positive frequencies and excitations.");
    const double maximumK = bem_detail::planningWavenumber(
        request.frequencies, request.interactionPlanFrequencies, request.medium.soundSpeed);
    const auto quality = inspectMesh(mesh);
    if (!quality.watertight() || quality.connectedComponents != 1 || quality.signedVolume <= 0)
        throw std::invalid_argument(
            "Baffled BEM requires one closed, outward-oriented air cavity.");
    const auto space = galerkin::buildSpace(mesh);
    auto exteriorMesh = plane.exteriorDrivenMesh;
    for (auto &p : exteriorMesh.vertices)
        p = sub(p, worldOrigin);
    const bool hasExterior = !exteriorMesh.triangles.empty();
    if (hasExterior && request.axisymmetric)
        throw std::invalid_argument("Separate exterior sources currently require full 3D BEM.");
    bem_detail::validateBaffleExteriorMesh(mesh, exteriorMesh,
                                         sub(plane.origin, worldOrigin), plane.normal);
    const auto exteriorSpace = hasExterior ? galerkin::buildSpace(exteriorMesh) : galerkin::Space{};
    const auto exteriorCount = exteriorSpace.panels.size();
    std::vector<unsigned> pressureDofs(space.vertices.size()), fluxDofs(space.panels.size());
    std::iota(pressureDofs.begin(), pressureDofs.end(), 0u);
    std::iota(fluxDofs.begin(), fluxDofs.end(), 0u);
    std::size_t pressureCount = space.vertices.size(), fluxCount = space.panels.size();
    std::optional<galerkin::AxisymmetricSpace> ringSpace;
    if (request.axisymmetric) {
        const auto &reduction = *request.axisymmetric;
        ringSpace = galerkin::buildAxisymmetricSpace(space, worldOrigin, reduction);
        if (dot(reduction.axis, plane.normal) < 1 - 1e-12)
            throw std::invalid_argument(
                "Axisymmetric baffle normal must follow the rotation axis.");
        pressureDofs = reduction.pressureDofs;
        fluxDofs = reduction.fluxDofs;
        pressureCount = ringSpace->meridian.points.size();
        fluxCount = ringSpace->elements.size();
    }
    const auto measures = surfaceMeasures(mesh, request.axisymmetric, worldOrigin);
    const auto planeOrigin = sub(plane.origin, worldOrigin);
    double extent{};
    for (auto p : mesh.vertices)
        extent = std::max(extent, length(sub(p, planeOrigin)));
    const double planeTolerance = 128 * std::numeric_limits<double>::epsilon() * extent;
    std::set<std::uint32_t> apertureFaces(plane.apertureFaceIds.begin(),
                                          plane.apertureFaceIds.end());
    if (apertureFaces.size() != plane.apertureFaceIds.size())
        throw std::invalid_argument("Baffled aperture face IDs must be unique.");
    for (auto p : mesh.vertices)
        if (dot(sub(p, planeOrigin), plane.normal) > planeTolerance)
            throw std::invalid_argument("Unendliche Schallwand: Dieser Pfad verlangt eine "
                                        "vollständig versenkte Luftkavität.");
    std::vector<int> apertureIndex(space.panels.size(), -1);
    std::vector<int> fluxApertures(fluxCount, -1), stripApertureState(fluxCount, -1);
    std::vector<Panel> panels;
    std::vector<double> scales(pressureCount), pressureMass(pressureCount), apertureAreas;
    if (ringSpace)
        pressureMass = ringSpace->lumpedMass;
    else
        pressureMass = space.lumpedMass;
    for (std::size_t dof = 0; dof < pressureCount; ++dof)
        scales[dof] = 1 / pressureMass[dof];
    std::set<std::uint32_t> foundFaces;
    std::map<std::uint32_t, BoundaryCondition> faceConditions;
    for (std::size_t t = 0; t < space.panels.size(); ++t) {
        const auto &source = space.panels[t];
        const auto [it, inserted] = faceConditions.emplace(source.faceId, source.boundary);
        if (!inserted && it->second != source.boundary)
            throw std::invalid_argument("Baffled face group has inconsistent boundary conditions.");
        Panel target{};
        for (unsigned i = 0; i < 3; ++i) {
            target.vertices[i] = source.vertices[i];
            target.nodes[i] = source.nodes[i];
            target.curls[i] = source.curls[i];
        }
        target.normal = source.normal;
        target.centroid = source.centroid;
        target.area = source.area;
        target.edge = source.maximumEdgeLength;
        if (!ringSpace)
            panels.push_back(target);
        const int apertureState = apertureFaces.contains(source.faceId) ? 1 : 0;
        auto &stripState = stripApertureState[fluxDofs[t]];
        if (stripState >= 0 && stripState != apertureState)
            throw std::invalid_argument(
                "Axisymmetric flux strip mixes aperture and physical surfaces.");
        stripState = apertureState;
        if (apertureFaces.contains(source.faceId)) {
            foundFaces.insert(source.faceId);
            if (source.boundary != BoundaryCondition::wall ||
                dot(source.normal, plane.normal) < 1 - 1e-12)
                throw std::invalid_argument("Aperturflächen müssen Wall-Tags und nach außen "
                                            "gerichtete Ebenennormalen besitzen.");
            for (auto p : source.vertices)
                if (std::abs(dot(sub(p, planeOrigin), plane.normal)) > planeTolerance)
                    throw std::invalid_argument("Baffled aperture is not on the specified plane.");
            auto &index = fluxApertures[fluxDofs[t]];
            if (index < 0) {
                index = static_cast<int>(apertureAreas.size());
                apertureAreas.push_back(0);
            }
            apertureIndex[t] = index;
            apertureAreas[index] += measures.areas[t];
        }
    }
    for (double area : apertureAreas)
        scales.push_back(1 / area);
    if (foundFaces != apertureFaces)
        throw std::invalid_argument("Baffled aperture references an unknown face.");
    std::vector<unsigned> aperturePanels(apertureAreas.size());
    if (hasExterior) {
        for (std::size_t t = 0; t < space.panels.size(); ++t)
            if (apertureIndex[t] >= 0)
                aperturePanels[apertureIndex[t]] = static_cast<unsigned>(t);
        for (const auto &panel : exteriorSpace.panels)
            faceConditions.emplace(panel.faceId, panel.boundary);
    }
    std::optional<bem_detail::CouplingBasis> couplingBasis;
    std::vector<double> panelAreas;
    auto excitations = request.excitations;
    const auto knownCount = excitations.size();
    std::set<std::uint32_t> portFaces;
    if (request.coupling) {
        for (const auto *element : lem::bemPorts(*request.coupling)) {
            const auto &port = std::get<lem::BemPort>(element->component);
            for (auto id : port.faceIds) {
                if (apertureFaces.contains(id))
                    throw std::invalid_argument(
                        "A LEM port must be disjoint from the Rayleigh aperture.");
                portFaces.insert(id);
            }
        }
        panelAreas = measures.areas;
        for (const auto &panel : exteriorSpace.panels)
            panelAreas.push_back(panel.area);
        const auto couplingMesh = bem_detail::appendBaffleExterior(mesh, exteriorMesh);
        couplingBasis = bem_detail::makeCouplingBasis(*request.coupling, couplingMesh, panelAreas,
                                                      bem_detail::CouplingDomain::baffledAir);
        excitations.insert(excitations.end(), couplingBasis->excitations.begin(),
                           couplingBasis->excitations.end());
    }
    const std::size_t count = scales.size(), panelCount = space.panels.size(),
                      rhsCount = excitations.size();
    if (count > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
        panelCount > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
        rhsCount > 65535 ||
        request.observers.size() > static_cast<std::size_t>(std::numeric_limits<unsigned>::max()))
        throw std::invalid_argument("Baffled dense solver index capacity exceeded.");
    std::vector<std::map<std::uint32_t, std::complex<double>>> velocities;
    for (std::size_t r = 0; r < excitations.size(); ++r) {
        const auto &excitation = excitations[r];
        const bool portBasis = couplingBasis && r >= knownCount;
        std::map<std::uint32_t, std::complex<double>> values;
        bool nonzero = false;
        for (const auto &face : excitation.faces) {
            if (!faceConditions.contains(face.faceId) ||
                (!portBasis && faceConditions.at(face.faceId) != BoundaryCondition::driven) ||
                portFaces.contains(face.faceId) != portBasis ||
                apertureFaces.contains(face.faceId) || !std::isfinite(face.normalVelocity.real()) ||
                !std::isfinite(face.normalVelocity.imag()) ||
                !values.emplace(face.faceId, face.normalVelocity).second)
                throw std::invalid_argument(
                    "Baffled excitation must reference unique physical Driven faces, "
                    "disjoint from LEM port and aperture faces.");
            nonzero = nonzero || std::abs(face.normalVelocity) > 0;
        }
        if (!nonzero && (!couplingBasis || excitation.faces.empty()))
            throw std::invalid_argument("Baffled excitation has no nonzero driven velocity.");
        velocities.push_back(std::move(values));
    }
    if (request.axisymmetric)
        for (const auto &values : velocities) {
            std::vector<std::complex<double>> shape(fluxCount);
            std::vector<bool> seen(fluxCount);
            for (std::size_t t = 0; t < panelCount; ++t) {
                const auto found = values.find(space.panels[t].faceId);
                const auto velocity =
                    found == values.end() ? std::complex<double>{} : found->second;
                const auto strip = fluxDofs[t];
                if (seen[strip] && shape[strip] != velocity)
                    throw std::invalid_argument(
                        "Axisymmetric excitation or LEM mode varies around a meridian strip.");
                seen[strip] = true;
                shape[strip] = velocity;
            }
        }
    std::vector<Vec3> observers;
    const Vec3 farOrigin = request.farField ? sub(request.farField->origin, worldOrigin) : Vec3{};
    if (request.farField &&
        (!finite(farOrigin) || !std::isfinite(request.farField->referenceDistance) ||
         request.farField->referenceDistance <= 0))
        throw std::invalid_argument("Baffled far-field origin and distance must be finite.");
    for (auto observer : request.observers) {
        const auto p = sub(observer, worldOrigin);
        if (!finite(p) || (request.farField && length(sub(p, farOrigin)) == 0))
            throw std::invalid_argument("Baffled observer is invalid.");
        if (!request.farField && dot(sub(p, planeOrigin), plane.normal) <= planeTolerance)
            throw std::invalid_argument(
                "Endliche Schallwand-Beobachter müssen oberhalb der Ebene liegen.");
        observers.push_back(p);
    }
    const auto &n = request.numerics;
    if (!std::isfinite(n.relativeResidualTolerance) || n.relativeResidualTolerance <= 0 ||
        n.maximumRefinementSteps > 20)
        throw std::invalid_argument("Invalid baffled residual settings.");
    const auto matrixElements = productSize(count, count),
               rhsElements = productSize(count, rhsCount);
    const auto minimumBytes = productSize(matrixElements, 2 * sizeof(Complex)) +
                              productSize(rhsElements, 4 * sizeof(Complex));
    std::size_t free{}, total{};
    check(cudaMemGetInfo(&free, &total), "Baffled memory query");
    if (minimumBytes > free - free / 10)
        throw std::runtime_error(
            "Unendliche Schallwand: dichte GPU-Matrix überschreitet den verfügbaren Speicher.");
    galerkin::PairQuadratureSettings qs;
    qs.regularOrder = n.quadratureOrder;
    qs.singularOrder = n.singularOrder;
    qs.nearDistanceRatio = n.nearDistanceRatio;
    qs.relativeTolerance = n.nearRelativeTolerance;
    qs.maximumSubdivisionDepth = n.nearMaximumDepth;
    std::shared_ptr<const galerkin::InteractionPlan> cartesianPlan;
    std::shared_ptr<const galerkin::RingInteractionPlan> ringPlan;
    double planMilliseconds{};
    bool planReused{};
    if (ringSpace) {
        const auto planned =
            planCache.acquireRing(mesh, *ringSpace, n, *request.axisymmetric, stop);
        ringPlan = planned.plan;
        planMilliseconds = planned.buildMilliseconds;
        planReused = planned.reused;
    } else {
        const auto planned = planCache.acquire(mesh, space, qs, maximumK, {}, {}, stop);
        cartesianPlan = planned.plan;
        planMilliseconds = planned.buildMilliseconds;
        planReused = planned.reused;
    }
    const galerkin::InteractionPlan emptyPlan;
    const auto &plan = cartesianPlan ? *cartesianPlan : emptyPlan;
    const auto regular = galerkin::triangleRule(n.quadratureOrder);
    const auto gauss = galerkin::gaussLegendre01(n.singularOrder);
    auto nearRules = galerkin::triangleRule(8);
    const auto high = galerkin::triangleRule(14), fieldRule = galerkin::triangleRule(14);
    const auto baseCount = nearRules.size();
    nearRules.insert(nearRules.end(), high.begin(), high.end());
    Device<Panel> dp(panels.size());
    dp.upload(panels);
    Device<int> da(apertureIndex.size());
    da.upload(apertureIndex);
    Device<double> ds(scales.size());
    ds.upload(scales);
    Device<galerkin::TriangleQuadraturePoint> dr(regular.size());
    dr.upload(regular);
    Device<galerkin::TriangleQuadraturePoint> dnr(nearRules.size());
    dnr.upload(nearRules);
    Device<galerkin::TriangleQuadraturePoint> dfr(fieldRule.size());
    dfr.upload(fieldRule);
    Device<galerkin::GaussPoint> dg(gauss.size());
    dg.upload(gauss);
    Device<galerkin::NearInteraction> dn(plan.near.size());
    dn.upload(plan.near);
    Device<unsigned> dno(plan.nearOffsets.size());
    dno.upload(plan.nearOffsets);
    Device<unsigned> dnt(plan.nearTrials.size());
    dnt.upload(plan.nearTrials);
    Device<galerkin::PairQuadraturePatch> dpatch(plan.patches.size());
    dpatch.upload(plan.patches);
    Device<Vec3> dob(observers.size());
    dob.upload(observers);
    std::unique_ptr<AxisymmetricAssembly> ringAssembly;
    std::vector<RingObservation> ringObservers;
    if (ringSpace) {
        std::vector<unsigned> representatives(fluxCount, std::numeric_limits<unsigned>::max()),
            faces;
        for (std::size_t t = 0; t < panelCount; ++t)
            if (representatives[fluxDofs[t]] == std::numeric_limits<unsigned>::max())
                representatives[fluxDofs[t]] = static_cast<unsigned>(t);
        for (const auto &e : ringSpace->elements)
            faces.push_back(e.faceId);
        ringAssembly = std::make_unique<AxisymmetricAssembly>(*ringPlan, fluxApertures,
                                                              representatives, faces);
        const auto axisOrigin = sub(ringSpace->meridian.origin, worldOrigin);
        const auto axis = ringSpace->meridian.axis;
        for (const auto &p : observers) {
            RingObservation value;
            if (request.farField) {
                const auto direction = mul(sub(p, farOrigin), 1 / length(sub(p, farOrigin)));
                value.directionZ = dot(direction, axis);
                value.directionR = length(sub(direction, mul(axis, value.directionZ)));
                value.originPhase = dot(direction, sub(axisOrigin, farOrigin));
                value.blocked = dot(direction, plane.normal) < 0;
            } else {
                const auto relative = sub(p, axisOrigin);
                value.axial = dot(relative, axis);
                value.radius = length(sub(relative, mul(axis, value.axial)));
            }
            ringObservers.push_back(value);
        }
    }
    Device<Complex> matrix(matrixElements), factor(matrixElements), rhs(rhsElements),
        solution(rhsElements), residual(rhsElements), correction(rhsElements),
        flux(productSize(panelCount, rhsCount)), fields(productSize(observers.size(), rhsCount));
    Device<int> pivots(count), info(1), error(1);
    Handles handles;
    int workspaceCount{};
    check(cusolverDnZgetrf_bufferSize(handles.solver, static_cast<int>(count),
                                      static_cast<int>(count), factor.data(),
                                      static_cast<int>(count), &workspaceCount),
          "Baffled LU workspace query");
    if (workspaceCount < 0)
        throw std::runtime_error("Invalid baffled LU workspace size.");
    Device<Complex> workspace(static_cast<std::size_t>(workspaceCount));
    BemResult output;
    output.device = CudaBemSolver::deviceDescription();
    output.triangleCount = panelCount;
    output.pressureDofCount = pressureCount;
    output.acceleration = BemAcceleration::dense;
    output.geometryMode =
        request.axisymmetric ? BemGeometryMode::axisymmetric : BemGeometryMode::threeDimensional;
    output.peakDeviceBytes =
        dp.bytes() + da.bytes() + ds.bytes() + dr.bytes() + dnr.bytes() + dfr.bytes() + dg.bytes() +
        dn.bytes() + dno.bytes() + dnt.bytes() + dpatch.bytes() + dob.bytes() + matrix.bytes() +
        factor.bytes() + rhs.bytes() + solution.bytes() + residual.bytes() + correction.bytes() +
        flux.bytes() + fields.bytes() + pivots.bytes() + info.bytes() + error.bytes() +
        workspace.bytes() + (ringAssembly ? ringAssembly->bytes() : 0);
    const unsigned pc = static_cast<unsigned>(pressureCount),
                   tc = static_cast<unsigned>(panelCount), dim = static_cast<unsigned>(count),
                   nrhs = static_cast<unsigned>(rhsCount);
    for (std::size_t fi = 0; fi < request.frequencies.size(); ++fi) {
        cancelled(stop);
        BemFrequencyResult frequency;
        frequency.interactionPlanMilliseconds = fi == 0 ? planMilliseconds : 0;
        frequency.interactionPlanReused = fi == 0 && planReused;
        frequency.frequency = request.frequencies[fi];
        const double omega = 2 * std::numbers::pi * frequency.frequency,
                     k = omega / request.medium.soundSpeed;
        std::vector<Complex> knownFlux(panelCount * rhsCount);
        for (std::size_t r = 0; r < rhsCount; ++r)
            for (std::size_t t = 0; t < panelCount; ++t) {
                const auto it = velocities[r].find(space.panels[t].faceId);
                if (it != velocities[r].end()) {
                    // Public positive velocity points from the solid INTO AIR.
                    const auto q =
                        std::complex<double>{0, request.medium.density * omega} * it->second;
                    knownFlux[t + r * panelCount] = {q.real(), q.imag()};
                }
            }
        const auto assemblyStart = Clock::now();
        const auto exteriorOperator = hasExterior
            ? bem_detail::baffleExteriorOperator(space, aperturePanels, exteriorSpace, n, k, stop)
            : bem_detail::BaffleExteriorOperator{};
        const auto exteriorPressure = [&](const std::vector<Complex> &values, std::size_t r,
                                         const auto &sourceVelocities) {
            std::vector<std::complex<double>> pressure(exteriorCount);
            for (std::size_t e = 0; e < exteriorCount; ++e) {
                for (std::size_t a = 0; a < aperturePanels.size(); ++a)
                    pressure[e] -= 2 * k * exteriorOperator.cross[e * aperturePanels.size() + a] *
                                   host(values[pressureCount + a + r * count]);
                for (std::size_t j = 0; j < exteriorCount; ++j) {
                    const auto found = sourceVelocities.find(exteriorSpace.panels[j].faceId);
                    if (found != sourceVelocities.end())
                        pressure[e] += std::complex<double>{0, 2 * request.medium.density * omega} *
                                       exteriorOperator.self[e * exteriorCount + j] * found->second;
                }
            }
            return pressure;
        };
        matrix.clear();
        rhs.clear();
        error.clear();
        flux.upload(knownFlux);
        if (ringAssembly)
            ringAssembly->assemble(pc, dim, tc, nrhs, k, n.nearRelativeTolerance,
                                   n.nearMaximumDepth, ds.data(), flux.data(), matrix.data(),
                                   rhs.data(), stop);
        else {
            for (unsigned first = 0; first < tc; first += 16) {
                cancelled(stop);
                assemblePairs<<<dim3(tc, std::min(16u, tc - first)), 64>>>(
                    dp.data(), tc, pc, dim, da.data(), ds.data(), flux.data(), nrhs, k, dr.data(),
                    static_cast<unsigned>(regular.size()), dg.data(), n.singularOrder, dn.data(),
                    dno.data(), dnt.data(), dpatch.data(), dnr.data(),
                    static_cast<unsigned>(baseCount), static_cast<unsigned>(high.size()),
                    matrix.data(), rhs.data(), error.data(), first);
                check(cudaDeviceSynchronize(), "Baffled pair assembly");
            }
            addInteriorMass<<<(tc + 127) / 128, 128>>>(dp.data(), tc, pc, dim, da.data(), ds.data(),
                                                       flux.data(), nrhs, k, matrix.data(),
                                                       rhs.data(), error.data());
            check(cudaDeviceSynchronize(), "Baffled mass assembly");
        }
        if (error.download()[0])
            throw std::runtime_error("Baffled Galerkin integration is not finite.");
        if (hasExterior) {
            auto exteriorRhs = rhs.download();
            for (std::size_t r = 0; r < rhsCount; ++r)
                for (std::size_t a = 0; a < aperturePanels.size(); ++a) {
                    std::complex<double> value{};
                    for (std::size_t e = 0; e < exteriorCount; ++e) {
                        const auto found = velocities[r].find(exteriorSpace.panels[e].faceId);
                        if (found != velocities[r].end())
                            value += exteriorSpace.panels[e].area *
                                exteriorOperator.cross[e * aperturePanels.size() + a] * found->second;
                    }
                    value *= std::complex<double>{0, 2 * request.medium.density * omega} *
                             scales[pressureCount + a];
                    auto &entry = exteriorRhs[pressureCount + a + r * count];
                    value += host(entry);
                    entry = {value.real(), value.imag()};
                }
            rhs.upload(exteriorRhs);
        }
        frequency.assemblyMilliseconds = elapsed(assemblyStart) + (fi == 0 ? planMilliseconds : 0);
        cancelled(stop);
        check(cudaMemcpy(factor.data(), matrix.data(), matrix.bytes(), cudaMemcpyDeviceToDevice),
              "Baffled matrix copy");
        check(cudaMemcpy(solution.data(), rhs.data(), rhs.bytes(), cudaMemcpyDeviceToDevice),
              "Baffled RHS copy");
        const auto factorStart = Clock::now();
        check(cusolverDnZgetrf(handles.solver, static_cast<int>(count), static_cast<int>(count),
                               factor.data(), static_cast<int>(count), workspace.data(),
                               pivots.data(), info.data()),
              "Baffled GPU LU");
        check(cudaDeviceSynchronize(), "Baffled LU synchronization");
        if (info.download()[0] != 0)
            throw std::runtime_error("Baffled GPU factorization is singular.");
        frequency.factorizationMilliseconds = elapsed(factorStart);
        const auto solveStart = Clock::now();
        auto triangular = [&](Complex *value, std::size_t activeCount) {
            check(cusolverDnZgetrs(handles.solver, CUBLAS_OP_N, static_cast<int>(count),
                                   static_cast<int>(activeCount), factor.data(),
                                   static_cast<int>(count), pivots.data(), value,
                                   static_cast<int>(count), info.data()),
                  "Baffled GPU triangular solve");
            check(cudaDeviceSynchronize(), "Baffled solve synchronization");
            if (info.download()[0])
                throw std::runtime_error("Baffled GPU triangular solve failed.");
        };
        triangular(solution.data(), rhsCount);
        frequency.triangularSolveMilliseconds = elapsed(solveStart);
        auto rhsHost = rhs.download();
        std::vector<double> norms(rhsCount), relative(rhsCount);
        unsigned refinements{};
        auto refine = [&](std::size_t activeCount) {
            const auto residualStart = Clock::now();
            std::fill(norms.begin(), norms.end(), 0);
            for (std::size_t r = 0; r < activeCount; ++r)
                for (std::size_t i = 0; i < count; ++i)
                    norms[r] =
                        std::hypot(norms[r], std::abs(host(rhsHost[i + r * count])) / scales[i]);
            unsigned steps{};
            for (;;) {
                cancelled(stop);
                check(
                    cudaMemcpy(residual.data(), rhs.data(), rhs.bytes(), cudaMemcpyDeviceToDevice),
                    "Baffled residual reset");
                const Complex one{1, 0}, minus{-1, 0};
                check(cublasZgemm(handles.blas, CUBLAS_OP_N, CUBLAS_OP_N, static_cast<int>(count),
                                  static_cast<int>(activeCount), static_cast<int>(count), &one,
                                  matrix.data(), static_cast<int>(count), solution.data(),
                                  static_cast<int>(count), &minus, residual.data(),
                                  static_cast<int>(count)),
                      "Baffled original residual");
                const auto values = residual.download();
                std::fill(relative.begin(), relative.end(), 0);
                for (std::size_t r = 0; r < activeCount; ++r) {
                    for (std::size_t i = 0; i < count; ++i)
                        relative[r] = std::hypot(relative[r],
                                                 std::abs(host(values[i + r * count])) / scales[i]);
                    relative[r] /= std::max(norms[r], std::numeric_limits<double>::min());
                    if (!std::isfinite(relative[r]))
                        throw std::runtime_error("Baffled residual is not finite.");
                }
                if (*std::max_element(relative.begin(), relative.end()) <=
                    n.relativeResidualTolerance)
                    break;
                if (steps >= n.maximumRefinementSteps)
                    throw std::runtime_error(
                        "Baffled original block residual exceeds the requested tolerance.");
                check(cudaMemcpy(correction.data(), residual.data(), residual.bytes(),
                                 cudaMemcpyDeviceToDevice),
                      "Baffled correction reset");
                triangular(correction.data(), activeCount);
                check(cublasZaxpy(handles.blas, static_cast<int>(count * activeCount), &minus,
                                  correction.data(), 1, solution.data(), 1),
                      "Baffled refinement update");
                ++refinements;
                ++steps;
            }
            frequency.residualCheckMilliseconds += elapsed(residualStart);
        };
        refine(rhsCount);

        auto actualExcitations = request.excitations;
        std::optional<lem::Solution> lemSolution;
        if (couplingBasis) {
            cancelled(stop);
            const auto basisPressure = solution.download();
            const auto portCount = couplingBasis->areas.size();
            lem::AcousticLoad load;
            load.impedance.resize(productSize(portCount, portCount));
            load.blockedPressure.resize(portCount);
            // The unit-Q velocity shape is power-dual to modal p. Keep every mutual term
            // and the independently prescribed, frequency-dependent blocked p0.
            for (std::size_t t = 0; t < panelCount; ++t) {
                const auto port = couplingBasis->panelPorts[t];
                if (port < 0)
                    continue;
                const double weight = panelAreas[t] * couplingBasis->panelVelocities[t];
                for (std::size_t r = 0; r < rhsCount; ++r) {
                    std::complex<double> mean{};
                    for (unsigned local = 0; local < 3; ++local) {
                        const auto node = space.panels[t].nodes[local];
                        const auto value = host(basisPressure[pressureDofs[node] + r * count]);
                        if (ringSpace)
                            mean +=
                                value * (measures.pressureIntegrals[t][local] / measures.areas[t]);
                        else
                            mean += value / 3.0;
                    }
                    if (r < knownCount)
                        load.blockedPressure[port] += weight * mean;
                    else
                        load.impedance[static_cast<std::size_t>(port) * portCount + r -
                                       knownCount] += weight * mean;
                }
            }
            for (std::size_t r = 0; r < rhsCount; ++r) {
                const auto pressure = exteriorPressure(basisPressure, r, velocities[r]);
                for (std::size_t e = 0; e < exteriorCount; ++e) {
                    const auto t = panelCount + e;
                    const auto port = couplingBasis->panelPorts[t];
                    if (port < 0)
                        continue;
                    const auto weighted = panelAreas[t] * couplingBasis->panelVelocities[t] * pressure[e];
                    if (r < knownCount)
                        load.blockedPressure[port] += weighted;
                    else
                        load.impedance[static_cast<std::size_t>(port) * portCount + r - knownCount] += weighted;
                }
            }
            lemSolution = lem::solve(*request.coupling, frequency.frequency, load,
                                     {request.medium.density, request.medium.soundSpeed}, stop);
            auto actualVelocities = bem_detail::coupledVelocities(*couplingBasis, *lemSolution);
            if (knownCount)
                actualVelocities.insert(actualVelocities.end(),
                                        request.excitations[0].faces.begin(),
                                        request.excitations[0].faces.end());
            actualExcitations = {{"LEM", std::move(actualVelocities)}};

            // Reconstruct pressure AND aperture flux in the original closed-air
            // block space (full 3D or explicit m=0). Audit/refine A*x=b directly
            // with the retained matrix; basis residuals do not estimate the
            // coupled residual or the discretization error against full 3D.
            std::vector<Complex> actualSolution(rhsElements), actualRhs(rhsElements);
            for (std::size_t i = 0; i < count; ++i) {
                std::complex<double> value{}, right{};
                for (std::size_t r = 0; r < rhsCount; ++r) {
                    const auto amplitude = r < knownCount
                                               ? std::complex<double>{1, 0}
                                               : lemSolution->ports[r - knownCount].volumeVelocity;
                    value += host(basisPressure[i + r * count]) * amplitude;
                    right += host(rhsHost[i + r * count]) * amplitude;
                }
                if (!std::isfinite(value.real()) || !std::isfinite(value.imag()) ||
                    !std::isfinite(right.real()) || !std::isfinite(right.imag()))
                    throw std::runtime_error("Coupled baffled solution is not finite.");
                actualSolution[i] = {value.real(), value.imag()};
                actualRhs[i] = {right.real(), right.imag()};
            }
            solution.upload(actualSolution);
            rhs.upload(actualRhs);
            rhsHost = std::move(actualRhs);
            refine(1);
        }
        frequency.solveMilliseconds = frequency.factorizationMilliseconds + elapsed(solveStart);
        const auto fieldStart = Clock::now();
        const auto outputCount = actualExcitations.size();
        const auto fieldCount = static_cast<unsigned>(outputCount);
        if (!observers.empty() && (!hasExterior || request.farField)) {
            if (ringAssembly)
                ringAssembly->field(ringObservers, request.farField.has_value(),
                                    request.farField ? request.farField->referenceDistance : 1, pc,
                                    dim, fieldCount, k, n.nearRelativeTolerance, n.nearMaximumDepth,
                                    solution.data(), fields.data(), stop);
            else {
                rayleighField<<<
                    dim3((static_cast<unsigned>(observers.size()) + 63) / 64, fieldCount), 64>>>(
                    dp.data(), tc, pc, dim, da.data(), solution.data(), fieldCount, dob.data(),
                    static_cast<unsigned>(observers.size()), request.farField.has_value(),
                    farOrigin, plane.normal,
                    request.farField ? request.farField->referenceDistance : 1, k, dfr.data(),
                    static_cast<unsigned>(fieldRule.size()), n.nearRelativeTolerance,
                    n.nearMaximumDepth, fields.data(), error.data());
                check(cudaDeviceSynchronize(), "Baffled Rayleigh field");
            }
            if (error.download()[0])
                throw std::runtime_error("Baffled Rayleigh field is not finite.");
        }
        const auto x = solution.download(), field = fields.download();
        frequency.fieldMilliseconds = elapsed(fieldStart);
        for (std::size_t r = 0; r < outputCount; ++r) {
            BemExcitationResult result;
            BemBoundaryData boundary;
            boundary.surfacePressure.reserve(panelCount);
            boundary.normalDerivative.reserve(panelCount);
            boundary.pressureTrace.reserve(space.vertices.size());
            result.name = actualExcitations[r].name;
            result.boundaryResidual = relative[r];
            result.refinementSteps = refinements;
            result.solvedVelocities = actualExcitations[r].faces;
            result.lemSolution = lemSolution;
            std::map<std::uint32_t, std::complex<double>> actualVelocities;
            for (const auto &face : result.solvedVelocities)
                actualVelocities.emplace(face.faceId, face.normalVelocity);
            std::map<std::uint32_t, std::pair<double, std::complex<double>>> ports;
            for (std::size_t vertex = 0; vertex < space.vertices.size(); ++vertex)
                boundary.pressureTrace.push_back(host(x[pressureDofs[vertex] + r * count]));
            for (std::size_t t = 0; t < panelCount; ++t) {
                std::complex<double> mean{};
                for (unsigned local = 0; local < 3; ++local) {
                    const auto node = space.panels[t].nodes[local];
                    if (ringSpace)
                        mean += boundary.pressureTrace[node] *
                                (measures.pressureIntegrals[t][local] / measures.areas[t]);
                    else
                        mean += boundary.pressureTrace[node] / 3.0;
                }
                boundary.surfacePressure.push_back(mean);
                if (apertureIndex[t] >= 0)
                    boundary.normalDerivative.push_back(k *
                                                      host(x[pc + apertureIndex[t] + r * count]));
                else {
                    const auto found = actualVelocities.find(space.panels[t].faceId);
                    boundary.normalDerivative.push_back(
                        found == actualVelocities.end()
                            ? std::complex<double>{}
                            : std::complex<double>{0, request.medium.density * omega} *
                                  found->second);
                }
                if (actualVelocities.contains(space.panels[t].faceId)) {
                    auto &port = ports[space.panels[t].faceId];
                    port.first += measures.areas[t];
                    port.second += measures.areas[t] * mean;
                }
            }
            boundary.exteriorSurfacePressure = exteriorPressure(x, r, actualVelocities);
            for (std::size_t e = 0; e < exteriorCount; ++e) {
                const auto &panel = exteriorSpace.panels[e];
                if (actualVelocities.contains(panel.faceId)) {
                    auto &port = ports[panel.faceId];
                    port.first += panel.area;
                    port.second += panel.area * boundary.exteriorSurfacePressure[e];
                }
            }
            if (couplingBasis) {
                auto allPressure = boundary.surfacePressure;
                allPressure.insert(allPressure.end(), boundary.exteriorSurfacePressure.begin(),
                                   boundary.exteriorSurfacePressure.end());
                bem_detail::validateCouplingClosure(*couplingBasis, *result.lemSolution, panelAreas,
                                                    allPressure);
            }
            result.setBoundary(std::move(boundary));
            if (hasExterior && !request.farField && !observers.empty()) {
                auto worldMesh = mesh;
                for (auto &p : worldMesh.vertices)
                    p = add(p, worldOrigin);
                const auto samples = evaluatePressureField(worldMesh, request, frequency.frequency,
                                                           result, request.observers, {}, stop);
                for (const auto &sample : samples) {
                    if (sample.status != FieldPointStatus::air)
                        throw std::runtime_error("Coupled exterior field quadrature failed.");
                    result.observerPressure.push_back(sample.pressure);
                }
            } else for (std::size_t o = 0; o < observers.size(); ++o) {
                auto value = host(field[o + r * observers.size()]);
                if (hasExterior) {
                    const auto direction = mul(sub(observers[o], farOrigin),
                                               1 / length(sub(observers[o], farOrigin)));
                    if (dot(direction, plane.normal) >= 0) {
                        std::complex<double> direct{};
                        for (const auto &panel : exteriorSpace.panels) {
                            const auto found = actualVelocities.find(panel.faceId);
                            if (found == actualVelocities.end())
                                continue;
                            for (const auto &point : fieldRule) {
                                Vec3 y{};
                                for (unsigned i = 0; i < 3; ++i)
                                    y = add(y, mul(panel.vertices[i], point.barycentric[i]));
                                direct += panel.area * point.weight / std::numbers::pi *
                                          std::exp(std::complex<double>{0, k * dot(direction, sub(y, farOrigin))}) *
                                          std::complex<double>{0, request.medium.density * omega} * found->second;
                            }
                        }
                        value += direct * std::exp(std::complex<double>{0, -k * request.farField->referenceDistance}) /
                                 request.farField->referenceDistance;
                    }
                }
                result.observerPressure.push_back(value);
            }
            for (const auto &[face, port] : ports) {
                const auto velocity = actualVelocities.at(face);
                if (std::abs(velocity) == 0)
                    continue;
                const auto mean = port.second / port.first, volume = velocity * port.first;
                result.impedances.push_back({face, port.first, mean, volume, mean / volume});
            }
            frequency.excitations.push_back(std::move(result));
        }
        output.frequencies.push_back(std::move(frequency));
        if (progress)
            progress({fi + 1, request.frequencies.size(), request.frequencies[fi]});
    }
    return output;
}
} // namespace hornsculpt::cuda_detail
