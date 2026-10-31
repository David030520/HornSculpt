#include "subdomain_solver.hpp"
#include <chrono>
#include <cmath>
#include <stdexcept>
#include <string>

namespace hornsculpt::cuda_detail {
namespace {
__device__ void addValue(cuDoubleComplex *destination, DoubleComplex value, int *error) {
    if (!isfinite(value.real) || !isfinite(value.imag)) {
        atomicExch(error, 1);
        return;
    }
    atomicAdd(&destination->x, value.real);
    atomicAdd(&destination->y, value.imag);
}
// Interior trace identity in baffled_solver.cu, using the same outward-air
// orientation for every domain, including the unbounded exterior. Interface
// q is represented in P1, independently of physical panel DP0 fluxes.
__global__ void assembleSamples(const SubdomainPanel *panels, const SubdomainPairSample *samples,
                                std::size_t count, std::size_t panelCount, unsigned dimension,
                                double k, unsigned excitationCount,
                                const cuDoubleComplex *knownFlux, const double *rowScale,
                                cuDoubleComplex *matrix, cuDoubleComplex *rhs, int *error) {
    const auto index = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (index >= count)
        return;
    const auto sample = samples[index];
    const auto &test = panels[sample.testPanel];
    const auto &trial = panels[sample.trialPanel];
    const auto &x = test.panel;
    const auto &y = trial.panel;
    const auto delta = pairDelta(x, y, sample.point.test, sample.point.trial);
    const double r = length(delta);
    if (!(r > 0) || !isfinite(r)) {
        atomicExch(error, 1);
        return;
    }
    double sine, cosine;
    sincos(k * r, &sine, &cosine);
    const double factor = sample.point.weight * 4 * x.area * y.area / (12.56637061435917295385 * r);
    const DoubleComplex g{cosine * factor, -sine * factor};
    const DoubleComplex radial{g.real - g.imag * k * r, g.imag + g.real * k * r};
    const bool coincident = sample.testPanel == sample.trialPanel;
    const auto d = times(radial, coincident ? 0 : dot(delta, y.normal) / (r * r));
    const auto kp = times(radial, coincident ? 0 : -dot(delta, x.normal) / (r * r));
    const double nn = dot(x.normal, y.normal);
    const DoubleComplex a{d.real - k * nn * g.imag, d.imag + k * nn * g.real};
    const DoubleComplex b{g.real - kp.imag / k, g.imag + kp.real / k};
    for (unsigned i = 0; i < 3; ++i) {
        const auto row = test.rows[i];
        const double scale = rowScale[row];
        for (unsigned j = 0; j < 3; ++j) {
            const double basis = sample.point.test[i] * sample.point.trial[j];
            const double curl = dot(x.curls[i], y.curls[j]) / k;
            const auto av = plus(times(a, basis), {g.imag * curl, -g.real * curl});
            addValue(matrix + row + static_cast<std::size_t>(trial.pressureColumns[j]) * dimension,
                     times(av, scale), error);
            if (trial.fluxColumns[j] >= 0) {
                // Unknown is q/k, so the operator column contains -k B.
                addValue(matrix + row + static_cast<std::size_t>(trial.fluxColumns[j]) * dimension,
                         times(b, -k * trial.fluxSign * basis * scale), error);
            }
        }
        if (trial.fluxColumns[0] < 0) {
            for (unsigned e = 0; e < excitationCount; ++e) {
                const auto q = promoted(knownFlux[sample.trialPanel + e * panelCount]);
                addValue(rhs + row + static_cast<std::size_t>(e) * dimension,
                         times(product(b, q), sample.point.test[i] * scale), error);
            }
        }
    }
}
__global__ void assembleMass(const SubdomainPanel *panels, std::size_t count, unsigned dimension,
                             double k, unsigned excitationCount, const cuDoubleComplex *knownFlux,
                             const double *rowScale, cuDoubleComplex *matrix, cuDoubleComplex *rhs,
                             int *error) {
    const auto t = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (t >= count)
        return;
    const auto &p = panels[t];
    for (unsigned i = 0; i < 3; ++i) {
        const auto row = p.rows[i];
        for (unsigned j = 0; j < 3; ++j) {
            const double m = p.panel.area * (i == j ? 2.0 : 1.0) / 12.0;
            addValue(matrix + row + static_cast<std::size_t>(p.pressureColumns[j]) * dimension,
                     {rowScale[row] * m / 2, 0}, error);
            if (p.fluxColumns[j] >= 0)
                addValue(matrix + row + static_cast<std::size_t>(p.fluxColumns[j]) * dimension,
                         {0, rowScale[row] * m * p.fluxSign / 2}, error);
        }
        if (p.fluxColumns[0] < 0) {
            for (unsigned e = 0; e < excitationCount; ++e) {
                const auto q = promoted(knownFlux[t + e * count]);
                const DoubleComplex coefficient{0, -rowScale[row] * p.panel.area / (6 * k)};
                addValue(rhs + row + static_cast<std::size_t>(e) * dimension,
                         product(coefficient, q), error);
            }
        }
    }
}
void checkLaunch(const char *operation) {
    const auto status = cudaGetLastError();
    if (status != cudaSuccess)
        throw std::runtime_error(std::string(operation) + ": " + cudaGetErrorString(status));
}
} // namespace
void assembleSubdomainSamples(const SubdomainPanel *panels, const SubdomainPairSample *samples,
                              std::size_t count, std::size_t panelCount, unsigned dimension,
                              double k, unsigned excitationCount, const cuDoubleComplex *knownFlux,
                              const double *scale, cuDoubleComplex *matrix, cuDoubleComplex *rhs,
                              int *error, cudaStream_t stream) {
    if (!count)
        return;
    assembleSamples<<<static_cast<unsigned>((count + 127) / 128), 128, 0, stream>>>(
        panels, samples, count, panelCount, dimension, k, excitationCount, knownFlux, scale, matrix,
        rhs, error);
    checkLaunch("Subdomain pair assembly");
}
void assembleSubdomainMass(const SubdomainPanel *panels, std::size_t count, unsigned dimension,
                           double k, unsigned excitationCount, const cuDoubleComplex *knownFlux,
                           const double *scale, cuDoubleComplex *matrix, cuDoubleComplex *rhs,
                           int *error, cudaStream_t stream) {
    if (!count)
        return;
    assembleMass<<<static_cast<unsigned>((count + 127) / 128), 128, 0, stream>>>(
        panels, count, dimension, k, excitationCount, knownFlux, scale, matrix, rhs, error);
    checkLaunch("Subdomain mass assembly");
}
} // namespace hornsculpt::cuda_detail

// Host-side assembly infrastructure. Every domain owns its own P1 test rows;
// equal global mesh vertices share pressure columns, while interface vertices
// receive additional P1 flux columns. Quadrature is streamed in bounded batches.
#include "../galerkin/near_projection.hpp"
#include "../galerkin/regular_qualification.hpp"
#include "../galerkin/space.hpp"
#include "../lem_coupling.hpp"
#include "../planning_frequencies.hpp"
#include <algorithm>
#include <complex>
#include <cusolverDn.h>
#include <hornsculpt/acoustic_partition.hpp>
#include <hornsculpt/bem.hpp>
#include <hornsculpt/field_evaluation.hpp>
#include <hornsculpt/meshing.hpp>
#include <hornsculpt/surface_measures.hpp>
#include <limits>
#include <map>
#include <numbers>
#include <numeric>
#include <set>
#include <vector>

namespace hornsculpt::cuda_detail {
namespace {
template <class T> class DomainDevice {
  public:
    explicit DomainDevice(std::size_t n) : n_(n) {
        if (n_ && cudaMalloc(reinterpret_cast<void **>(&p_), n_ * sizeof(T)) != cudaSuccess)
            throw std::runtime_error("Subdomain CUDA allocation failed.");
    }
    ~DomainDevice() {
        if (p_)
            cudaFree(p_);
    }
    DomainDevice(const DomainDevice &) = delete;
    T *data() const {
        return p_;
    }
    void upload(const std::vector<T> &v) {
        if (v.size() > n_ ||
            cudaMemcpy(p_, v.data(), v.size() * sizeof(T), cudaMemcpyHostToDevice) != cudaSuccess)
            throw std::runtime_error("Subdomain CUDA upload failed.");
    }
    std::vector<T> download() const {
        std::vector<T> v(n_);
        if (cudaMemcpy(v.data(), p_, n_ * sizeof(T), cudaMemcpyDeviceToHost) != cudaSuccess)
            throw std::runtime_error("Subdomain CUDA download failed.");
        return v;
    }
    void zero() {
        if (cudaMemset(p_, 0, n_ * sizeof(T)) != cudaSuccess)
            throw std::runtime_error("Subdomain CUDA clear failed.");
    }

  private:
    T *p_{};
    std::size_t n_{};
};
struct DomainHost {
    const BemSubdomain *domain{};
    SurfaceMesh mesh;
    galerkin::Space space;
    std::vector<std::uint32_t> globalVertices;
    std::vector<std::uint32_t> globalTriangles;
    unsigned rowOffset{}, panelOffset{};
};
} // namespace

namespace {
BemResult solveCoupledSubdomains(const SurfaceMesh &, Vec3, const BemRequest &, std::stop_token,
                                 const std::function<void(const BemProgress &)> &);
}
BemResult solveSubdomains(const SurfaceMesh &mesh, Vec3 worldOrigin, const BemRequest &request,
                          std::stop_token stop, std::function<void(const BemProgress &)> progress,
                          bool allowZero) {
    if (request.coupling)
        return solveCoupledSubdomains(mesh, worldOrigin, request, stop, progress);
    if (stop.stop_requested())
        throw std::runtime_error("Subdomain solve cancelled.");
    if (progress)
        progress({0, request.frequencies.size(), 0, BemSolvePhase::preparing});
    if (!request.subdomains || request.infiniteBaffle || request.axisymmetric ||
        request.numerics.acceleration != BemAcceleration::dense)
        throw std::invalid_argument("Subdomain BEM requires 3D dense mode without baffle or LEM.");
    const auto &model = *request.subdomains;
    validateAcousticPartition(mesh, model);
    if (model.domains.size() < 2 || model.interfaces.empty() || request.excitations.empty() ||
        request.frequencies.empty() || !std::isfinite(request.medium.density) ||
        !(request.medium.density > 0) || !std::isfinite(request.medium.soundSpeed) ||
        !(request.medium.soundSpeed > 0) ||
        (request.numerics.quadratureOrder != 2 && request.numerics.quadratureOrder != 4) ||
        request.numerics.singularOrder < 2 || request.numerics.singularOrder > 20 ||
        !std::isfinite(request.numerics.nearDistanceRatio) ||
        request.numerics.nearDistanceRatio < 0 || request.numerics.nearDistanceRatio > 10 ||
        !std::isfinite(request.numerics.nearRelativeTolerance) ||
        request.numerics.nearRelativeTolerance < 1e-10 ||
        request.numerics.nearRelativeTolerance > .01 || request.numerics.nearMaximumDepth > 10 ||
        !std::isfinite(request.numerics.relativeResidualTolerance) ||
        request.numerics.relativeResidualTolerance < 1e-7 ||
        request.numerics.relativeResidualTolerance > .01)
        throw std::invalid_argument("Invalid subdomain BEM request.");
    (void)bem_detail::planningWavenumber(request.frequencies, request.interactionPlanFrequencies,
                                         request.medium.soundSpeed);
    std::set<double> frequencies;
    for (const auto frequency : request.frequencies)
        if (!std::isfinite(frequency) || frequency <= 0 || !frequencies.insert(frequency).second ||
            !std::isfinite(2 * std::numbers::pi * frequency / request.medium.soundSpeed))
            throw std::invalid_argument(
                "Subdomain frequencies must be finite, positive and unique.");
    if (mesh.vertices.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
        request.excitations.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        throw std::invalid_argument("Subdomain system exceeds cuSOLVER index capacity.");
    std::map<std::uint32_t, std::size_t> triangleInterface;
    std::map<std::uint32_t, unsigned> interfaceNode;
    std::vector<std::vector<std::uint32_t>> interfaceVertices(model.interfaces.size());
    for (std::size_t n = 0; n < model.interfaces.size(); ++n) {
        std::set<std::uint32_t> vertices;
        for (auto t : model.interfaces[n].triangleIndices) {
            triangleInterface.emplace(t, n);
            for (auto v : mesh.triangles[t].vertices)
                vertices.insert(v);
        }
        for (auto v : vertices) {
            if (!interfaceNode.emplace(v, static_cast<unsigned>(interfaceNode.size())).second)
                throw std::invalid_argument(
                    "Junctions shared by multiple interfaces are not supported.");
            interfaceVertices[n].push_back(v);
        }
    }
    std::vector<DomainHost> domains;
    unsigned rowCount = 0, panelCount = 0;
    for (const auto &domain : model.domains) {
        DomainHost h;
        h.domain = &domain;
        h.rowOffset = rowCount;
        h.panelOffset = panelCount;
        std::map<std::uint32_t, std::uint32_t> local;
        for (const auto &entry : domain.panels) {
            auto tri = mesh.triangles[entry.triangleIndex];
            if (entry.normalSign < 0)
                std::swap(tri.vertices[1], tri.vertices[2]);
            for (auto &v : tri.vertices) {
                if (v >= mesh.vertices.size())
                    throw std::invalid_argument("Invalid BEM vertex index.");
                const auto [pos, inserted] =
                    local.emplace(v, static_cast<std::uint32_t>(local.size()));
                if (inserted) {
                    h.globalVertices.push_back(v);
                    h.mesh.vertices.push_back(mesh.vertices[v]);
                }
                v = pos->second;
            }
            h.mesh.triangles.push_back(tri);
            h.globalTriangles.push_back(entry.triangleIndex);
        }
        h.space = galerkin::buildSpace(h.mesh);
        const auto intMax = static_cast<std::size_t>(std::numeric_limits<int>::max());
        if (h.space.vertices.size() > intMax - rowCount ||
            h.space.panels.size() > intMax - panelCount)
            throw std::invalid_argument("Subdomain panel or row count exceeds index capacity.");
        rowCount += static_cast<unsigned>(h.space.vertices.size());
        panelCount += static_cast<unsigned>(h.space.panels.size());
        domains.push_back(std::move(h));
    }
    if (interfaceNode.size() >
        static_cast<std::size_t>(std::numeric_limits<int>::max()) - mesh.vertices.size())
        throw std::invalid_argument("Subdomain unknown count exceeds index capacity.");
    const auto dimension = static_cast<unsigned>(mesh.vertices.size() + interfaceNode.size());
    if (dimension != rowCount || !dimension)
        throw std::invalid_argument(
            "Subdomain P1 interface topology has unsupported shared vertices.");
    std::vector<SubdomainPanel> panels(panelCount);
    std::vector<double> scales(rowCount);
    for (const auto &h : domains) {
        for (unsigned v = 0; v < h.space.vertices.size(); ++v)
            scales[h.rowOffset + v] = 1 / h.space.lumpedMass[v];
        for (unsigned t = 0; t < h.space.panels.size(); ++t) {
            const auto &src = h.space.panels[t];
            auto &dst = panels[h.panelOffset + t];
            dst.panel.centroid = src.centroid;
            dst.panel.normal = src.normal;
            dst.panel.area = src.area;
            dst.panel.edge = src.maximumEdgeLength;
            const auto found = triangleInterface.find(h.globalTriangles[t]);
            const bool isInterface = found != triangleInterface.end();
            dst.fluxSign =
                isInterface && model.interfaces[found->second].firstDomainId != h.domain->id ? -1
                                                                                             : 1;
            for (unsigned j = 0; j < 3; ++j) {
                dst.panel.vertices[j] = src.vertices[j];
                dst.panel.nodes[j] = src.nodes[j];
                dst.panel.curls[j] = src.curls[j];
                dst.rows[j] = static_cast<int>(h.rowOffset + src.nodes[j]);
                dst.pressureColumns[j] = static_cast<int>(h.globalVertices[src.nodes[j]]);
                dst.fluxColumns[j] =
                    isInterface ? static_cast<int>(mesh.vertices.size() +
                                                   interfaceNode.at(h.globalVertices[src.nodes[j]]))
                                : -1;
            }
        }
    }
    const unsigned nrhs = static_cast<unsigned>(request.excitations.size());
    const auto checkedMultiply = [](std::size_t a, std::size_t b) {
        if (a && b > std::numeric_limits<std::size_t>::max() / a)
            throw std::invalid_argument("Subdomain allocation size overflows.");
        return a * b;
    };
    const auto checkedAdd = [](std::size_t a, std::size_t b) {
        if (b > std::numeric_limits<std::size_t>::max() - a)
            throw std::invalid_argument("Subdomain allocation size overflows.");
        return a + b;
    };
    constexpr std::size_t batchCapacity = 32768;
    auto bytes = checkedMultiply(sizeof(SubdomainPanel), panelCount);
    bytes = checkedAdd(bytes, checkedMultiply(sizeof(double), scales.size()));
    bytes = checkedAdd(
        bytes, checkedMultiply(sizeof(cuDoubleComplex), checkedMultiply(dimension, dimension)));
    bytes = checkedAdd(bytes,
                       checkedMultiply(sizeof(cuDoubleComplex), checkedMultiply(dimension, nrhs)));
    bytes = checkedAdd(bytes,
                       checkedMultiply(sizeof(cuDoubleComplex), checkedMultiply(panelCount, nrhs)));
    bytes = checkedAdd(bytes, checkedMultiply(sizeof(SubdomainPairSample), batchCapacity));
    bytes =
        checkedAdd(bytes, checkedMultiply(sizeof(int), static_cast<std::size_t>(dimension) + 2));
    std::size_t freeBytes{}, totalBytes{};
    if (cudaMemGetInfo(&freeBytes, &totalBytes) != cudaSuccess)
        throw std::runtime_error("Cannot query GPU memory for subdomain BEM.");
    if (bytes > freeBytes)
        throw std::runtime_error("Insufficient GPU memory for subdomain BEM buffers.");
    DomainDevice<SubdomainPanel> dp(panelCount);
    dp.upload(panels);
    DomainDevice<double> ds(scales.size());
    ds.upload(scales);
    DomainDevice<cuDoubleComplex> matrix(static_cast<std::size_t>(dimension) * dimension),
        rhs(static_cast<std::size_t>(dimension) * nrhs);
    DomainDevice<cuDoubleComplex> dq(static_cast<std::size_t>(panelCount) * nrhs);
    DomainDevice<SubdomainPairSample> db(batchCapacity);
    DomainDevice<int> error(1), info(1), pivots(dimension);
    cusolverDnHandle_t solver{};
    if (cusolverDnCreate(&solver) != CUSOLVER_STATUS_SUCCESS)
        throw std::runtime_error("Subdomain cuSOLVER creation failed.");
    struct HandleGuard {
        cusolverDnHandle_t handle;
        ~HandleGuard() {
            cusolverDnDestroy(handle);
        }
    } guard{solver};
    BemResult result;
    result.triangleCount = mesh.triangles.size();
    result.pressureDofCount = mesh.vertices.size();
    result.device = CudaBemSolver::deviceDescription();
    result.acceleration = BemAcceleration::dense;
    galerkin::PairQuadratureSettings qs;
    qs.regularOrder = request.numerics.quadratureOrder;
    qs.singularOrder = request.numerics.singularOrder;
    qs.nearDistanceRatio = request.numerics.nearDistanceRatio;
    qs.relativeTolerance = request.numerics.nearRelativeTolerance;
    qs.maximumSubdivisionDepth = request.numerics.nearMaximumDepth;
    const auto low = galerkin::triangleRule(8), high = galerkin::triangleRule(14);
    using Clock = std::chrono::steady_clock;
    const auto elapsed = [](Clock::time_point begin) {
        return std::chrono::duration<double, std::milli>(Clock::now() - begin).count();
    };
    std::size_t completed = 0;
    for (double frequency : request.frequencies) {
        const auto assemblyStart = Clock::now();
        if (progress)
            progress({completed, request.frequencies.size(), frequency,
                      BemSolvePhase::nearConstruction});
        if (stop.stop_requested())
            throw std::runtime_error("Subdomain solve cancelled.");
        if (!(frequency > 0) || !std::isfinite(frequency))
            throw std::invalid_argument("Invalid BEM frequency.");
        const double omega = 2 * std::numbers::pi * frequency,
                     k = omega / request.medium.soundSpeed;
        std::vector<cuDoubleComplex> q(static_cast<std::size_t>(panelCount) * nrhs);
        for (unsigned e = 0; e < nrhs; ++e) {
            std::map<std::uint32_t, std::complex<double>> values;
            for (const auto &v : request.excitations[e].faces) {
                if (!std::isfinite(v.normalVelocity.real()) ||
                    !std::isfinite(v.normalVelocity.imag()) ||
                    !values.emplace(v.faceId, v.normalVelocity).second)
                    throw std::invalid_argument("Invalid subdomain excitation.");
            }
            std::set<std::uint32_t> usedFaces;
            bool nonzero = false;
            for (const auto &h : domains)
                for (unsigned t = 0; t < h.space.panels.size(); ++t) {
                    const auto g = h.globalTriangles[t];
                    if (triangleInterface.contains(g))
                        continue;
                    const auto &triangle = mesh.triangles[g];
                    const auto it = values.find(triangle.faceId);
                    if (it == values.end())
                        continue;
                    if (triangle.boundary != BoundaryCondition::driven)
                        throw std::invalid_argument(
                            "Subdomain excitation references a non-driven face.");
                    usedFaces.insert(triangle.faceId);
                    nonzero = nonzero || std::abs(it->second) > 0;
                    const auto value =
                        std::complex<double>{0, request.medium.density * omega} * it->second;
                    q[h.panelOffset + t + static_cast<std::size_t>(e) * panelCount] =
                        make_cuDoubleComplex(value.real(), value.imag());
                }
            if (usedFaces.size() != values.size() || usedFaces.empty() || (!allowZero && !nonzero))
                throw std::invalid_argument(
                    "Subdomain excitation contains unknown faces or is empty.");
        }
        dq.upload(q);
        matrix.zero();
        rhs.zero();
        error.zero();
        std::vector<SubdomainPairSample> batch;
        batch.reserve(batchCapacity);
        std::vector<std::vector<galerkin::RegularPanelBounds>> bounds;
        unsigned minimumDegree = 4;
        for (const auto &h : domains) {
            auto &domainBounds = bounds.emplace_back();
            for (const auto &panel : h.space.panels) {
                domainBounds.push_back(galerkin::makeRegularPanelBounds(panel));
                minimumDegree = std::max(
                    minimumDegree,
                    galerkin::minimumEndpointDegree(k * panel.maximumEdgeLength, qs.regularOrder));
            }
        }
        const auto regular4 = galerkin::triangleRule(4), regular6 = galerkin::triangleRule(6),
                   regular8 = galerkin::triangleRule(8);
        const auto flush = [&]() {
            if (batch.empty())
                return;
            db.upload(batch);
            assembleSubdomainSamples(dp.data(), db.data(), batch.size(), panelCount, dimension, k,
                                     nrhs, dq.data(), ds.data(), matrix.data(), rhs.data(),
                                     error.data(), nullptr);
            if (cudaDeviceSynchronize() != cudaSuccess || error.download()[0])
                throw std::runtime_error("Subdomain quadrature assembly failed.");
            batch.clear();
        };
        const auto push = [&](unsigned test, unsigned trial,
                              const galerkin::PairQuadraturePoint &point) {
            batch.push_back({test, trial, point});
            if (batch.size() == batchCapacity)
                flush();
        };
        for (std::size_t domainIndex = 0; domainIndex < domains.size(); ++domainIndex) {
            const auto &h = domains[domainIndex];
            for (unsigned i = 0; i < h.space.panels.size(); ++i)
                for (unsigned j = 0; j < h.space.panels.size(); ++j) {
                    if (stop.stop_requested())
                        throw std::runtime_error("Subdomain solve cancelled.");
                    const auto &test = h.space.panels[i], &trial = h.space.panels[j];
                    const auto topology = galerkin::classifyPair(test, trial).topology;
                    const auto degree = topology == galerkin::PairTopology::disjoint
                                            ? galerkin::qualifiedRegularDegree(
                                                  bounds[domainIndex][i], bounds[domainIndex][j], k,
                                                  qs.relativeTolerance, minimumDegree)
                                            : 0;
                    const bool near = topology == galerkin::PairTopology::disjoint &&
                                      galerkin::isNearPair(test, trial, qs.nearDistanceRatio);
                    if (degree && !near) {
                        const auto &rule = degree == 4   ? regular4
                                           : degree == 6 ? regular6
                                                         : regular8;
                        for (const auto &x : rule)
                            for (const auto &y : rule) {
                                galerkin::PairQuadraturePoint point{};
                                for (unsigned z = 0; z < 3; ++z) {
                                    point.test[z] = x.barycentric[z];
                                    point.trial[z] = y.barycentric[z];
                                }
                                point.weight = x.weight * y.weight;
                                push(h.panelOffset + i, h.panelOffset + j, point);
                            }
                        continue;
                    }
                    const auto rule = galerkin::makePairQuadrature(
                        test, trial, qs, k, stop,
                        topology == galerkin::PairTopology::disjoint && !degree);
                    for (const auto &point : rule.points)
                        push(h.panelOffset + i, h.panelOffset + j, point);
                    for (const auto &patch : rule.patches) {
                        const auto &tr = patch.degree == 14 ? high : low;
                        const auto count =
                            galerkin::nearPatchPointCount(patch, static_cast<unsigned>(tr.size()));
                        for (unsigned x = 0; x < count; ++x)
                            push(h.panelOffset + i, h.panelOffset + j,
                                 galerkin::mapNearPatchPoint(
                                     tr.data(), static_cast<unsigned>(tr.size()), patch, x,
                                     panels[h.panelOffset + i].panel,
                                     panels[h.panelOffset + j].panel));
                    }
                }
        }
        flush();
        assembleSubdomainMass(dp.data(), panelCount, dimension, k, nrhs, dq.data(), ds.data(),
                              matrix.data(), rhs.data(), error.data(), nullptr);
        if (cudaDeviceSynchronize() != cudaSuccess || error.download()[0])
            throw std::runtime_error("Subdomain mass assembly failed.");
        const double assemblyMs = elapsed(assemblyStart);
        if (stop.stop_requested())
            throw std::runtime_error("Subdomain solve cancelled.");
        const auto original = matrix.download(), right = rhs.download();
        int workspaceSize{};
        if (cusolverDnZgetrf_bufferSize(solver, dimension, dimension, matrix.data(), dimension,
                                        &workspaceSize) != CUSOLVER_STATUS_SUCCESS ||
            workspaceSize <= 0)
            throw std::runtime_error("Subdomain LU workspace failed.");
        if (cudaMemGetInfo(&freeBytes, &totalBytes) != cudaSuccess ||
            checkedMultiply(sizeof(cuDoubleComplex), static_cast<std::size_t>(workspaceSize)) >
                freeBytes)
            throw std::runtime_error("Insufficient GPU memory for subdomain LU workspace.");
        DomainDevice<cuDoubleComplex> workspace(static_cast<std::size_t>(workspaceSize));
        const auto ownedBytes =
            checkedAdd(bytes, checkedMultiply(sizeof(cuDoubleComplex),
                                              static_cast<std::size_t>(workspaceSize)));
        result.peakDeviceBytes = std::max(result.peakDeviceBytes, ownedBytes);
        const auto factorStart = Clock::now();
        if (cusolverDnZgetrf(solver, dimension, dimension, matrix.data(), dimension,
                             workspace.data(), pivots.data(),
                             info.data()) != CUSOLVER_STATUS_SUCCESS ||
            cudaDeviceSynchronize() != cudaSuccess || info.download()[0] != 0)
            throw std::runtime_error("Subdomain LU factorization failed.");
        if (stop.stop_requested())
            throw std::runtime_error("Subdomain solve cancelled.");
        const double factorMs = elapsed(factorStart);
        const auto solveStart = Clock::now();
        if (cusolverDnZgetrs(solver, CUBLAS_OP_N, dimension, nrhs, matrix.data(), dimension,
                             pivots.data(), rhs.data(), dimension,
                             info.data()) != CUSOLVER_STATUS_SUCCESS ||
            cudaDeviceSynchronize() != cudaSuccess || info.download()[0] != 0)
            throw std::runtime_error("Subdomain LU solve failed.");
        if (stop.stop_requested())
            throw std::runtime_error("Subdomain solve cancelled.");
        const double solveMs = elapsed(solveStart);
        const auto solved = rhs.download();
        BemFrequencyResult output;
        output.frequency = frequency;
        output.assemblyMilliseconds = assemblyMs;
        output.factorizationMilliseconds = factorMs;
        output.triangularSolveMilliseconds = solveMs;
        output.solveMilliseconds = factorMs + solveMs;
        if (progress)
            progress(
                {completed, request.frequencies.size(), frequency, BemSolvePhase::residualAudit});
        for (unsigned e = 0; e < nrhs; ++e) {
            BemExcitationResult excitation;
            excitation.name = request.excitations[e].name;
            const auto residualStart = Clock::now();
            double residual = 0, norm = 0;
            for (unsigned i = 0; i < dimension; ++i) {
                if (stop.stop_requested())
                    throw std::runtime_error("Subdomain solve cancelled.");
                std::complex<double> sum{};
                for (unsigned j = 0; j < dimension; ++j) {
                    const auto a = original[i + static_cast<std::size_t>(j) * dimension],
                               x = solved[j + static_cast<std::size_t>(e) * dimension];
                    sum += std::complex<double>{a.x, a.y} * std::complex<double>{x.x, x.y};
                }
                const auto b = right[i + static_cast<std::size_t>(e) * dimension];
                residual += std::norm((sum - std::complex<double>{b.x, b.y}) / scales[i]);
                norm += (b.x * b.x + b.y * b.y) / (scales[i] * scales[i]);
            }
            output.residualCheckMilliseconds += elapsed(residualStart);
            excitation.boundaryResidual = std::sqrt(residual / std::max(norm, 1e-300));
            if (!std::isfinite(excitation.boundaryResidual) ||
                excitation.boundaryResidual > request.numerics.relativeResidualTolerance)
                throw std::runtime_error("Subdomain direct residual exceeds tolerance.");
            std::vector<BemSubdomainTrace> domainTraces;
            for (const auto &h : domains) {
                BemSubdomainTrace trace;
                trace.domainId = h.domain->id;
                trace.vertexIndices = h.globalVertices;
                for (auto v : h.globalVertices) {
                    const auto z = solved[v + static_cast<std::size_t>(e) * dimension];
                    trace.pressureTrace.emplace_back(z.x, z.y);
                }
                trace.triangleIndices = h.globalTriangles;
                for (unsigned t = 0; t < h.space.panels.size(); ++t) {
                    const auto &panel = panels[h.panelOffset + t];
                    std::complex<double> derivative{};
                    if (panel.fluxColumns[0] >= 0)
                        for (unsigned j = 0; j < 3; ++j) {
                            const auto z = solved[panel.fluxColumns[j] +
                                                  static_cast<std::size_t>(e) * dimension];
                            derivative += panel.fluxSign * k * std::complex<double>{z.x, z.y} / 3.;
                        }
                    else {
                        const auto z =
                            q[h.panelOffset + t + static_cast<std::size_t>(e) * panelCount];
                        derivative = {z.x, z.y};
                    }
                    trace.normalDerivative.push_back(derivative);
                }
                for (std::size_t i = 0; i < model.interfaces.size(); ++i) {
                    const auto &it = model.interfaces[i];
                    if (it.firstDomainId != h.domain->id && it.secondDomainId != h.domain->id)
                        continue;
                    BemSubdomainTrace::InterfaceFlux flux;
                    flux.interfaceId = it.id;
                    flux.vertexIndices = interfaceVertices[i];
                    const double sign = it.firstDomainId == h.domain->id ? 1 : -1;
                    for (auto v : flux.vertexIndices) {
                        const auto z = solved[mesh.vertices.size() + interfaceNode.at(v) +
                                              static_cast<std::size_t>(e) * dimension];
                        flux.normalDerivative.emplace_back(sign * k * z.x, sign * k * z.y);
                    }
                    trace.interfaceFluxes.push_back(std::move(flux));
                }
                domainTraces.push_back(std::move(trace));
            }
            excitation.setSubdomainTraces(std::move(domainTraces));
            if (progress)
                progress({completed, request.frequencies.size(), frequency,
                          BemSolvePhase::fieldEvaluation});
            const auto fieldStart = Clock::now();
            // Finite field uses the shared adaptive evaluator, including its
            // geometric boundary exclusion and convergence mask.
            if (!request.farField && !request.observers.empty()) {
                auto worldMesh = mesh;
                for (auto &v : worldMesh.vertices)
                    v = add(v, worldOrigin);
                PressureFieldSettings settings;
                settings.relativeQuadratureTolerance =
                    std::min(request.numerics.nearRelativeTolerance,
                             .025 * request.numerics.relativeResidualTolerance);
                settings.maximumQuadratureDepth = request.numerics.nearMaximumDepth;
                const auto field = evaluatePressureField(worldMesh, request, frequency, excitation,
                                                         request.observers, settings, stop);
                for (const auto &sample : field) {
                    if (sample.status != FieldPointStatus::air)
                        throw std::runtime_error(
                            "Subdomain observation is excluded or quadrature did not converge.");
                    excitation.observerPressure.push_back(sample.pressure);
                }
            }
            // Compare adjacent quadrature orders against the same complex
            // observer-vector criterion as the single-domain CUDA path.
            if (request.farField && !request.observers.empty()) {
                const DomainHost *selected = nullptr;
                for (const auto &h : domains)
                    if (h.domain->kind == BemDomainKind::exterior) {
                        if (selected)
                            throw std::invalid_argument("Far field requires one exterior domain.");
                        selected = &h;
                    }
                if (!selected)
                    throw std::invalid_argument("Far field requires an exterior domain.");
                const auto origin = sub(request.farField->origin, worldOrigin);
                const auto distance = request.farField->referenceDistance;
                if (!std::isfinite(distance) || !(distance > 0))
                    throw std::invalid_argument("Invalid far-field reference distance.");
                std::vector<Vec3> directions;
                directions.reserve(request.observers.size());
                for (const auto worldPoint : request.observers) {
                    auto direction = sub(sub(worldPoint, worldOrigin), origin);
                    const auto directionLength = length(direction);
                    if (!std::isfinite(directionLength) || !(directionLength > 0))
                        throw std::invalid_argument("Invalid far-field observation.");
                    directions.push_back(mul(direction, 1 / directionLength));
                }
                const auto evaluate = [&](unsigned degree) {
                    const auto rule = galerkin::triangleRule(degree);
                    std::vector<std::complex<double>> values(directions.size());
                    for (std::size_t observer = 0; observer < directions.size(); ++observer) {
                        const auto direction = directions[observer];
                        for (unsigned t = 0; t < selected->space.panels.size(); ++t) {
                            if (stop.stop_requested())
                                throw std::runtime_error("Subdomain solve cancelled.");
                            const auto &panel = selected->space.panels[t];
                            const auto &sp = panels[selected->panelOffset + t];
                            for (const auto &sample : rule) {
                                Vec3 y{};
                                std::complex<double> pressure{}, flux{};
                                for (unsigned j = 0; j < 3; ++j) {
                                    y = add(y, mul(panel.vertices[j], sample.barycentric[j]));
                                    const auto p = solved[sp.pressureColumns[j] +
                                                          static_cast<std::size_t>(e) * dimension];
                                    pressure +=
                                        sample.barycentric[j] * std::complex<double>{p.x, p.y};
                                    if (sp.fluxColumns[j] >= 0) {
                                        const auto v =
                                            solved[sp.fluxColumns[j] +
                                                   static_cast<std::size_t>(e) * dimension];
                                        flux += sample.barycentric[j] * sp.fluxSign * k *
                                                std::complex<double>{v.x, v.y};
                                    }
                                }
                                if (sp.fluxColumns[0] < 0) {
                                    const auto v = q[selected->panelOffset + t +
                                                     static_cast<std::size_t>(e) * panelCount];
                                    flux = {v.x, v.y};
                                }
                                const auto kernel =
                                    std::polar(1.0 / (4 * std::numbers::pi * distance),
                                               -k * distance + k * dot(direction, sub(y, origin)));
                                const auto derivative =
                                    std::complex<double>{0, k * dot(direction, panel.normal)};
                                values[observer] += kernel * (flux - derivative * pressure) *
                                                    (2 * panel.area * sample.weight);
                            }
                        }
                    }
                    return values;
                };
                unsigned degree = minimumDegree;
                auto observed = evaluate(degree);
                for (;;) {
                    const auto next = galerkin::refinedEndpointDegree(degree);
                    auto refined = evaluate(next);
                    double differenceSquared = 0, normSquared = 0;
                    for (std::size_t i = 0; i < refined.size(); ++i) {
                        differenceSquared += std::norm(refined[i] - observed[i]);
                        normSquared += std::norm(refined[i]);
                    }
                    const double indicator =
                        normSquared > 0          ? std::sqrt(differenceSquared / normSquared)
                        : differenceSquared == 0 ? 0
                                                 : std::numeric_limits<double>::infinity();
                    excitation.fieldQuadratureComparisonDegree = degree;
                    excitation.fieldQuadratureVerificationDegree = next;
                    excitation.fieldQuadratureDiscrepancy = indicator;
                    observed = std::move(refined);
                    if (std::isfinite(indicator) &&
                        indicator <= .025 * request.numerics.relativeResidualTolerance)
                        break;
                    if (degree == 12)
                        throw std::runtime_error("Subdomain complex observer-field quadrature "
                                                 "12/14 did not converge; refine the mesh.");
                    degree = next;
                }
                excitation.observerPressure = std::move(observed);
            }
            output.fieldMilliseconds += elapsed(fieldStart);
            excitation.solvedVelocities = request.excitations[e].faces;
            for (const auto &face : request.excitations[e].faces) {
                double area = 0;
                std::complex<double> integratedPressure{};
                for (const auto &h : domains)
                    for (unsigned t = 0; t < h.space.panels.size(); ++t) {
                        const auto global = h.globalTriangles[t];
                        if (triangleInterface.contains(global) ||
                            mesh.triangles[global].faceId != face.faceId)
                            continue;
                        const auto &panel = h.space.panels[t];
                        std::complex<double> mean{};
                        for (unsigned j = 0; j < 3; ++j) {
                            const auto z = solved[h.globalVertices[panel.nodes[j]] +
                                                  static_cast<std::size_t>(e) * dimension];
                            mean += std::complex<double>{z.x, z.y} / 3.;
                        }
                        area += panel.area;
                        integratedPressure += panel.area * mean;
                    }
                const auto volume = area * face.normalVelocity;
                if (area > 0 && std::abs(volume) > 0) {
                    const auto mean = integratedPressure / area;
                    excitation.impedances.push_back(
                        {face.faceId, area, mean, volume, mean / volume});
                }
            }
            output.excitations.push_back(std::move(excitation));
        }
        result.frequencies.push_back(std::move(output));
        ++completed;
        if (progress)
            progress({completed, request.frequencies.size(), frequency,
                      BemSolvePhase::frequencyComplete});
    }
    return result;
}

namespace {
BemResult solveCoupledSubdomains(const SurfaceMesh &mesh, Vec3 origin, const BemRequest &request,
                                 std::stop_token stop,
                                 const std::function<void(const BemProgress &)> &progress) {
    if (!request.subdomains || request.infiniteBaffle || request.axisymmetric ||
        request.numerics.acceleration != BemAcceleration::dense || !request.excitations.empty())
        throw std::invalid_argument(
            "Coupled subdomains require dense 3D, ports and no independent excitation.");
    if (stop.stop_requested())
        throw std::runtime_error("Subdomain coupling cancelled.");
    const auto &model = *request.subdomains;
    validateAcousticPartition(mesh, model);
    (void)bem_detail::planningWavenumber(request.frequencies, request.interactionPlanFrequencies,
                                         request.medium.soundSpeed);
    std::set<double> frequencies;
    for (auto f : request.frequencies)
        if (!frequencies.insert(f).second)
            throw std::invalid_argument("Coupled subdomain frequencies must be unique.");
    const auto measures = surfaceMeasures(mesh);
    const auto basis = bem_detail::makeCouplingBasis(*request.coupling, mesh, measures.areas);
    const auto ports = lem::bemPorts(*request.coupling);
    std::set<std::uint32_t> interfaceTriangles, interfaceFaces;
    for (const auto &interface : model.interfaces)
        for (auto t : interface.triangleIndices) {
            interfaceTriangles.insert(t);
            interfaceFaces.insert(mesh.triangles[t].faceId);
        }
    std::map<std::uint32_t, std::set<std::uint32_t>> faceDomains;
    std::map<std::uint32_t, std::uint32_t> physicalOwner;
    for (const auto &domain : model.domains)
        for (const auto &entry : domain.panels) {
            if (interfaceTriangles.contains(entry.triangleIndex))
                continue;
            const auto face = mesh.triangles[entry.triangleIndex].faceId;
            faceDomains[face].insert(domain.id);
            physicalOwner.emplace(entry.triangleIndex, domain.id);
        }
    for (const auto *element : ports) {
        const auto &port = std::get<lem::BemPort>(element->component);
        if (port.kind != lem::BemPortKind::drivenSurface)
            throw std::invalid_argument(
                "LEM interface ports are not supported with BEM subdomains.");
        std::optional<std::uint32_t> owner;
        for (auto id : port.faceIds) {
            if (interfaceFaces.contains(id) || !faceDomains.contains(id) ||
                faceDomains.at(id).size() != 1)
                throw std::invalid_argument(
                    "A subdomain LEM port needs physical faces in one unambiguous domain.");
            const auto domain = *faceDomains.at(id).begin();
            if (owner && *owner != domain)
                throw std::invalid_argument("One LEM port cannot span BEM subdomains.");
            owner = domain;
        }
    }
    const auto pressureByPanel = [&](const BemExcitationResult &excitation) {
        std::vector<std::complex<double>> result(mesh.triangles.size());
        std::map<std::uint32_t, const BemSubdomainTrace *> traces;
        for (const auto &trace : excitation.subdomainTraces())
            traces.emplace(trace.domainId, &trace);
        std::map<std::uint32_t, std::map<std::uint32_t, std::complex<double>>> nodes;
        for (const auto &[id, trace] : traces) {
            if (trace->vertexIndices.size() != trace->pressureTrace.size())
                throw std::logic_error("Invalid subdomain pressure trace.");
            for (std::size_t n = 0; n < trace->vertexIndices.size(); ++n)
                nodes[id].emplace(trace->vertexIndices[n], trace->pressureTrace[n]);
        }
        for (const auto &[triangle, owner] : physicalOwner) {
            if (basis.panelPorts[triangle] < 0)
                continue;
            for (auto vertex : mesh.triangles[triangle].vertices)
                result[triangle] += nodes.at(owner).at(vertex) / 3.;
        }
        return result;
    };
    BemResult combined;
    std::size_t completed = 0;
    const auto count = ports.size();
    for (double frequency : request.frequencies) {
        if (stop.stop_requested())
            throw std::runtime_error("Subdomain coupling cancelled.");
        if (progress)
            progress({completed, request.frequencies.size(), frequency,
                      BemSolvePhase::nearConstruction});
        auto basisRequest = request;
        basisRequest.coupling.reset();
        basisRequest.frequencies = {frequency};
        basisRequest.excitations = basis.excitations;
        basisRequest.observers.clear();
        basisRequest.farField.reset();
        auto basisResult = solveSubdomains(mesh, origin, basisRequest, stop, {}, false);
        if (basisResult.frequencies.size() != 1 ||
            basisResult.frequencies[0].excitations.size() != count)
            throw std::logic_error("BEM coupling basis result count mismatch.");
        lem::AcousticLoad load;
        load.impedance.resize(count * count);
        for (std::size_t j = 0; j < count; ++j) {
            const auto pressure = pressureByPanel(basisResult.frequencies[0].excitations[j]);
            for (std::size_t t = 0; t < mesh.triangles.size(); ++t) {
                const int port = basis.panelPorts[t];
                if (port >= 0)
                    load.impedance[static_cast<std::size_t>(port) * count + j] +=
                        pressure[t] * measures.areas[t] * basis.panelVelocities[t];
            }
        }
        auto solution = lem::solve(*request.coupling, frequency, load,
                                   {request.medium.density, request.medium.soundSpeed}, stop);
        auto finalRequest = request;
        finalRequest.coupling.reset();
        finalRequest.frequencies = {frequency};
        finalRequest.excitations = {{"LEM", bem_detail::coupledVelocities(basis, solution)}};
        auto finalResult = solveSubdomains(mesh, origin, finalRequest, stop, {}, true);
        if (finalResult.frequencies.size() != 1 ||
            finalResult.frequencies[0].excitations.size() != 1)
            throw std::logic_error("Coupled subdomain final result count mismatch.");
        auto &excitation = finalResult.frequencies[0].excitations[0];
        const auto pressure = pressureByPanel(excitation);
        bem_detail::validateCouplingClosure(basis, solution, measures.areas, pressure);
        excitation.lemSolution = std::move(solution);
        if (combined.frequencies.empty()) {
            combined.device = finalResult.device;
            combined.triangleCount = finalResult.triangleCount;
            combined.pressureDofCount = finalResult.pressureDofCount;
            combined.acceleration = finalResult.acceleration;
            combined.geometryMode = finalResult.geometryMode;
        }
        combined.peakDeviceBytes = std::max(
            {combined.peakDeviceBytes, basisResult.peakDeviceBytes, finalResult.peakDeviceBytes});
        auto &frame = finalResult.frequencies[0];
        const auto &basisFrame = basisResult.frequencies[0];
        frame.assemblyMilliseconds += basisFrame.assemblyMilliseconds;
        frame.solveMilliseconds += basisFrame.solveMilliseconds;
        frame.factorizationMilliseconds += basisFrame.factorizationMilliseconds;
        frame.triangularSolveMilliseconds += basisFrame.triangularSolveMilliseconds;
        frame.residualCheckMilliseconds += basisFrame.residualCheckMilliseconds;
        combined.frequencies.push_back(std::move(frame));
        ++completed;
        if (progress)
            progress({completed, request.frequencies.size(), frequency,
                      BemSolvePhase::frequencyComplete});
    }
    return combined;
}
} // namespace
} // namespace hornsculpt::cuda_detail
