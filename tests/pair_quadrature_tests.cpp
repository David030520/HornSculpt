#include "near_laplace_reference.hpp"
#include "projection_point_reference.hpp"
#include <galerkin/near_projection.hpp>
#include <galerkin/pair_quadrature.hpp>
#include <galerkin/space.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <thread>

using namespace hornsculpt;
using namespace hornsculpt::galerkin;
namespace {
void require(bool value, const char *message) {
    if (!value)
        throw std::runtime_error(message);
}
void close(double actual, double expected, double tolerance, const char *message) {
    if (std::abs(actual - expected) > tolerance) {
        std::cerr << std::setprecision(17) << message << ": " << actual << " versus " << expected
                  << '\n';
        throw std::runtime_error(message);
    }
}
Vec3 subtract(Vec3 a, Vec3 b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}
double length(Vec3 a) {
    return std::sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
}
Vec3 pointDifference(const Panel &test, const double *tx, const Panel &trial, const double *ty) {
    Vec3 delta{};
    // Independent translation-invariant identity using all nine vertex pairs.
    for (unsigned i = 0; i < 3; ++i)
        for (unsigned j = 0; j < 3; ++j) {
            const auto difference = subtract(test.vertices[i], trial.vertices[j]);
            delta.x += difference.x * tx[i] * ty[j];
            delta.y += difference.y * tx[i] * ty[j];
            delta.z += difference.z * tx[i] * ty[j];
        }
    return delta;
}
Panel panel(const std::array<Vec3, 3> &vertices, std::array<std::uint32_t, 3> nodes) {
    SurfaceMesh mesh;
    mesh.vertices.assign(vertices.begin(), vertices.end());
    mesh.triangles.push_back({{0, 1, 2}, 1});
    auto p = buildSpace(mesh).panels.front();
    p.nodes = nodes;
    return p;
}
template <class Visitor>
void eachPoint(const Panel &test, const Panel &trial, const PairQuadratureRule &rule,
               Visitor visitor) {
    for (const auto &p : rule.points)
        visitor(p);
    for (const auto &patch : rule.patches) {
        const auto triangle = triangleRule(patch.degree);
        for (unsigned i = 0; i < nearPatchPointCount(patch, static_cast<unsigned>(triangle.size()));
             ++i)
            visitor(mapNearPatchPoint(triangle.data(), static_cast<unsigned>(triangle.size()),
                                      patch, i, test, trial));
    }
}
std::complex<double> integral(const Panel &test, const Panel &trial, const PairQuadratureRule &rule,
                              double k = 0) {
    std::complex<double> value{};
    eachPoint(test, trial, rule, [&](const auto &q) {
        const double r = length(pointDifference(test, q.test, trial, q.trial));
        require(r > 0, "Quadrature samples a singular point");
        value += q.weight * std::exp(std::complex<double>(0, -k * r)) / r;
    });
    return value * (4 * test.area * trial.area);
}
double determinant(std::array<std::array<double, 4>, 4> matrix) {
    double result = 1;
    for (unsigned i = 0; i < 4; ++i) {
        unsigned pivot = i;
        for (unsigned j = i + 1; j < 4; ++j)
            if (std::abs(matrix[j][i]) > std::abs(matrix[pivot][i]))
                pivot = j;
        if (pivot != i) {
            std::swap(matrix[pivot], matrix[i]);
            result = -result;
        }
        result *= matrix[i][i];
        for (unsigned j = i + 1; j < 4; ++j) {
            const double factor = matrix[j][i] / matrix[i][i];
            for (unsigned k = i + 1; k < 4; ++k)
                matrix[j][k] -= factor * matrix[i][k];
        }
    }
    return result;
}
void forcedAdaptiveFarTests() {
    const auto a = panel({{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}}, {0, 1, 2});
    const auto b = panel({{{0, 0, 3}, {1, 0, 3}, {0, 1, 3}}}, {3, 4, 5});
    PairQuadratureSettings settings;
    settings.nearDistanceRatio = 0;
    settings.relativeTolerance = 1e-4;
    const auto regular = makePairQuadrature(a, b, settings, 1.0);
    const auto adaptive = makePairQuadrature(a, b, settings, 1.0, {}, true);
    settings.relativeTolerance = 1e-6;
    const auto refined = makePairQuadrature(a, b, settings, 1.0, {}, true);
    require(!regular.points.empty() && regular.patches.empty(), "Far shortcut was not exercised");
    require(adaptive.points.empty() && !adaptive.patches.empty() && !refined.patches.empty(),
            "Forced adaptive fallback was not exercised");
    PairQuadratureRule high;
    const auto reference = triangleRule(14);
    for (const auto &x : reference)
        for (const auto &y : reference) {
            PairQuadraturePoint point{};
            for (unsigned j = 0; j < 3; ++j) {
                point.test[j] = x.barycentric[j];
                point.trial[j] = y.barycentric[j];
            }
            point.weight = x.weight * y.weight;
            high.points.push_back(point);
        }
    // All P1xP1 moments of G and its three spatial derivatives. Their linear
    // combinations cover the weak single/double/adjoint layers and the
    // regularized hypersingular terms for a planar panel pair.
    using Complex = std::complex<double>;
    const auto moments = [&](const PairQuadratureRule &rule) {
        std::array<Complex, 36> result{};
        eachPoint(a, b, rule, [&](const PairQuadraturePoint &q) {
            const auto d = pointDifference(a, q.test, b, q.trial);
            const double r = length(d);
            const Complex g = std::exp(Complex{0, -r}) / r;
            const Complex gradient = -(Complex{1, r} * g) / (r * r);
            const std::array<Complex, 4> kernel{g, gradient * d.x, gradient * d.y, gradient * d.z};
            for (unsigned i = 0; i < 3; ++i)
                for (unsigned j = 0; j < 3; ++j)
                    for (unsigned c = 0; c < 4; ++c)
                        result[4 * (3 * i + j) + c] +=
                            q.weight * q.test[i] * q.trial[j] * kernel[c];
        });
        return result;
    };
    const auto expected = moments(high), coarse = moments(adaptive), fine = moments(refined);
    const auto error = [&](const auto &actual) {
        double difference = 0, norm = 0;
        for (std::size_t i = 0; i < expected.size(); ++i) {
            difference += std::norm(actual[i] - expected[i]);
            norm += std::norm(expected[i]);
        }
        return std::sqrt(difference / norm);
    };
    const auto coarseError = error(coarse), fineError = error(fine);
    require(coarseError < 1e-4 && fineError < 1e-6,
            "Adaptive far-pair P1 Green/gradient moments miss the degree-14 reference");
    require(fineError <= coarseError * 1.1 + 1e-13,
            "Tighter far-pair tolerance worsened the P1 moment error");
}
void partitionTests() {
    for (auto topology : {PairTopology::disjoint, PairTopology::shared_vertex,
                          PairTopology::shared_edge, PairTopology::coincident}) {
        const auto rule = referencePairRule(topology, 4);
        double mass = 0;
        double moments[9]{};
        for (const auto &q : rule) {
            mass += q.weight;
            double tx = 0, ty = 0;
            require(q.weight > 0, "Nonpositive transformed weight");
            for (unsigned i = 0; i < 3; ++i) {
                require(q.test[i] >= 0 && q.test[i] <= 1 && q.trial[i] >= 0 && q.trial[i] <= 1,
                        "Transformed point outside triangle");
                tx += q.test[i];
                ty += q.trial[i];
            }
            close(tx, 1, 2e-15, "Test partition of unity");
            close(ty, 1, 2e-15, "Trial partition of unity");
            for (unsigned i = 0; i < 3; ++i)
                for (unsigned j = 0; j < 3; ++j)
                    moments[3 * i + j] += q.weight * q.test[i] * q.trial[j];
        }
        close(mass, .25, 2e-14, "Pair domains do not partition the product triangle");
        for (auto m : moments)
            close(m, 1.0 / 36, 2e-14, "Transformed P1 product moment is wrong");
        for (unsigned sector = 0; sector < pairSectorCount(topology); ++sector) {
            const std::array<double, 4> parameters{.37, .43, .61, .23};
            std::array<std::array<double, 4>, 4> jacobian{};
            for (unsigned column = 0; column < 4; ++column) {
                auto lo = parameters, hi = parameters;
                lo[column] -= 1e-6;
                hi[column] += 1e-6;
                const auto a = canonicalPairTransform(topology, sector, lo[0], lo[1], lo[2], lo[3]);
                const auto b = canonicalPairTransform(topology, sector, hi[0], hi[1], hi[2], hi[3]);
                const std::array<double, 4> av{1 - a.test[0], a.test[2], 1 - a.trial[0],
                                               a.trial[2]},
                    bv{1 - b.test[0], b.test[2], 1 - b.trial[0], b.trial[2]};
                for (unsigned row = 0; row < 4; ++row)
                    jacobian[row][column] = (bv[row] - av[row]) / 2e-6;
            }
            const auto q = canonicalPairTransform(topology, sector, parameters[0], parameters[1],
                                                  parameters[2], parameters[3]);
            close(std::abs(determinant(jacobian)), q.weight, 2e-10,
                  "Four-dimensional Jacobian is incorrect");
        }
    }
    for (unsigned degree : {2, 4, 6, 8}) {
        const auto rule = triangleRule(degree);
        double mass = 0;
        for (const auto &q : rule)
            mass += q.weight;
        close(mass, .5, 3e-15, "Triangle rule measure");
        for (unsigned a = 0; a <= degree; ++a)
            for (unsigned b = 0; b + a <= degree; ++b) {
                double actual = 0;
                for (const auto &q : rule)
                    actual += q.weight * std::pow(q.barycentric[1], static_cast<int>(a)) *
                              std::pow(q.barycentric[2], static_cast<int>(b));
                const double exact =
                    std::tgamma(a + 1.0) * std::tgamma(b + 1.0) / std::tgamma(a + b + 3.0);
                close(actual, exact, 2e-15,
                      "Triangle rule lacks its documented polynomial exactness");
            }
    }
}
void singularTests() {
    const auto a = panel({{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}}, {0, 1, 2});
    // Independent closed form from integrating radial factors and the three
    // difference directions of a unit right triangle, without a numerical oracle.
    const double exact = (2 + std::sqrt(2.0)) * std::asinh(1.0) / 3;
    double previous = 1;
    for (unsigned order : {2, 4, 8, 16}) {
        PairQuadratureSettings settings;
        settings.singularOrder = order;
        const auto rule = makePairQuadrature(a, a, settings);
        const double value = integral(a, a, rule).real(), error = std::abs(value - exact);
        std::cout << "Self Laplace order " << order << ": " << value << " absolute error " << error
                  << '\n';
        require(error < previous, "Self singular quadrature did not converge");
        previous = error;
        if (order == 4)
            require(error < 2e-4, "Production singular order 4 is inaccurate");
    }
    require(previous < 2e-12, "Singular integral misses the independent analytic value");
    const auto edge = panel({{{1, 0, 0}, {0, 0, 0}, {0, 0, 1}}}, {1, 0, 3});
    const auto vertex = panel({{{0, 0, 0}, {-1, 0, 1}, {0, -1, 1}}}, {0, 4, 5});
    for (const auto &b : {edge, vertex}) {
        PairQuadratureSettings settings;
        settings.singularOrder = 12;
        const auto reference = integral(a, b, makePairQuadrature(a, b, settings), .7);
        const auto reversed = integral(b, a, makePairQuadrature(b, a, settings), .7);
        require(std::abs(reference - reversed) < 2e-9, "Singular reciprocity does not converge");
        settings.singularOrder = 4;
        const auto production = integral(a, b, makePairQuadrature(a, b, settings), .7);
        require(std::abs(production - reference) / std::abs(reference) < .003,
                "Production touching-panel rule is inaccurate");
        auto permuted = b;
        std::swap(permuted.nodes[0], permuted.nodes[2]);
        std::swap(permuted.vertices[0], permuted.vertices[2]);
        const auto reordered = integral(a, permuted, makePairQuadrature(a, permuted, settings), .7);
        require(std::abs(production - reordered) < 2e-13,
                "Shared-node permutations change the integral");
    }
    auto scaled = a;
    for (auto &p : scaled.vertices) {
        p.x *= 3;
        p.y *= 3;
        p.z *= 3;
    }
    scaled.area *= 9;
    scaled.maximumEdgeLength *= 3;
    PairQuadratureSettings settings;
    settings.singularOrder = 6;
    const auto original = integral(a, a, makePairQuadrature(a, a, settings), .7),
               rescaled =
                   integral(scaled, scaled, makePairQuadrature(scaled, scaled, settings), .7 / 3);
    require(std::abs(rescaled - 27.0 * original) < 3e-12, "Single-layer scaling is wrong");
}
void nearTests() {
    const auto a = panel({{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}}, {0, 1, 2});
    const auto b = panel({{{.13, .07, .4}, {1.13, .07, .4}, {.13, 1.07, .4}}}, {3, 4, 5});
    require(isNearPair(a, b), "Close disjoint pair not classified near");
    require(!isNearPair(a, b, 0), "Explicit diagnostic near disable was ignored");
    PairQuadratureSettings disabled;
    disabled.nearDistanceRatio = 0;
    const auto uncorrected = makePairQuadrature(a, b, disabled, 1.3);
    require(uncorrected.patches.empty() && uncorrected.points.size() == 36 &&
                uncorrected.estimatedRelativeError < 0,
            "Diagnostic rule falsely claims controlled near accuracy");
    PairQuadratureSettings loose;
    loose.relativeTolerance = 1e-3;
    loose.maximumSubdivisionDepth = 10;
    auto strict = loose;
    strict.relativeTolerance = 1e-5;
    const auto coarse = makePairQuadrature(a, b, loose, 1.3),
               fine = makePairQuadrature(a, b, strict, 1.3);
    require(!fine.patches.empty() && fine.points.empty(), "Near rule is not compact");
    require(fine.estimatedRelativeError <= strict.relativeTolerance,
            "Near rule exceeded its error budget");
    require(fine.patches.size() > coarse.patches.size(),
            "Tighter near tolerance did not refine the integration");
    require(std::abs(integral(a, b, coarse, 1.3) - integral(a, b, fine, 1.3)) /
                    std::abs(integral(a, b, fine, 1.3)) <
                1e-3,
            "Embedded near-pair error control is inaccurate");
    double mass = 0, compensation = 0;
    eachPoint(a, b, fine, [&](const auto &q) {
        const double corrected = q.weight - compensation, next = mass + corrected;
        compensation = (next - mass) - corrected;
        mass = next;
    });
    close(mass, .25, 2e-13, "Adaptive subdivision lost area");
    auto insufficient = strict;
    insufficient.maximumSubdivisionDepth = 0;
    bool rejected = false;
    try {
        makePairQuadrature(a, b, insufficient, 1.3);
    } catch (const std::runtime_error &) {
        rejected = true;
    }
    require(rejected, "Unresolved near interaction was silently accepted");
    std::cout << "Near pair: " << coarse.patches.size() << " / " << fine.patches.size()
              << " compact leaves, maximum depth " << fine.maximumDepth << ", estimate "
              << fine.estimatedRelativeError << '\n';
}

void nearCornerTests() {
    const auto a = panel({{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}}, {0, 1, 2});
    const auto b = panel({{{-.01, 0, 0}, {-1.01, 0, 0}, {-.01, -1, 0}}}, {3, 4, 5});
    PairQuadratureSettings settings;
    settings.relativeTolerance = 1e-4;
    settings.maximumSubdivisionDepth = 6;
    const auto controlled = makePairQuadrature(a, b, settings, .7);
    require(controlled.maximumDepth <= 6 &&
                controlled.estimatedRelativeError <= settings.relativeTolerance,
            "Near-corner error budget or subdivision limit changed");
    // Independent tensor Gauss rules, without the adaptive degree-6/8 patch logic.
    PairQuadratureRule reference16, reference20;
    reference16.points = referencePairRule(PairTopology::disjoint, 16);
    reference20.points = referencePairRule(PairTopology::disjoint, 20);
    const auto medium = integral(a, b, reference16, .7),
               reference = integral(a, b, reference20, .7), actual = integral(a, b, controlled, .7);
    require(std::abs(medium - reference) / std::abs(reference) < 1e-5,
            "Near-corner independent reference did not converge");
    require(std::abs(actual - reference) / std::abs(reference) < settings.relativeTolerance,
            "Global near-corner estimator accepted an inaccurate integral");
    std::stop_source cancelled;
    cancelled.request_stop();
    bool rejected = false;
    try {
        (void)makePairQuadrature(a, b, settings, .7, cancelled.get_token());
    } catch (const std::runtime_error &) {
        rejected = true;
    }
    require(rejected, "Cancelled near-pair integration continued");
    // An explicitly unresolved accuracy request must fail. The former Cartesian
    // failure at depth six is no longer an invariant of this geometry: the
    // projected integration coordinates can now resolve it at that depth.
    const auto parallel = panel({{{0, 0, .001}, {1, 0, .001}, {0, 1, .001}}}, {3, 4, 5});
    auto insufficient = settings;
    insufficient.maximumSubdivisionDepth = 0;
    insufficient.relativeTolerance = 1e-12;
    rejected = false;
    try {
        (void)makePairQuadrature(a, parallel, insufficient, .7);
    } catch (const std::runtime_error &) {
        rejected = true;
    }
    require(rejected, "Unresolved close parallel pair was accepted by weakening accuracy");
    std::stop_source during;
    std::jthread trigger([&] {
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
        during.request_stop();
    });
    rejected = false;
    try {
        (void)makePairQuadrature(a, parallel, settings, .7, during.get_token());
    } catch (const std::runtime_error &error) {
        rejected = std::string(error.what()).find("cancelled") != std::string::npos;
    }
    require(rejected, "Active near-pair refinement ignored cancellation");
    std::cout << "Near-corner controlled integral: " << controlled.patches.size()
              << " patches; independent relative error "
              << std::abs(actual - reference) / std::abs(reference) << '\n';
}

using OperatorMoments = std::array<std::complex<double>, 24>;
OperatorMoments operatorMoments(const Panel &test, const Panel &trial,
                                const PairQuadratureRule &rule, double k) {
    OperatorMoments result{};
    auto dot = [](Vec3 a, Vec3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; };
    eachPoint(test, trial, rule, [&](const auto &q) {
        const auto delta = pointDifference(test, q.test, trial, q.trial);
        const double r = length(delta);
        require(r > 0, "Translated valid geometry produced a zero quadrature distance");
        const auto green = std::exp(std::complex<double>(0, -k * r)) *
                           (q.weight * test.area * trial.area / (std::numbers::pi * r));
        const auto radial = green * std::complex<double>(1, k * r) / (r * r);
        const auto source = radial * dot(delta, trial.normal);
        const auto target = -radial * dot(delta, test.normal);
        for (unsigned i = 0; i < 3; ++i) {
            result[i] += green * q.test[i];
            result[12 + i] += target * q.test[i];
            for (unsigned j = 0; j < 3; ++j) {
                result[3 + i * 3 + j] += source * q.test[i] * q.trial[j];
                result[15 + i * 3 + j] +=
                    green * (dot(test.curls[i], trial.curls[j]) -
                             k * k * dot(test.normal, trial.normal) * q.test[i] * q.trial[j]);
            }
        }
    });
    return result;
}

PairQuadratureRule subdividedTensorReference(unsigned order) {
    PairQuadratureRule result;
    const auto base = referencePairRule(PairTopology::disjoint, order);
    // Independent uniform subdivision and tensor Gauss integration. Neither
    // the production error estimator nor its promoted triangle rule is used.
    const double children[4][3][3] = {{{1, 0, 0}, {.5, .5, 0}, {.5, 0, .5}},
                                      {{.5, .5, 0}, {0, 1, 0}, {0, .5, .5}},
                                      {{.5, 0, .5}, {0, .5, .5}, {0, 0, 1}},
                                      {{.5, .5, 0}, {0, .5, .5}, {.5, 0, .5}}};
    result.points.reserve(16 * base.size());
    for (const auto &test : children)
        for (const auto &trial : children)
            for (const auto &point : base) {
                PairQuadraturePoint mapped{};
                mapped.weight = point.weight / 16;
                for (unsigned i = 0; i < 3; ++i)
                    for (unsigned j = 0; j < 3; ++j) {
                        mapped.test[i] += point.test[j] * test[j][i];
                        mapped.trial[i] += point.trial[j] * trial[j][i];
                    }
                result.points.push_back(mapped);
            }
    return result;
}

void promotedNearTests() {
    for (const unsigned degree : {12U, 14U}) {
        const auto rule = triangleRule(degree);
        for (unsigned px = 0; px <= degree; ++px)
            for (unsigned py = 0; py + px <= degree; ++py) {
                double actual = 0;
                for (const auto &point : rule) {
                    require(point.weight > 0, "Promoted triangle rule has a nonpositive weight");
                    actual += point.weight * std::pow(point.barycentric[1], px) *
                              std::pow(point.barycentric[2], py);
                }
                const double expected =
                    std::tgamma(px + 1.) * std::tgamma(py + 1.) / std::tgamma(px + py + 3.);
                require(std::abs(actual / expected - 1) < 5e-13,
                        "Promoted triangle rule fails polynomial exactness");
            }
    }
    const auto a = panel({{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}}, {0, 1, 2});
    const auto b = panel({{{.13, .07, .1}, {1.13, .07, .1}, {.13, 1.07, .1}}}, {3, 4, 5});
    PairQuadratureSettings settings;
    settings.relativeTolerance = 1e-4;
    settings.maximumSubdivisionDepth = 6;
    const auto controlled = makePairQuadrature(a, b, settings, .7);
    require(controlled.maximumDepth <= 6 &&
                controlled.estimatedRelativeError <= settings.relativeTolerance &&
                std::any_of(controlled.patches.begin(), controlled.patches.end(),
                            [](const auto &patch) { return patch.degree == 14; }),
            "Close parallel panels did not use controlled polynomial refinement");
    const auto low = operatorMoments(a, b, subdividedTensorReference(12), .7),
               high = operatorMoments(a, b, subdividedTensorReference(16), .7),
               actual = operatorMoments(a, b, controlled, .7);
    double worst = 0;
    for (const auto [start, count] :
         std::array<std::array<unsigned, 2>, 4>{{{0, 3}, {3, 9}, {12, 3}, {15, 9}}}) {
        double error = 0, convergence = 0, reference = 0;
        for (unsigned i = start; i < start + count; ++i) {
            error += std::norm(actual[i] - high[i]);
            convergence += std::norm(low[i] - high[i]);
            reference += std::norm(high[i]);
        }
        require(std::sqrt(convergence / reference) < 1e-5,
                "Independent close-parallel operator reference did not converge");
        const double relative = std::sqrt(error / reference);
        require(relative < settings.relativeTolerance,
                "Promoted near quadrature exceeds the requested operator error");
        worst = std::max(worst, relative);
    }
    double mass = 0;
    eachPoint(a, b, controlled, [&](const auto &point) { mass += point.weight; });
    close(mass, .25, 2e-12, "Mixed-degree adaptive integration lost reference area");
    std::cout << "Close parallel panels: " << controlled.patches.size()
              << " patches, unchanged depth/tolerance 6/1e-4, independent S/D/K'/H error " << worst
              << '\n';
}

void cachedProjectionPointParity(const Panel &a, const Panel &b) {
    PairQuadraturePatch patch{};
    patch.test[0] = patch.test[4] = patch.test[8] = 1;
    patch.areaFraction = 1;
    for (unsigned child = 0; child < 2; ++child) {
        if (child) {
            const double vertices[9]{.5, .5, 0, 0, 1, 0, 0, .5, .5};
            std::copy_n(vertices, 9, patch.test);
            patch.areaFraction = .25;
        }
        for (unsigned order : {14u, 20u}) {
            const auto gauss = gaussLegendre01(order);
            for (unsigned intervals : {1u, 2u})
                for (const auto &x : triangleRule(8)) {
                    const auto projected = makeProjectionOuter(x, patch, a, b);
                    for (unsigned edge = 0; edge < 3; ++edge) {
                        const auto geometry = makeProjectionEdge(projected, edge, b);
                        for (unsigned iu = 0; iu < intervals; ++iu)
                            for (const auto &gu : gauss) {
                                const GaussPoint u{(iu + gu.position) / intervals,
                                                   gu.weight / intervals};
                                const auto ray = makeProjectionRay(projected, geometry, u);
                                for (unsigned iv = 0; iv < intervals; ++iv)
                                    for (const auto &gv : gauss) {
                                        const GaussPoint v{(iv + gv.position) / intervals,
                                                           gv.weight / intervals};
                                        const auto cached = mapProjectionRay(projected, ray, v);
                                        const auto direct =
                                            hornsculpt::testing::directProjectionPoint(
                                                x, edge, u, v, patch, a, b);
                                        for (unsigned j = 0; j < 3; ++j)
                                            require(cached.test[j] == direct.test[j] &&
                                                        cached.trial[j] == direct.trial[j],
                                                    "Projection caching moved a quadrature node.");
                                        require(cached.weight == direct.weight,
                                                "Projection caching changed a signed weight.");
                                    }
                            }
                    }
                }
        }
    }
}

void checkProjectionPair(const Panel &a, const Panel &b, const char *label) {
    cachedProjectionPointParity(a, b);
    PairQuadratureSettings settings;
    settings.relativeTolerance = 1e-4;
    settings.maximumSubdivisionDepth = 8;
    constexpr double k = 293.093;
    std::cout << "Actual Coax near pair " << label << ": starting controlled rule" << std::endl;
    const auto controlled = makePairQuadrature(a, b, settings, k);
    std::cout << "Actual Coax controlled: " << controlled.patches.size() << " patches, depth "
              << controlled.maximumDepth << ", estimate " << controlled.estimatedRelativeError
              << std::endl;
    require(!controlled.patches.empty() &&
                std::all_of(controlled.patches.begin(), controlled.patches.end(),
                            [](const auto &p) { return p.projection == 1; }) &&
                controlled.patches.size() < 4096 && controlled.maximumDepth <= 8 &&
                controlled.estimatedRelativeError <= settings.relativeTolerance,
            "Strongly near projected pair lost bounded controlled integration.");
    auto strict = settings;
    strict.relativeTolerance = 5e-6;
    const auto refined = makePairQuadrature(a, b, strict, k);
    std::cout << "Actual Coax refined: " << refined.patches.size() << " patches, estimate "
              << refined.estimatedRelativeError << std::endl;
    const auto actual = operatorMoments(a, b, controlled, k),
               reference = operatorMoments(a, b, refined, k);
    double worst{};
    for (const auto block :
         std::array<std::array<unsigned, 2>, 4>{{{0, 3}, {3, 9}, {12, 3}, {15, 9}}}) {
        double error{}, norm{};
        for (unsigned i = block[0]; i < block[0] + block[1]; ++i) {
            error += std::norm(actual[i] - reference[i]);
            norm += std::norm(reference[i]);
        }
        const double difference = std::sqrt(error / norm);
        require(difference < settings.relativeTolerance,
                "Projected near S/D/K'/H refinement exceeds its error budget.");
        worst = std::max(worst, difference);
    }
    double mass{};
    eachPoint(a, b, controlled, [&](const auto &p) { mass += p.weight; });
    close(mass, .25, 1e-6, "Signed projection fans do not cover the original triangle pair.");
    auto independentLow = hornsculpt::testing::semiAnalyticLaplacePair(a, b, 4);
    std::array<double, 36> independentHigh{};
    for (unsigned level = 5;; ++level) {
        independentHigh = hornsculpt::testing::semiAnalyticLaplacePair(a, b, level);
        double convergence{};
        for (const auto [start, count] :
             std::array<std::array<unsigned, 2>, 2>{{{0, 9}, {9, 27}}}) {
            double norm{}, difference{};
            for (unsigned i = start; i < start + count; ++i) {
                norm += independentHigh[i] * independentHigh[i];
                difference += std::pow(independentHigh[i] - independentLow[i], 2);
            }
            convergence = std::max(convergence, std::sqrt(difference / norm));
        }
        std::cout << "Independent analytic reference " << label << " level " << level
                  << ", uniform refinement difference " << convergence << std::endl;
        if (convergence < 2e-6)
            break;
        require(level < 8, "Independent semi-analytic near-pair reference did not converge.");
        independentLow = independentHigh;
    }
    std::array<double, 36> staticMoments{};
    const double scale = std::max(a.maximumEdgeLength, b.maximumEdgeLength);
    eachPoint(a, b, controlled, [&](const auto &p) {
        const auto delta = pointDifference(a, p.test, b, p.trial);
        const double radius = length(delta), derivative = -scale / (radius * radius * radius);
        const double kernels[]{1 / radius, derivative * delta.x, derivative * delta.y,
                               derivative * delta.z};
        for (unsigned i = 0; i < 3; ++i)
            for (unsigned j = 0; j < 3; ++j)
                for (unsigned c = 0; c < 4; ++c)
                    staticMoments[9 * c + 3 * i + j] +=
                        p.weight * p.test[i] * p.trial[j] * kernels[c];
    });
    double independentError{};
    for (const auto [start, count] : std::array<std::array<unsigned, 2>, 2>{{{0, 9}, {9, 27}}}) {
        double referenceNorm{}, convergence{}, error{};
        for (unsigned i = start; i < start + count; ++i) {
            referenceNorm += independentHigh[i] * independentHigh[i];
            convergence += std::pow(independentHigh[i] - independentLow[i], 2);
            error += std::pow(staticMoments[i] - independentHigh[i], 2);
        }
        require(std::sqrt(convergence / referenceNorm) < 2e-6,
                "Independent semi-analytic near-pair reference did not converge.");
        const double relative = std::sqrt(error / referenceNorm);
        require(relative < settings.relativeTolerance,
                "Projected near-pair moments fail the independent analytic Laplace reference.");
        independentError = std::max(independentError, relative);
    }
    std::cout << "Actual Coax near pair: " << controlled.patches.size() << '/'
              << refined.patches.size() << " projected patches, depth " << controlled.maximumDepth
              << ", estimate " << controlled.estimatedRelativeError
              << ", operator refinement difference " << worst << ", mass " << std::setprecision(17)
              << mass << ", independent static error " << independentError << '\n';
}

void projectionQuadratureTests() {
    for (unsigned order : {8u, 14u, 20u}) {
        const auto generated = gaussLegendre01(order);
        for (unsigned i = 0; i < order; ++i) {
            const auto stored = projection_detail::gauss(order, i);
            close(stored.position, generated[i].position, 2e-14,
                  "Shared projection Gauss position differs");
            close(stored.weight, generated[i].weight, 2e-14,
                  "Shared projection Gauss weight differs");
        }
    }
    // Actual disjoint Frame/Aperture panels from the validated rotational Coax
    // integration mesh (64 sectors, 200 meridian vertices). These immutable SI
    // fixtures exhausted 50,000 Cartesian patches at unchanged depth8/tol1e-4.
    const auto flatTest = panel({{{.07875513710298755, 0, -.00007329972059773033},
                                  {.07837590959186806, .007719353325097316, -.00007329972059773033},
                                  {.07979800834873546, 0, -.00005705693629905962}}},
                                {0, 1, 2});
    const auto flatTrial = panel({{{.08088055202128128, 0, -.0000104113930859},
                                   {.08049109005639521, .007927680417402253, -.0000104113930859},
                                   {.07953254282092659, 0, -.0000104113930859}}},
                                 {3, 4, 5});
    checkProjectionPair(flatTest, flatTrial, "1664/24704");
    // The normal dot product is -0.88818 (about28degrees from antiparallel).
    // A parallel-normal heuristic excluded this valid one-sided projection.
    const auto slopedTest =
        panel({{{.0779794520310366, .007680305663559467, -.00027914495337881377},
                {.07837590959186806, .007719353325097316, -.00007329972059773033},
                {.07875513710298755, 0, -.00007329972059773033}}},
              {0, 1, 2});
    const auto slopedTrial = panel({{{.07953254282092659, 0, -.0000104113930859},
                                     {.07914957188878863, .007795552410445549, -.0000104113930859},
                                     {.0781845336205719, 0, -.0000104113930859}}},
                                   {3, 4, 5});
    checkProjectionPair(slopedTest, slopedTrial, "1537/24576");
    // One Frame vertex lies exactly on the aperture plane; all open-interior
    // test points remain strictly below it. The panels have NO shared nodes.
    const auto boundaryTest =
        panel({{{.08518012928817169, 0, -.000022577166473138084},
                {.08476996368355154, .008349112685728839, -.000022577166473138084},
                {.0862725888227, 0, -.0000104113930859}}},
              {0, 1, 2});
    const auto boundaryTrial =
        panel({{{.08585716272682153, .008456192445229068, -.0000104113930859},
                {.08451564455921498, .008324064438272365, -.0000104113930859},
                {.08492457962234534, 0, -.0000104113930859}}},
              {3, 4, 5});
    checkProjectionPair(boundaryTest, boundaryTrial, "2432/25217");
}

Panel recentered(const Panel &source, Vec3 origin) {
    auto vertices = source.vertices;
    for (auto &vertex : vertices)
        vertex = subtract(vertex, origin);
    return panel(vertices, source.nodes);
}

void translationTests() {
    PairQuadratureSettings settings;
    settings.relativeTolerance = 1e-4;
    settings.maximumSubdivisionDepth = 6;
    const double k = 2 * std::numbers::pi * 16000 / 343;
    for (const double offset : {1e9, 1e12}) {
        const Vec3 origin{offset, -offset, offset};
        const auto a = panel({{{origin.x, origin.y, origin.z},
                               {origin.x + .008, origin.y, origin.z},
                               {origin.x, origin.y + .008, origin.z}}},
                             {0, 1, 2});
        const auto edge = panel(
            {{a.vertices[1], a.vertices[0], {origin.x, origin.y, origin.z + .008}}}, {1, 0, 3});
        const auto near = panel({{{origin.x + .0007, origin.y + .0003, origin.z + .006},
                                  {origin.x + .0007, origin.y + .0083, origin.z + .006},
                                  {origin.x + .0087, origin.y + .0003, origin.z + .006}}},
                                {3, 4, 5});
        for (const auto &b : {a, edge, near}) {
            // Recenter the ALREADY rounded input coordinates. Comparing against
            // the unshifted ideal fixture would confuse input quantization with
            // a loss of numerical translation invariance.
            const auto localA = recentered(a, origin), localB = recentered(b, origin);
            const auto shifted = makePairQuadrature(a, b, settings, k),
                       local = makePairQuadrature(localA, localB, settings, k);
            require(shifted.topology == local.topology, "Translation changed pair topology");
            const auto actual = operatorMoments(a, b, shifted, k),
                       expected = operatorMoments(localA, localB, local, k);
            for (const auto [start, count] :
                 std::array<std::array<unsigned, 2>, 4>{{{0, 3}, {3, 9}, {12, 3}, {15, 9}}}) {
                double differenceSquared = 0, referenceSquared = 0;
                for (unsigned index = start; index < start + count; ++index) {
                    differenceSquared += std::norm(actual[index] - expected[index]);
                    referenceSquared += std::norm(expected[index]);
                }
                if (referenceSquared == 0)
                    require(differenceSquared == 0, "Translated zero operator became nonzero");
                else
                    require(std::sqrt(differenceSquared / referenceSquared) < 2e-12,
                            "CAD translation changes S/D/K'/H operator moments");
            }
            if (shifted.topology == PairTopology::disjoint)
                require(!shifted.patches.empty() &&
                            shifted.estimatedRelativeError <= settings.relativeTolerance,
                        "Translated near-pair quadrature lost its error control");
        }
    }
}
void workspaceTests() {
    const auto a = panel({{{0, 0, 0}, {.01, 0, 0}, {0, .01, 0}}}, {0, 1, 2});
    const auto close =
        panel({{{.001, .001, .002}, {.011, .001, .002}, {.001, .011, .002}}}, {3, 4, 5});
    const auto far = panel({{{0, 0, .1}, {.01, 0, .1}, {0, .01, .1}}}, {6, 7, 8});
    PairQuadratureSettings settings;
    settings.relativeTolerance = 1e-4;
    for (const auto &b : {a, close, far}) {
        const auto original = makePairQuadrature(a, b, settings, 20);
        require(original.peakWorkspaceBytes > 0, "Quadrature omitted actual workspace accounting");
        auto bounded = settings;
        bounded.maximumWorkspaceBytes = original.peakWorkspaceBytes;
        const auto repeated = makePairQuadrature(a, b, bounded, 20);
        require(repeated.peakWorkspaceBytes == original.peakWorkspaceBytes &&
                    repeated.points.size() == original.points.size() &&
                    repeated.patches.size() == original.patches.size(),
                "Exact workspace budget changed quadrature");
        bounded.maximumWorkspaceBytes--;
        bool rejected = false;
        try {
            (void)makePairQuadrature(a, b, bounded, 20);
        } catch (const std::runtime_error &) {
            rejected = true;
        }
        require(rejected, "Quadrature exceeded its owned workspace budget");
    }
}
} // namespace
int main() {
    try {
        forcedAdaptiveFarTests();
        partitionTests();
        singularTests();
        nearTests();
        nearCornerTests();
        promotedNearTests();
        projectionQuadratureTests();
        translationTests();
        workspaceTests();
        std::cout << "Triangle-pair quadrature tests passed.\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
