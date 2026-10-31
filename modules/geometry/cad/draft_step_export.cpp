#include <hornsculpt/draft_step_export.hpp>

#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeSolid.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepGProp.hxx>
#include <BRepLib.hxx>
#include <BRep_Builder.hxx>
#include <DESTEP_Parameters.hxx>
#include <GProp_GProps.hxx>
#include <GeomAPI_Interpolate.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Precision.hxx>
#include <STEPControl_Writer.hxx>
#include <Standard_Failure.hxx>
#include <StepData_StepModel.hxx>
#include <TColStd_Array1OfInteger.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <TColStd_HArray1OfBoolean.hxx>
#include <TColStd_HArray1OfReal.hxx>
#include <TColgp_Array1OfVec.hxx>
#include <TColgp_Array2OfPnt.hxx>
#include <TColgp_HArray1OfPnt.hxx>
#include <TopExp.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Shell.hxx>
#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <gp_Pln.hxx>
#include <limits>
#include <numbers>
#include <numeric>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace hornsculpt {
namespace {
constexpr double fullTurn = 2 * std::numbers::pi;
constexpr std::size_t maximumAngularSamples = 512;
constexpr std::size_t maximumAxialSamples = 1025;
constexpr std::size_t maximumSurfaceSamples = 65536;
constexpr unsigned maximumRefinementPasses = 48;
// Independent positions between interpolation nodes, including both cell ends.
constexpr std::array checkFractions{0., .173, .419, .733, 1.};

void checkStop(std::stop_token stop) {
    if (stop.stop_requested())
        throw std::runtime_error("STEP-Export abgebrochen.");
}

gp_Pnt hornPoint(const DraftGeometry &geometry, double z, double phi, double thickness = 0) {
    // Evaluate the seam using one identical point, avoiding two rounded sin(2pi) endpoints.
    if (phi == fullTurn)
        phi = 0;
    const double r = geometry.radius(z, phi) + thickness;
    return {r * std::cos(phi), r * std::sin(phi), z};
}

std::vector<double> axialNodes(const DraftGeometry &geometry) {
    std::vector<double> nodes;
    for (unsigned i = 0; i <= 16; ++i)
        nodes.push_back(geometry.length() * i / 16);
    const auto &p = geometry.parameters();
    const auto fraction = [&](double value) { nodes.push_back(geometry.length() * value); };
    if (p.transition.method == TransitionMethod::lut)
        for (const auto &knot : p.transition.knots)
            fraction(knot.fraction);
    else {
        fraction(p.transition.startFraction);
        fraction(p.transition.endFraction);
    }
    if (p.guideEnabled)
        nodes.push_back(geometry.guidePosition());
    if (p.mouthMorph.enabled) {
        fraction(p.mouthMorph.startFraction);
        fraction(p.mouthMorph.endFraction);
    }
    std::sort(nodes.begin(), nodes.end());
    nodes.erase(std::unique(nodes.begin(), nodes.end()), nodes.end());
    return nodes;
}

std::vector<std::vector<double>> axialSegments(const DraftGeometry &geometry) {
    const auto nodes = axialNodes(geometry);
    const auto boundaries = geometry.axialBreakpoints();
    std::vector<std::vector<double>> segments;
    for (std::size_t i = 0; i + 1 < boundaries.size(); ++i) {
        std::vector<double> segment;
        const double first = boundaries[i], last = boundaries[i + 1];
        for (const double z : nodes)
            if (z >= first && z <= last)
                segment.push_back(z);
        // Short adapters and straight slots need their own interpolation grid;
        // a global grid can miss the entire physical segment.
        for (unsigned part = 0; part <= 4; ++part)
            segment.push_back(std::lerp(first, last, part / 4.));
        std::sort(segment.begin(), segment.end());
        segment.erase(std::unique(segment.begin(), segment.end()), segment.end());
        segments.push_back(std::move(segment));
    }
    return segments;
}

std::vector<double> angularNodes() {
    std::vector<double> result;
    for (unsigned i = 0; i <= 32; ++i)
        result.push_back(fullTurn * i / 32);
    return result;
}

void refine(std::vector<double> &nodes, const std::vector<bool> &intervals) {
    std::vector<double> refined;
    for (std::size_t i = 0; i + 1 < nodes.size(); ++i) {
        refined.push_back(nodes[i]);
        if (intervals[i]) {
            const double middle = std::midpoint(nodes[i], nodes[i + 1]);
            if (middle == nodes[i] || middle == nodes[i + 1])
                throw std::runtime_error("STEP-Approximation erreicht die Zahlengenauigkeit.");
            refined.push_back(middle);
        }
    }
    refined.push_back(nodes.back());
    nodes = std::move(refined);
}

using CurveFunction = std::function<gp_Pnt(double)>;
using EndpointTangents = std::array<std::optional<gp_Vec>, 2>;

Handle(Geom_BSplineCurve) interpolate(const std::vector<double> &nodes, const CurveFunction &point,
                                      bool periodic, std::stop_token stop,
                                      const EndpointTangents &endTangents = {}) {
    const auto count = static_cast<int>(nodes.size()) - (periodic ? 1 : 0);
    Handle(TColgp_HArray1OfPnt) points = new TColgp_HArray1OfPnt(1, count);
    Handle(TColStd_HArray1OfReal) parameters =
        new TColStd_HArray1OfReal(1, static_cast<int>(nodes.size()));
    for (int i = 1; i <= count; ++i) {
        checkStop(stop);
        points->SetValue(i, point(nodes[static_cast<std::size_t>(i - 1)]));
    }
    for (int i = 1; i <= parameters->Length(); ++i)
        parameters->SetValue(i, nodes[static_cast<std::size_t>(i - 1)]);
    // This is an interpolation coincidence threshold, not the export accuracy.
    GeomAPI_Interpolate fit(points, parameters, periodic, 1e-12);
    if (periodic) {
        TColgp_Array1OfVec tangents(1, count);
        Handle(TColStd_HArray1OfBoolean) flags = new TColStd_HArray1OfBoolean(1, count);
        flags->Init(false);
        // Opposite-side superellipse profiles are C1 at the cardinal axes.
        // Keep explicit tangency knots there instead of imposing C2 across
        // the source curvature jumps. A centered angular derivative also
        // includes rotated guides and mouth transformations.
        constexpr double derivativeStep = 1e-5; // radians, centered O(h^2)
        for (int i = 1; i <= count; ++i) {
            const double angle = nodes[static_cast<std::size_t>(i - 1)];
            const double quadrant = angle / (std::numbers::pi / 2);
            if (std::abs(quadrant - std::round(quadrant)) > 1e-13)
                continue;
            tangents.SetValue(i,
                              gp_Vec(point(angle - derivativeStep), point(angle + derivativeStep)) /
                                  (2 * derivativeStep));
            flags->SetValue(i, true);
        }
        fit.Load(tangents, flags, false);
    } else if (endTangents[0] || endTangents[1]) {
        TColgp_Array1OfVec tangents(1, count);
        Handle(TColStd_HArray1OfBoolean) flags = new TColStd_HArray1OfBoolean(1, count);
        flags->Init(false);
        for (std::size_t endpoint = 0; endpoint < 2; ++endpoint)
            if (endTangents[endpoint]) {
                const int index = endpoint == 0 ? 1 : count;
                tangents.SetValue(index, *endTangents[endpoint]);
                flags->SetValue(index, true);
            }
        fit.Load(tangents, flags, false);
    }
    checkStop(stop);
    fit.Perform();
    checkStop(stop);
    if (!fit.IsDone())
        throw std::runtime_error("STEP-BSpline-Interpolation ist fehlgeschlagen.");
    return fit.Curve();
}

bool sameBasis(const Handle(Geom_BSplineCurve) & a, const Handle(Geom_BSplineCurve) & b) {
    if (a->Degree() != b->Degree() || a->NbPoles() != b->NbPoles() ||
        a->NbKnots() != b->NbKnots() || a->IsPeriodic() != b->IsPeriodic() || a->IsRational() ||
        b->IsRational())
        return false;
    for (int i = 1; i <= a->NbKnots(); ++i)
        if (a->Knot(i) != b->Knot(i) || a->Multiplicity(i) != b->Multiplicity(i))
            return false;
    return true;
}

gp_Vec axialTangent(const DraftGeometry &geometry, double z, double phi, AxialDerivativeSide side) {
    const double slope = geometry.radiusDerivative(z, phi, side);
    if (!std::isfinite(slope))
        throw std::runtime_error("STEP-Abschnittsgrenze besitzt keine endliche Tangente.");
    return {slope * std::cos(phi), slope * std::sin(phi), 1};
}

std::vector<gp_Vec> axialTangentPoles(const DraftGeometry &geometry,
                                    const std::vector<double> &angles,
                                    const Handle(Geom_BSplineCurve) &ring, double z,
                                    double thickness, AxialDerivativeSide side,
                                    std::stop_token stop) {
    if (z == 0 || z == geometry.length())
        return {};
    // Angular interpolation is linear in positions and its prescribed angular
    // tangents. Interpolate a shifted ring to obtain derivative control poles
    // without trying to interpolate a degenerate zero-radius derivative ring.
    constexpr double derivativeScaleMetres = 1e-4;
    const auto shifted = interpolate(
        angles,
        [&](double phi) {
            return hornPoint(geometry, z, phi, thickness)
                .Translated(axialTangent(geometry, z, phi, side) * derivativeScaleMetres);
        },
        true, stop);
    if (!sameBasis(ring, shifted))
        throw std::runtime_error("STEP-Randtangenten besitzen keine gemeinsame Splinebasis.");
    std::vector<gp_Vec> result;
    for (int pole = 1; pole <= ring->NbPoles(); ++pole)
        result.push_back(gp_Vec(ring->Pole(pole), shifted->Pole(pole)) / derivativeScaleMetres);
    return result;
}

Handle(Geom_BSplineSurface)
    interpolateSurface(const DraftGeometry &geometry, const std::vector<double> &angles,
                       const std::vector<double> &axial, double thickness, std::stop_token stop) {
    std::vector<Handle(Geom_BSplineCurve)> rings;
    for (const double z : axial) {
        rings.push_back(interpolate(
            angles, [&](double phi) { return hornPoint(geometry, z, phi, thickness); }, true,
            stop));
        if (!sameBasis(rings.front(), rings.back()))
            throw std::runtime_error("STEP-Querschnitte besitzen keine gemeinsame Splinebasis.");
    }
    std::vector<Handle(Geom_BSplineCurve)> columns;
    const auto firstTangents = axialTangentPoles(geometry, angles, rings.front(), axial.front(),
                                                thickness, AxialDerivativeSide::after, stop);
    const auto lastTangents = axialTangentPoles(geometry, angles, rings.back(), axial.back(),
                                               thickness, AxialDerivativeSide::before, stop);
    const int uPoles = rings.front()->NbPoles();
    for (int i = 1; i <= uPoles; ++i) {
        checkStop(stop);
        const auto point = [&](double z) {
            const auto row = static_cast<std::size_t>(
                std::lower_bound(axial.begin(), axial.end(), z) - axial.begin());
            return rings[row]->Pole(i);
        };
        EndpointTangents tangents;
        if (!firstTangents.empty())
            tangents[0] = firstTangents[static_cast<std::size_t>(i - 1)];
        if (!lastTangents.empty())
            tangents[1] = lastTangents[static_cast<std::size_t>(i - 1)];
        columns.push_back(interpolate(axial, point, false, stop, tangents));
        if (!sameBasis(columns.front(), columns.back()))
            throw std::runtime_error("STEP-Laengsprofile besitzen keine gemeinsame Splinebasis.");
    }
    const auto &u = rings.front();
    const auto &v = columns.front();
    TColgp_Array2OfPnt poles(1, uPoles, 1, v->NbPoles());
    for (int i = 1; i <= uPoles; ++i)
        for (int j = 1; j <= v->NbPoles(); ++j)
            poles.SetValue(i, j, columns[static_cast<std::size_t>(i - 1)]->Pole(j));
    TColStd_Array1OfReal uKnots(1, u->NbKnots()), vKnots(1, v->NbKnots());
    TColStd_Array1OfInteger uMults(1, u->NbKnots()), vMults(1, v->NbKnots());
    u->Knots(uKnots);
    v->Knots(vKnots);
    u->Multiplicities(uMults);
    v->Multiplicities(vMults);
    return new Geom_BSplineSurface(poles, uKnots, vKnots, uMults, vMults, u->Degree(), v->Degree(),
                                   u->IsPeriodic(), false);
}

struct SurfacePatch {
    Handle(Geom_BSplineSurface) inner, outer;
};

struct SurfaceFit {
    std::vector<SurfacePatch> patches;
    double maximumDeviation{};
};

SurfaceFit fitSurfaces(const DraftGeometry &geometry, const DraftStepExportOptions &options,
                       std::stop_token stop) {
    auto angles = angularNodes();
    auto segments = axialSegments(geometry);
    const bool outer = options.kind == DraftStepExportKind::solid;
    for (unsigned pass = 0; pass < maximumRefinementPasses; ++pass) {
        checkStop(stop);
        std::size_t axialSamples = 0;
        for (const auto &segment : segments)
            axialSamples += segment.size();
        if (angles.size() > maximumAngularSamples + 1 || axialSamples > maximumAxialSamples ||
            angles.size() * axialSamples > maximumSurfaceSamples)
            throw std::runtime_error("STEP-Approximation ueberschreitet das Kontrollpunktbudget.");
        SurfaceFit fit;
        std::vector<bool> angularRefine(angles.size() - 1);
        bool changed = false;
        // Every patch uses the same angular basis, so adjacent fitted rings
        // coincide exactly. Sewing must not conceal the approximation error
        // of independently fitted rings by enlarging the BRep tolerance.
        for (auto &axial : segments) {
            SurfacePatch patch;
            patch.inner = interpolateSurface(geometry, angles, axial, 0, stop);
            if (outer)
                patch.outer = interpolateSurface(geometry, angles, axial,
                                                 geometry.parameters().wallThickness, stop);
            std::vector<bool> axialRefine(axial.size() - 1);
            for (const auto &surface : {patch.inner, patch.outer}) {
                if (surface.IsNull())
                    continue;
                const double offset =
                    surface == patch.outer ? geometry.parameters().wallThickness : 0;
                for (std::size_t i = 0; i + 1 < angles.size(); ++i)
                    for (std::size_t j = 0; j + 1 < axial.size(); ++j) {
                        checkStop(stop);
                        double angularError = 0, axialError = 0, mixedError = 0;
                        for (const double u : checkFractions)
                            for (const double v : checkFractions) {
                                const double phi = std::lerp(angles[i], angles[i + 1], u);
                                const double z = std::lerp(axial[j], axial[j + 1], v);
                                gp_Pnt actual;
                                gp_Vec du, dv;
                                surface->D1(phi, z, actual, du, dv);
                                const auto expected = hornPoint(geometry, z, phi, offset);
                                double error = actual.Distance(expected);
                                const gp_Vec radial(std::cos(phi), std::sin(phi), 0);
                                const double orientation = du.Crossed(dv).Dot(radial);
                                if (!std::isfinite(error) || !std::isfinite(orientation) ||
                                    orientation <= 0)
                                    error = std::numeric_limits<double>::infinity();
                                fit.maximumDeviation = std::max(fit.maximumDeviation, error);
                                if (v == 0 || v == 1)
                                    angularError = std::max(angularError, error);
                                if (u == 0 || u == 1)
                                    axialError = std::max(axialError, error);
                                mixedError = std::max(mixedError, error);
                            }
                        // Reserve half of the requested accuracy for unsampled
                        // positions and downstream STEP roundoff. This remains a
                        // sampled check, not a claimed global approximation theorem.
                        const double limit = options.toleranceMetres * .5;
                        bool angular = angularError > limit, axialSplit = axialError > limit;
                        if (mixedError > limit && !angular && !axialSplit)
                            angular = axialSplit = true;
                        angularRefine[i] = angularRefine[i] || angular;
                        axialRefine[j] = axialRefine[j] || axialSplit;
                        changed = changed || angular || axialSplit;
                    }
            }
            refine(axial, axialRefine);
            fit.patches.push_back(std::move(patch));
        }
        if (!changed)
            return fit;
        refine(angles, angularRefine);
    }
    throw std::runtime_error("STEP-Approximation erreicht die angeforderte Genauigkeit nicht.");
}

TopoDS_Wire wire(const Handle(Geom_Curve) & curve) {
    BRepBuilderAPI_MakeEdge edge(curve);
    if (!edge.IsDone())
        throw std::runtime_error("STEP-Kurve kann nicht als BRep-Kante dargestellt werden.");
    BRepBuilderAPI_MakeWire result(edge.Edge());
    if (!result.IsDone())
        throw std::runtime_error("STEP-Kontur kann nicht als Draht dargestellt werden.");
    return result.Wire();
}

void requireValid(const TopoDS_Shape &shape) {
    if (shape.IsNull() || !BRepCheck_Analyzer(shape, true).IsValid())
        throw std::runtime_error("STEP-Export erzeugt keine gueltige BRep-Topologie.");
}

TopoDS_Shape surfaceShape(const DraftGeometry &geometry, const SurfaceFit &fit, bool solid,
                          std::stop_token stop) {
    const double topologyTolerance = Precision::Confusion();
    BRepBuilderAPI_Sewing sewing(topologyTolerance);
    for (const auto &patch : fit.patches) {
        auto inner = BRepBuilderAPI_MakeFace(patch.inner, topologyTolerance).Face();
        inner.Reverse();
        if (!solid && fit.patches.size() == 1) {
            TopoDS_Shell shell;
            BRep_Builder builder;
            builder.MakeShell(shell);
            builder.Add(shell, inner);
            requireValid(shell);
            return shell;
        }
        sewing.Add(inner);
        if (solid)
            sewing.Add(BRepBuilderAPI_MakeFace(patch.outer, topologyTolerance).Face());
    }
    // Only the physical entry and mouth have annular material closures.
    // Internal slot boundaries remain open, shared edges between patches.
    if (solid)
        for (const double z : {0., geometry.length()}) {
            const auto &patch = z == 0 ? fit.patches.front() : fit.patches.back();
            const auto outer = wire(patch.outer->VIso(z));
            auto hole = wire(patch.inner->VIso(z));
            hole.Reverse();
            BRepBuilderAPI_MakeFace cap(gp_Pln(gp_Pnt(0, 0, z), gp_Dir(0, 0, 1)), outer, true);
            cap.Add(hole);
            if (!cap.IsDone())
                throw std::runtime_error("STEP-Ringabschluss kann nicht gebildet werden.");
            auto face = cap.Face();
            if (z == 0)
                face.Reverse();
            sewing.Add(face);
        }
    checkStop(stop);
    sewing.Perform();
    checkStop(stop);
    const auto sewed = sewing.SewedShape();
    if (sewed.ShapeType() != TopAbs_SHELL || (solid && sewing.NbFreeEdges() != 0) ||
        sewing.NbMultipleEdges() != 0)
        throw std::runtime_error(
            "STEP-Hornwand konnte nicht zu einer zusammenhaengenden Schale verbunden werden.");
    if (!solid) {
        requireValid(sewed);
        return sewed;
    }
    BRepBuilderAPI_MakeSolid makeSolid(TopoDS::Shell(sewed));
    if (!makeSolid.IsDone())
        throw std::runtime_error("STEP-Hornwand konnte nicht als Volumenkoerper gebildet werden.");
    auto result = makeSolid.Solid();
    if (!BRepLib::OrientClosedSolid(result))
        throw std::runtime_error("STEP-Hornwand besitzt keine eindeutige Volumenorientierung.");
    requireValid(result);
    GProp_GProps properties;
    BRepGProp::VolumeProperties(result, properties);
    if (!std::isfinite(properties.Mass()) || properties.Mass() <= 0)
        throw std::runtime_error("STEP-Hornwand besitzt kein positives Materialvolumen.");
    return result;
}

Handle(Geom_BSplineCurve)
    fitCurve(std::vector<double> nodes, const CurveFunction &point, bool periodic, double tolerance,
             double &maximum, std::stop_token stop, const EndpointTangents &tangents = {}) {
    for (unsigned pass = 0; pass < maximumRefinementPasses; ++pass) {
        checkStop(stop);
        if (nodes.size() > (periodic ? maximumAngularSamples + 1 : maximumAxialSamples))
            throw std::runtime_error("STEP-Kurve ueberschreitet das Kontrollpunktbudget.");
        const auto curve = interpolate(nodes, point, periodic, stop, tangents);
        std::vector<bool> refineIntervals(nodes.size() - 1);
        double measured = 0;
        bool changed = false;
        for (std::size_t i = 0; i + 1 < nodes.size(); ++i) {
            checkStop(stop);
            for (const double fraction : checkFractions) {
                const double parameter = std::lerp(nodes[i], nodes[i + 1], fraction);
                gp_Pnt fitted;
                gp_Vec tangent;
                curve->D1(parameter, fitted, tangent);
                double error = fitted.Distance(point(parameter));
                if (!std::isfinite(error))
                    throw std::runtime_error("STEP-Kurve erzeugt nicht endliche Punkte.");
                const double orientation =
                    periodic ? fitted.X() * tangent.Y() - fitted.Y() * tangent.X() : tangent.Z();
                if (!std::isfinite(orientation) || orientation <= 0)
                    error = std::numeric_limits<double>::infinity();
                measured = std::max(measured, error);
                refineIntervals[i] = refineIntervals[i] || error > tolerance * .5;
            }
            changed = changed || refineIntervals[i];
        }
        if (!changed) {
            maximum = std::max(maximum, measured);
            return curve;
        }
        refine(nodes, refineIntervals);
    }
    throw std::runtime_error("STEP-Kurve erreicht die angeforderte Genauigkeit nicht.");
}

TopoDS_Shape contours(const DraftGeometry &geometry, const DraftStepExportOptions &options,
                      DraftStepExportResult &report, std::stop_token stop) {
    TopoDS_Compound result;
    BRep_Builder builder;
    builder.MakeCompound(result);
    const auto add = [&](const Handle(Geom_BSplineCurve) & curve) {
        const auto contour = wire(curve);
        requireValid(contour);
        builder.Add(result, contour);
        ++report.curveCount;
    };
    const auto segments = axialSegments(geometry);
    for (unsigned side = 0; side < 4; ++side) {
        const double phi = side * std::numbers::pi / 2;
        for (const auto &segment : segments) {
            EndpointTangents tangents;
            if (segment.front() != 0)
                tangents[0] = axialTangent(geometry, segment.front(), phi, AxialDerivativeSide::after);
            if (segment.back() != geometry.length())
                tangents[1] = axialTangent(geometry, segment.back(), phi, AxialDerivativeSide::before);
            add(fitCurve(
                segment, [&](double z) { return hornPoint(geometry, z, phi); }, false,
                options.toleranceMetres, report.maximumDeviationMetres, stop, tangents));
        }
    }
    for (std::size_t section = 0; section < options.sectionCount; ++section) {
        const double z = geometry.length() * section / (options.sectionCount - 1);
        add(fitCurve(
            angularNodes(), [&](double phi) { return hornPoint(geometry, z, phi); }, true,
            options.toleranceMetres, report.maximumDeviationMetres, stop));
    }
    requireValid(result);
    return result;
}
} // namespace

DraftStepExportResult exportDraftStep(const DraftGeometry &geometry,
                                      const DraftStepExportOptions &options, std::stop_token stop) {
    checkStop(stop);
    if (!std::isfinite(options.toleranceMetres) ||
        options.toleranceMetres < 10 * Precision::Confusion() ||
        (options.kind != DraftStepExportKind::solid &&
         options.kind != DraftStepExportKind::innerSurface &&
         options.kind != DraftStepExportKind::sections) ||
        (options.kind == DraftStepExportKind::sections &&
         (options.sectionCount < 2 || options.sectionCount > 101)))
        throw std::invalid_argument(
            "STEP: Genauigkeit mindestens 0,001 mm und 2..101 Querschnitte erforderlich.");
    if (options.kind == DraftStepExportKind::solid &&
        geometry.parameters().wallThickness <= 10 * Precision::Confusion())
        throw std::invalid_argument("STEP: Die Wand ist fuer die BRep-Topologietoleranz zu duenn.");
    const auto &slot = geometry.parameters().slot;
    // A periodic smooth surface cannot preserve moving sharp rectangle corners.
    // Reject this case instead of silently rounding the user's zero-radius edge.
    if (slot.enabled && slot.contour == SlotContour::roundedRectangle && slot.cornerRadius == 0)
        throw std::invalid_argument(
            "STEP: Scharfe Rechteck-Spalte mit Eckradius 0 benoetigen mitwandernde "
            "Kantennähte. Dieser Export unterstuetzt positive Eckradien; Radius im "
            "Entwurf ausdrücklich einstellen. Vorschau und Netz bleiben verfügbar.");
    try {
        DraftStepExportResult result;
        TopoDS_Shape shape;
        if (options.kind == DraftStepExportKind::sections)
            shape = contours(geometry, options, result, stop);
        else {
            const auto fit = fitSurfaces(geometry, options, stop);
            result.maximumDeviationMetres = fit.maximumDeviation;
            shape = surfaceShape(geometry, fit, options.kind == DraftStepExportKind::solid, stop);
            TopTools_IndexedMapOfShape faces;
            TopExp::MapShapes(shape, TopAbs_FACE, faces);
            result.surfaceCount = static_cast<std::size_t>(faces.Extent());
        }
        checkStop(stop);
        STEPControl_Writer writer;
        writer.SetShapeProcessFlags(ShapeProcess::OperationsFlags{});
        writer.Model()->SetLocalLengthUnit(1000.0);
        writer.Model()->SetWriteLengthUnit(1000.0);
        writer.SetTolerance(Precision::Confusion());
        DESTEP_Parameters parameters;
        parameters.WriteUnit = UnitsMethods_LengthUnit_Meter;
        parameters.WriteTessellated = DESTEP_Parameters::RWMode_Tessellated_Off;
        parameters.WriteProductName = "HornSculpt Draft";
        const auto mode = options.kind == DraftStepExportKind::sections
                              ? STEPControl_GeometricCurveSet
                              : STEPControl_AsIs;
        if (writer.Transfer(shape, mode, parameters) != IFSelect_RetDone)
            throw std::runtime_error("STEP-Modelluebertragung ist fehlgeschlagen.");
        checkStop(stop);
        std::ostringstream output;
        if (writer.WriteStream(output) != IFSelect_RetDone || !output)
            throw std::runtime_error("STEP-Ausgabe ist fehlgeschlagen.");
        checkStop(stop);
        result.data = output.str();
        return result;
    } catch (const Standard_Failure &error) {
        throw std::runtime_error(std::string("STEP-Export: ") + error.GetMessageString());
    }
}
} // namespace hornsculpt
