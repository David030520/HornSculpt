#include <hornsculpt/draft_step_export.hpp>

#include <BRepAdaptor_Surface.hxx>
#include <BRepBndLib.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepClass3d_SolidClassifier.hxx>
#include <BRep_Tool.hxx>
#include <Bnd_Box.hxx>
#include <Geom_BSplineSurface.hxx>
#include <STEPControl_Reader.hxx>
#include <TopExp.hxx>
#include <TopTools_IndexedDataMapOfShapeListOfShape.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Solid.hxx>
#include <chrono>
#include <cmath>
#include <iostream>
#include <numbers>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <vector>

using namespace hornsculpt;
namespace {
void require(bool condition, const char *message) {
    if (!condition)
        throw std::runtime_error(message);
}
TopoDS_Shape read(const DraftStepExportResult &result) {
    require(result.data.starts_with("ISO-10303-21;"), "Missing STEP header.");
    require(result.data.find("FACETED_BREP") == std::string::npos &&
                result.data.find("TRIANGULATED_FACE") == std::string::npos,
            "Draft export contains tessellation instead of analytic spline BRep.");
    std::istringstream input(result.data);
    STEPControl_Reader reader;
    reader.SetShapeProcessFlags(ShapeProcess::OperationsFlags{});
    require(reader.ReadStream("draft.step", input) == IFSelect_RetDone,
            "Exported STEP cannot be read.");
    reader.SetSystemLengthUnit(1000.);
    require(reader.TransferRoots() != 0, "Exported STEP has no transferable roots.");
    auto shape = reader.OneShape();
    require(!shape.IsNull() && BRepCheck_Analyzer(shape, true).IsValid(),
            "STEP roundtrip damaged BRep validity.");
    return shape;
}

gp_Pnt expected(const DraftGeometry &geometry, double z, double phi, double offset = 0) {
    const double r = geometry.radius(z, phi) + offset;
    return {r * std::cos(phi), r * std::sin(phi), z};
}

void surfaces(const DraftGeometry &geometry, DraftStepExportKind kind) {
    DraftStepExportOptions options;
    options.kind = kind;
    const auto exported = exportDraftStep(geometry, options);
    require(exported.maximumDeviationMetres <= options.toleranceMetres * .5,
            "Exporter accepted its own excessive interpolation discrepancy.");
    const auto shape = read(exported);
    TopTools_IndexedMapOfShape faces, solids;
    TopExp::MapShapes(shape, TopAbs_FACE, faces);
    TopExp::MapShapes(shape, TopAbs_SOLID, solids);
    const bool material = kind == DraftStepExportKind::solid;
    const auto boundaries = geometry.axialBreakpoints();
    const auto patches = static_cast<unsigned>(boundaries.size() - 1);
    require(faces.Extent() == static_cast<int>(material ? 2 * patches + 2 : patches),
            "Unexpected wall/ring surface count.");
    require(solids.Extent() == (material ? 1 : 0), "Open surface/solid mode changed.");
    unsigned splineFaces = 0;
    std::vector<Handle(Geom_Surface)> innerSurfaces(patches);
    for (int i = 1; i <= faces.Extent(); ++i) {
        const auto face = TopoDS::Face(faces(i));
        const BRepAdaptor_Surface adaptor(face);
        if (adaptor.GetType() == GeomAbs_Plane)
            continue;
        require(adaptor.GetType() == GeomAbs_BSplineSurface, "Horn wall is not a BSpline surface.");
        ++splineFaces;
        const auto surface = BRep_Tool::Surface(face);
        const double firstZ = adaptor.FirstVParameter(), lastZ = adaptor.LastVParameter();
        std::size_t segment = 0;
        while (segment < patches && std::abs(firstZ - boundaries[segment]) > 1e-9)
            ++segment;
        require(segment < patches && std::abs(lastZ - boundaries[segment + 1]) < 1e-9,
                "Spline surface crosses a physical axial segment boundary.");
        const double midZ = (firstZ + lastZ) * .5;
        const double innerDistance =
            surface->Value(.5, midZ).Distance(expected(geometry, midZ, .5));
        const double offset = innerDistance > geometry.parameters().wallThickness * .5
                                  ? geometry.parameters().wallThickness
                                  : 0;
        if (offset == 0)
            innerSurfaces[segment] = surface;
        // A different dense grid from the exporter's acceptance samples checks
        // shape fidelity after STEP serialization and SI-unit transfer.
        for (unsigned a = 0; a < 79; ++a)
            for (unsigned b = 0; b < 23; ++b) {
                const double phi = 2 * std::numbers::pi * (a + .371) / 79;
                const double z = std::lerp(firstZ, lastZ, (b + .619) / 23);
                require(surface->Value(phi, z).Distance(expected(geometry, z, phi, offset)) <=
                            options.toleranceMetres,
                        "Reimported BSpline differs from the final Draft geometry.");
            }
        if (std::abs(lastZ - geometry.length()) < 1e-9)
            for (const double fraction : {1. - 1e-4, 1. - 1e-8, 1. - 1e-12})
                require(surface->Value(.71, geometry.length() * fraction)
                                .Distance(expected(geometry, geometry.length() * fraction, .71,
                                                   offset)) <= options.toleranceMetres,
                        "Near-vertical mouth endpoint lost requested shape fidelity.");
        for (unsigned axis = 0; axis < 4; ++axis) {
            const double phi = axis * std::numbers::pi / 2;
            gp_Pnt a, b;
            gp_Vec au, av, bu, bv;
            surface->D1(phi + 1e-8, midZ, a, au, av);
            surface->D1(phi == 0 ? 2 * std::numbers::pi - 1e-8 : phi - 1e-8, midZ, b, bu, bv);
            require(a.Distance(b) < 1e-7 && au.Angle(bu) < 1e-4,
                    "Cardinal or periodic seam has a visible tangent discontinuity.");
        }
    }
    require(splineFaces == (material ? 2 * patches : patches), "Missing spline wall.");
    for (const auto &surface : innerSurfaces)
        require(!surface.IsNull(), "Missing inner surface segment.");
    TopTools_IndexedDataMapOfShapeListOfShape edgeFaces;
    TopExp::MapShapesAndAncestors(shape, TopAbs_EDGE, TopAbs_FACE, edgeFaces);
    for (std::size_t boundary = 1; boundary + 1 < boundaries.size(); ++boundary) {
        unsigned joinedRings = 0;
        for (int edge = 1; edge <= edgeFaces.Extent(); ++edge) {
            double first, last;
            const auto curve = BRep_Tool::Curve(TopoDS::Edge(edgeFaces.FindKey(edge)), first, last);
            if (curve.IsNull())
                continue;
            const auto a = curve->Value(first), b = curve->Value(last);
            if (a.Distance(b) > 1e-8 || std::abs(a.Z() - boundaries[boundary]) > 1e-8 ||
                std::abs(curve->Value((first + last) * .5).Z() - a.Z()) > 1e-8)
                continue;
            require(edgeFaces.FindFromIndex(edge).Extent() == 2,
                    "An internal slot ring is not shared by exactly two wall faces.");
            ++joinedRings;
        }
        require(joinedRings == (material ? 2u : 1u), "An internal slot ring was capped or lost.");
    }
    if (geometry.parameters().slot.enabled) {
        const auto &slot = geometry.parameters().slot;
        const double exit = slot.adapterLength + slot.straightLength;
        for (std::size_t boundary = 1; boundary + 1 < boundaries.size(); ++boundary) {
            const double z = boundaries[boundary];
            const bool sharp = std::abs(z - exit) < 1e-9 &&
                               std::isnan(geometry.radiusDerivative(z, 0));
            for (const double phi : {0., .71, 1.73, 3.2}) {
                gp_Pnt before, after;
                gp_Vec uBefore, vBefore, uAfter, vAfter;
                innerSurfaces[boundary - 1]->D1(phi, z, before, uBefore, vBefore);
                innerSurfaces[boundary]->D1(phi, z, after, uAfter, vAfter);
                require(before.Distance(after) < 1e-8, "Slot boundary is not position-continuous.");
                if (sharp) {
                    if (phi == 0)
                        require(vBefore.Angle(vAfter) > .1,
                                "Sharp slot edge was smoothed by STEP fitting.");
                } else
                    require(vBefore.Angle(vAfter) < 1e-5,
                            "Smooth slot boundary acquired a tangent jump.");
            }
        }
    }
    Bnd_Box bounds;
    BRepBndLib::AddOptimal(shape, bounds, false, false);
    double x0, y0, z0, x1, y1, z1;
    bounds.Get(x0, y0, z0, x1, y1, z1);
    std::cout << "STEP z bounds=" << z0 << ".." << z1 << " m; expected 0.." << geometry.length()
              << "; sampled error=" << exported.maximumDeviationMetres << " m\n";
    require(std::abs(z0) <= 2e-7 && std::abs(z1 - geometry.length()) <= 2e-7,
            "STEP metre units or the physical axial endpoints changed.");
    if (material) {
        const auto solid = TopoDS::Solid(solids(1));
        for (std::size_t segment = 0; segment < patches; ++segment) {
            const double z = (boundaries[segment] + boundaries[segment + 1]) * .5;
            BRepClass3d_SolidClassifier air(solid, gp_Pnt(0, 0, z), 1e-8);
            require(air.State() == TopAbs_OUT, "The horn air passage was capped or filled.");
            for (double phi : {0., 1.7, 3.4, 5.1}) {
                const auto wall = expected(geometry, z, phi, geometry.parameters().wallThickness * .5);
                BRepClass3d_SolidClassifier materialPoint(solid, wall, 1e-8);
                require(materialPoint.State() == TopAbs_IN, "Radial wall material is absent.");
            }
        }
    }
}

void curves(const DraftGeometry &geometry) {
    DraftStepExportOptions options;
    options.kind = DraftStepExportKind::sections;
    options.sectionCount = 5;
    const auto exported = exportDraftStep(geometry, options);
    const auto boundaries = geometry.axialBreakpoints();
    const auto profiles = 4 * (boundaries.size() - 1);
    require(exported.curveCount == profiles + 5 && exported.surfaceCount == 0,
            "Wrong cardinal-profile/section counts.");
    const auto shape = read(exported);
    TopTools_IndexedMapOfShape edges, faces;
    TopExp::MapShapes(shape, TopAbs_EDGE, edges);
    TopExp::MapShapes(shape, TopAbs_FACE, faces);
    require(edges.Extent() == static_cast<int>(profiles + 5) && faces.IsEmpty(),
            "Curve-only STEP contains surfaces or missing curves.");
    unsigned closed = 0, open = 0;
    for (int i = 1; i <= edges.Extent(); ++i) {
        double first, last;
        const auto curve = BRep_Tool::Curve(TopoDS::Edge(edges(i)), first, last);
        require(!curve.IsNull(), "STEP edge has no 3D curve.");
        if (curve->Value(first).Distance(curve->Value(last)) < 1e-8) {
            ++closed;
            const double z = curve->Value(first).Z();
            for (unsigned sample = 0; sample < 97; ++sample) {
                const double phi = 2 * std::numbers::pi * (sample + .41) / 97;
                require(curve->Value(phi).Distance(expected(geometry, z, phi)) <=
                            options.toleranceMetres,
                        "Closed STEP section differs from the composed contour.");
            }
        } else {
            ++open;
            const auto mouth = curve->Value(last);
            const double phi = std::atan2(mouth.Y(), mouth.X());
            std::size_t segment = 0;
            while (segment + 1 < boundaries.size() && std::abs(first - boundaries[segment]) > 1e-9)
                ++segment;
            require(segment + 1 < boundaries.size() && std::abs(last - boundaries[segment + 1]) < 1e-9,
                    "STEP longitudinal curve crosses a physical slot corner.");
            for (unsigned sample = 0; sample < 67; ++sample) {
                const double z = std::lerp(first, last, (sample + .39) / 67);
                require(curve->Value(z).Distance(expected(geometry, z, phi)) <=
                            options.toleranceMetres,
                        "STEP cardinal profile differs from the composed longitudinal curve.");
            }
        }
    }
    require(closed == 5 && open == profiles, "STEP sections lost closure or cardinal profiles.");
}

void invalidAndCancelled() {
    const DraftGeometry geometry;
    DraftStepExportOptions options;
    options.toleranceMetres = 0;
    bool rejected = false;
    try {
        (void)exportDraftStep(geometry, options);
    } catch (const std::invalid_argument &) {
        rejected = true;
    }
    require(rejected, "Invalid STEP tolerance accepted.");
    std::stop_source stop;
    stop.request_stop();
    rejected = false;
    try {
        (void)exportDraftStep(geometry, {}, stop.get_token());
    } catch (const std::runtime_error &) {
        rejected = true;
    }
    require(rejected, "Pre-cancelled STEP export started work.");
    std::stop_source activeStop;
    DraftParameters slow;
    slow.profile.terminationQ = slow.verticalProfile.terminationQ = 1;
    const DraftGeometry verticalMouth(slow);
    std::jthread cancel([&] {
        std::this_thread::sleep_for(std::chrono::milliseconds(5));
        activeStop.request_stop();
    });
    rejected = false;
    try {
        (void)exportDraftStep(verticalMouth, {}, activeStop.get_token());
    } catch (const std::runtime_error &error) {
        rejected = std::string(error.what()).find("abgebrochen") != std::string::npos;
    }
    require(rejected, "Running STEP approximation did not observe cancellation.");
}
} // namespace

int main() {
    try {
        const DraftGeometry standard;
        surfaces(standard, DraftStepExportKind::solid);
        DraftParameters thin;
        thin.wallThickness = .0001;
        surfaces(DraftGeometry(thin), DraftStepExportKind::solid);
        DraftParameters verticalMouth;
        verticalMouth.profile.terminationQ = verticalMouth.verticalProfile.terminationQ = 1;
        surfaces(DraftGeometry(verticalMouth), DraftStepExportKind::innerSurface);
        DraftParameters p;
        p.asymmetricEnabled = true;
        p.leftAngle = .43;
        p.bottomAngle = .28;
        p.leftProfile.throatAngle = .1;
        p.bottomProfile.terminationScale = .18;
        p.guideEnabled = true;
        p.guide.position = .65;
        p.guide.width = .3;
        p.guide.aspectRatio = .8;
        p.guide.rotation = .19;
        p.mouthMorph.enabled = true;
        p.mouthMorph.startFraction = .8;
        p.mouthMorph.target.width = .38;
        p.mouthMorph.target.aspectRatio = .7;
        const DraftGeometry composed(p);
        surfaces(composed, DraftStepExportKind::innerSurface);
        surfaces(composed, DraftStepExportKind::solid);
        curves(composed);
        DraftParameters slot;
        slot.slot.enabled = true;
        slot.slot.transitionLength = 0;
        slot.transition.startFraction = .3;
        slot.profile.throatAngle = .3;
        slot.verticalProfile.throatAngle = .15;
        const DraftGeometry sharpSlot(slot);
        surfaces(sharpSlot, DraftStepExportKind::innerSurface);
        surfaces(sharpSlot, DraftStepExportKind::solid);
        curves(sharpSlot);
        slot.slot.transitionLength = .008;
        const DraftGeometry smoothSlot(slot);
        surfaces(smoothSlot, DraftStepExportKind::innerSurface);
        surfaces(smoothSlot, DraftStepExportKind::solid);
        curves(smoothSlot);
        slot.slot.contour = SlotContour::roundedRectangle;
        slot.slot.cornerRadius = .002;
        slot.slot.exit = SlotExit::circularFillet;
        slot.slot.exitRadius = .002;
        const DraftGeometry roundedSlot(slot);
        surfaces(roundedSlot, DraftStepExportKind::innerSurface);
        surfaces(roundedSlot, DraftStepExportKind::solid);
        curves(roundedSlot);
        slot.slot.cornerRadius = 0;
        slot.slot.exitRadius = 0;
        const DraftGeometry sharpRectangle(slot);
        for (const auto kind : {DraftStepExportKind::innerSurface, DraftStepExportKind::solid,
                                DraftStepExportKind::sections}) {
            DraftStepExportOptions unsupported;
            unsupported.kind = kind;
            bool rejected = false;
            try {
                (void)exportDraftStep(sharpRectangle, unsupported);
            } catch (const std::invalid_argument &error) {
                rejected = std::string(error.what()).find("Eckradius 0") != std::string::npos;
            }
            require(rejected, "STEP silently smoothed a zero-radius rectangle corner.");
        }
        invalidAndCancelled();
        std::cout << "Draft STEP solid, inner surface, profile/section, SI, seam, and cancellation "
                     "checks passed.\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
