#include <BRepBuilderAPI_Transform.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRepTools.hxx>
#include <TopExp.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <cmath>
#include <hornsculpt/axisymmetric_geometry.hpp>
#include <hornsculpt/axisymmetric_mesh.hpp>
#include <hornsculpt/cad_geometry.hpp>
#include <hornsculpt/draft_geometry.hpp>
#include <hornsculpt/meshing.hpp>
#include <iostream>
#include <numbers>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>

using namespace hornsculpt;
namespace {
void require(bool value, const char *message) {
    if (!value)
        throw std::runtime_error(message);
}
template <class F> void rejects(F f) {
    bool failed = false;
    try {
        f();
    } catch (const std::exception &) {
        failed = true;
    }
    require(failed, "Invalid axisymmetric input accepted.");
}
std::shared_ptr<const CadGeometry> cad(const TopoDS_Shape &shape) {
    std::ostringstream stream;
    BRepTools::Write(shape, stream, Standard_False, Standard_False, TopTools_FormatVersion_CURRENT);
    TopTools_IndexedMapOfShape faces;
    TopExp::MapShapes(shape, TopAbs_FACE, faces);
    CadBrepPart part{stream.str(),
                     std::vector<std::uint32_t>(static_cast<std::size_t>(faces.Extent()))};
    std::iota(part.faceIds.begin(), part.faceIds.end(), 101u);
    return CadGeometry::fromBrepParts(std::span{&part, 1}, 1);
}
void testCylinder() {
    const auto shape = BRepPrimAPI_MakeCylinder(.04, .03).Shape();
    const auto source = cad(shape);
    MeridianSampling settings;
    settings.maximumLength = .007;
    settings.chordTolerance = .00002;
    const auto meridian = makeAxisymmetricMeridian(*source, {}, {0, 0, 1}, settings);
    const auto coarse = buildAxisymmetricMesh(meridian, 16),
               fine = buildAxisymmetricMesh(meridian, 64);
    const double exact = std::numbers::pi * .04 * .04 * .03;
    const double a = std::abs(inspectMesh(coarse.mesh).signedVolume - exact),
                 b = std::abs(inspectMesh(fine.mesh).signedVolume - exact);
    require(b < a * .07 && b / exact < .002, "Rotational volume does not converge.");
    std::set<std::uint32_t> used;
    for (const auto &e : meridian.segments)
        used.insert(e.faceId);
    require(used.size() == 3, "CAD face identities were lost.");
    require(fine.pressureDofs.size() == fine.mesh.vertices.size() &&
                fine.fluxDofs.size() == fine.mesh.triangles.size(),
            "Rotational maps have wrong dimension.");
    require(fine.testPanels.size() * 64 == fine.mesh.triangles.size(),
            "Rotational test sector incomplete.");
    gp_Trsf placement;
    placement.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(1, 1, 0)), .7);
    placement.SetTranslationPart(gp_Vec(.13, -.07, .09));
    const auto moved = cad(BRepBuilderAPI_Transform(shape, placement, Standard_True).Shape());
    const auto origin = gp_Pnt(0, 0, 0).Transformed(placement);
    const auto axis = gp_Pnt(0, 0, 1).Transformed(placement).XYZ() - origin.XYZ();
    const auto rotated = makeAxisymmetricMeridian(*moved, {origin.X(), origin.Y(), origin.Z()},
                                                  {axis.X(), axis.Y(), axis.Z()}, settings);
    const auto q = inspectMesh(buildAxisymmetricMesh(rotated, 32).mesh);
    require(q.watertight() && q.connectedComponents == 1, "Rotated axis lost topology.");
    rejects([&] { makeAxisymmetricMeridian(*source, {.001, 0, 0}, {0, 0, 1}, settings); });
    rejects([&] { makeAxisymmetricMeridian(*source, {}, {0, 1, 0}, settings); });
    rejects([&] { buildAxisymmetricMesh(meridian, 15); });
    auto reversed = meridian;
    for (auto &e : reversed.segments)
        std::swap(e.first, e.second);
    rejects([&] { buildAxisymmetricMesh(reversed, 16); });
    std::stop_source cancellation;
    cancellation.request_stop();
    rejects([&] { buildAxisymmetricMesh(meridian, 16, cancellation.get_token()); });
}
void testSphereAndAsymmetry() {
    const auto sphere = cad(BRepPrimAPI_MakeSphere(.03).Shape());
    MeridianSampling settings;
    settings.maximumLength = .005;
    settings.chordTolerance = .00002;
    const auto m = makeAxisymmetricMeridian(*sphere, {}, {0, 0, 1}, settings);
    const auto q = inspectMesh(buildAxisymmetricMesh(m, 32).mesh);
    require(q.watertight() && q.connectedComponents == 1, "Sphere axis poles did not close.");
    const auto box = cad(BRepPrimAPI_MakeBox(.03, .03, .03).Shape());
    rejects([&] { makeAxisymmetricMeridian(*box, {}, {0, 0, 1}, settings); });
}
void testDraft() {
    DraftParameters p;
    p.verticalProfile = p.profile;
    p.verticalAngle = p.horizontalAngle;
    p.entry.height = p.entry.width;
    const DraftGeometry source(p);
    MeridianSampling settings;
    settings.maximumLength = .015;
    settings.chordTolerance = .0001;
    const auto m = makeAxisymmetricMeridian(source, settings);
    const auto mesh = buildAxisymmetricMesh(m, 16);
    require(inspectMesh(mesh.mesh).watertight(), "Circular Draft generator did not close.");
    std::set<std::uint32_t> faces;
    for (const auto &s : m.segments)
        faces.insert(s.faceId);
    require(faces.size() == 5, "Circular Draft lost physical surfaces.");
    p.verticalAngle *= .8;
    const DraftGeometry directional(p);
    rejects([&] { makeAxisymmetricMeridian(directional, settings); });
    p.verticalAngle = p.horizontalAngle;
    p.asymmetricEnabled = true;
    p.leftProfile = p.bottomProfile = p.profile;
    p.leftAngle = p.bottomAngle = p.horizontalAngle;
    // Enabling independent controls alone does not destroy actual symmetry.
    require(!makeAxisymmetricMeridian(DraftGeometry(p), settings).segments.empty(),
            "Four identical directional profiles must retain rotational eligibility.");
    p.leftAngle *= .8;
    rejects([&] { makeAxisymmetricMeridian(DraftGeometry(p), settings); });
    p.leftAngle = p.horizontalAngle;
    p.bottomProfile.terminationScale *= 1.5;
    rejects([&] { makeAxisymmetricMeridian(DraftGeometry(p), settings); });
}
} // namespace
int main() {
    try {
        testCylinder();
        testSphereAndAsymmetry();
        testDraft();
        std::cout << "Axisymmetric geometry/mesh checks passed.\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
