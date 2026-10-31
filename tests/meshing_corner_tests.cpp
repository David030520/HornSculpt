#include "../modules/meshing/src/surface_charts.hpp"
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepTools.hxx>
#include <TopExp.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <hornsculpt/cad_geometry.hpp>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <numbers>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>

using namespace hornsculpt;
namespace {
void require(bool condition, const char *message) {
    if (!condition)
        throw std::runtime_error(message);
}
template <class Function> std::string rejection(Function operation) {
    try {
        operation();
    } catch (const std::exception &error) {
        return error.what();
    }
    throw std::runtime_error("An invalid retained planar corner was accepted.");
}
std::shared_ptr<const CadGeometry> prism(double angleDegrees) {
    const double angle = angleDegrees * std::numbers::pi / 180.;
    BRepBuilderAPI_MakePolygon polygon;
    polygon.Add(gp_Pnt(0, 0, 0));
    polygon.Add(gp_Pnt(.02, 0, 0));
    polygon.Add(gp_Pnt(.02 * std::cos(angle), .02 * std::sin(angle), 0));
    polygon.Close();
    const auto shape =
        BRepPrimAPI_MakePrism(BRepBuilderAPI_MakeFace(polygon.Wire()).Face(), gp_Vec(0, 0, .01))
            .Shape();
    std::ostringstream data;
    BRepTools::Write(shape, data, Standard_False, Standard_False, TopTools_FormatVersion_CURRENT);
    TopTools_IndexedMapOfShape faces;
    TopExp::MapShapes(shape, TopAbs_FACE, faces);
    std::vector<std::uint32_t> ids(static_cast<std::size_t>(faces.Extent()));
    std::iota(ids.begin(), ids.end(), 1);
    const CadBrepPart part{data.str(), ids};
    return CadGeometry::fromBrepParts(std::span{&part, 1}, 1);
}
// The optional marker proves that a non-blocking case reaches front geometry,
// without asking an incompatible cusp to produce a successful quality mesh.
class CornerProbe final : public SurfaceGeometry {
  public:
    explicit CornerProbe(std::shared_ptr<const CadGeometry> source, bool stopAtFront = false)
        : records(source->boundaryCorners().begin(), source->boundaryCorners().end()),
          source_(std::move(source)), stopAtFront_(stopAtFront) {}
    std::span<const BoundaryCorner> boundaryCorners() const noexcept override {
        ++queries;
        return records;
    }
    SurfaceProjection projectSurface(std::uint32_t face, Vec3 p) const override {
        checkFront();
        return source_->projectSurface(face, p);
    }
    SurfaceProjection projectSupportingSurface(std::uint32_t face, Vec3 p) const override {
        checkFront();
        return source_->projectSupportingSurface(face, p);
    }
    Vec3 projectBoundary(std::span<const std::uint32_t> faces, Vec3 p) const override {
        checkFront();
        return source_->projectBoundary(faces, p);
    }
    Vec3 projectJunction(std::span<const std::uint32_t> faces, Vec3 p) const override {
        checkFront();
        return source_->projectJunction(faces, p);
    }
    mutable std::atomic<unsigned> queries{};
    std::vector<BoundaryCorner> records;

  private:
    void checkFront() const {
        if (stopAtFront_ && queries.load())
            throw std::runtime_error("Test front geometry reached.");
    }
    std::shared_ptr<const CadGeometry> source_;
    bool stopAtFront_;
};
MeshingOptions options(std::shared_ptr<const SurfaceGeometry> geometry) {
    MeshingOptions result;
    result.targetEdgeLength = .008;
    result.adaptive.enabled = false;
    result.geometry = std::move(geometry);
    return result;
}
void publicEntryLimits() {
    const auto source = prism(12);
    const auto preview = source->previewMesh(.00025);
    require(inspectMesh(preview).watertight(), "The planar wedge source is not closed.");
    auto probe = std::make_shared<CornerProbe>(source);
    auto limits = options(probe);
    const auto angleError = rejection([&] { (void)remeshIsotropic(preview, limits); });
    require(angleError.find("Flaeche ") != std::string::npos &&
                angleError.find("Mindestwinkel 20") != std::string::npos && probe->queries == 1,
            "Public meshing did not reject the retained planar line corner before the front.");
    const auto atlas = detail::makeSurfaceAtlas(preview, limits);
    std::size_t caps = 0;
    for (const auto &record : source->boundaryCorners())
        if (record.constraint && record.constraint->loopNormal.z > .9 &&
            std::hypot(record.position.x, record.position.y) < 1e-12) {
            probe->records = {record};
            require(rejection([&] {
                        detail::validateAtlasBoundaryCorners(atlas, limits);
                    }).find("Mindestwinkel 20") != std::string::npos,
                    "The reversed solid cap lost its retained native line proof.");
            ++caps;
        }
    require(caps == 2, "Both oppositely oriented prism caps were not tested.");
    probe->records.assign(source->boundaryCorners().begin(), source->boundaryCorners().end());
    limits.quality.minimumAngleDegrees = 1;
    const auto aspectError = rejection([&] { (void)remeshIsotropic(preview, limits); });
    require(aspectError.find("Seitenverhaeltnis 4") != std::string::npos,
            "The necessary aspect constraint was ignored at a retained planar corner.");
    limits.quality.maximumAspectRatio = 10;
    const auto accepted = remeshIsotropic(preview, limits);
    require(accepted.quality.watertight() && accepted.quality.minimumAngleDegrees >= 1 &&
                accepted.quality.maximumAspectRatio <= 10,
            "Generous requested quality limits were artificially blocked.");
    limits.quality.enabled = false;
    const auto count = probe->queries.load();
    const auto unrestricted = remeshIsotropic(preview, limits);
    require(unrestricted.quality.watertight() && probe->queries == count,
            "Disabled shape quality still consumed the early corner gate.");
}
void retainedEvidenceOnly() {
    const auto source = prism(30);
    auto limits = options(source);
    limits.quality.minimumAngleDegrees = 40;
    auto atlas = detail::makeSurfaceAtlas(source->previewMesh(.00025), limits);
    auto probe = std::make_shared<CornerProbe>(source);
    limits.geometry = probe;
    const auto found =
        std::find_if(probe->records.begin(), probe->records.end(), [](const auto &c) {
            return c.constraint && c.constraint->loopNormal.z > .9 && c.position.z > .009 &&
                   std::hypot(c.position.x, c.position.y) < 1e-12;
        });
    require(found != probe->records.end(), "Planar constraint evidence is absent.");
    const auto record = *found;
    // Synthetic provider records exercise trust/retention checks, not CAD math.
    probe->records = {record};
    auto &corner = probe->records.front();
    require(rejection([&] {
                detail::validateAtlasBoundaryCorners(atlas, limits);
            }).find("Mindestwinkel 40") != std::string::npos,
            "The retention regression does not start with a genuinely rejected corner.");
    const auto outgoing = corner.constraint->boundaryEnds[1];
    corner.constraint->boundaryEnds[1] = corner.constraint->boundaryEnds[0];
    detail::validateAtlasBoundaryCorners(atlas, limits);
    corner.constraint->boundaryEnds[1] = outgoing;
    corner.position.x += .001;
    detail::validateAtlasBoundaryCorners(atlas, limits);
    corner = record;
    corner.faceId = 9999;
    corner.position.x = std::numeric_limits<double>::quiet_NaN();
    detail::validateAtlasBoundaryCorners(atlas, limits);
    corner = record;
    corner.position.x = std::numeric_limits<double>::infinity();
    require(rejection([&] {
                detail::validateAtlasBoundaryCorners(atlas, limits);
            }).find("nicht endlichen Eckpunkt") != std::string::npos,
            "Corrupt active corner data were ignored.");
    corner = record;
    corner.constraint->loopNormal.z = std::numeric_limits<double>::quiet_NaN();
    require(rejection([&] {
                detail::validateAtlasBoundaryCorners(atlas, limits);
            }).find("Ecknormalenrichtung") != std::string::npos,
            "Corrupt active constraint data were ignored.");
    corner = record;
    corner.constraint.reset();
    detail::validateAtlasBoundaryCorners(atlas, limits);
    corner = record;
    corner.constraint->loopNormal = {-record.constraint->loopNormal.x,
                                     -record.constraint->loopNormal.y,
                                     -record.constraint->loopNormal.z};
    detail::validateAtlasBoundaryCorners(atlas, limits); // reflex complement is not acute
}
void curvedCuspIsNotAnImpossibilityProof() {
    std::ifstream stream(HS_CUSP_CORNER_FIXTURE_PATH, std::ios::binary);
    require(static_cast<bool>(stream), "The real curved cusp fixture is missing.");
    std::ostringstream bytes;
    bytes << stream.rdbuf();
    const CadBrepPart part{bytes.str(), {1, 2, 3, 4, 5, 6}};
    const auto source = CadGeometry::fromBrepParts(std::span{&part, 1}, 1);
    const Vec3 cusp{-.117357879186514, .41604817913815501, -.24851544106888598};
    for (const auto &corner : source->boundaryCorners())
        if (std::hypot(corner.position.x - cusp.x, corner.position.y - cusp.y,
                       corner.position.z - cusp.z) < 1e-10) {
            require(!corner.constraint, "A curved tangent cusp became a straight planar proof.");
        }
    auto probe = std::make_shared<CornerProbe>(source, true);
    const auto error =
        rejection([&] { (void)remeshIsotropic(source->previewMesh(.00025), options(probe)); });
    require(error == "Test front geometry reached." && probe->queries == 1,
            "A curved cusp was blocked as impossible instead of reaching the front.");

    const auto native = source->meshingSource(.00025);
    auto limits = options(source);
    limits.targetEdgeLength = .04;
    limits.parameterization = native.parameters;
    const auto refinementError = rejection([&] { (void)remeshIsotropic(native.mesh, limits); });
    require(
        (refinementError.find("shared-boundary refinement budget on face 1") != std::string::npos ||
         refinementError.find(
             "shared-boundary refinement stopped without quality progress on face 1") !=
             std::string::npos) &&
            refinementError.find("min angle ") != std::string::npos &&
            refinementError.find("max aspect ") != std::string::npos &&
            refinementError.find("Worst rejected cell:") != std::string::npos,
        "Stopping shared refinement lost the offending chart's measured quality and cell.");
}
} // namespace
int main() {
    try {
        publicEntryLimits();
        retainedEvidenceOnly();
        curvedCuspIsNotAnImpossibilityProof();
        std::cout << "Retained planar-line quality constraints, inactive/unknown/reflex evidence "
                     "and non-blocking curved cusp PASS.\n";
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
