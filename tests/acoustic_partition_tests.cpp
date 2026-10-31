#include "acoustic_partition_fixture.hpp"
#include "fixture_geometry.hpp"
#include <algorithm>
#include <cmath>
#include <hornsculpt/acoustic_partition.hpp>
#include <iostream>
#include <numbers>
#include <set>
#include <stdexcept>

using namespace hornsculpt;
namespace {
void require(bool condition, const char *message) {
    if (!condition)
        throw std::runtime_error(message);
}
template <class F> void rejects(F function) {
    try {
        function();
    } catch (const std::invalid_argument &) {
        return;
    }
    throw std::runtime_error("Invalid acoustic partition accepted.");
}

} // namespace
int main() {
    try {
        // One finite shell has an inward-oriented inner component and an
        // outward-oriented enclosing component. Their net signed volume stays
        // positive even if the inner component is incorrectly reversed.
        auto hollow = testing::sphere(0);
        auto outer = testing::sphere(0);
        const auto innerCount = static_cast<std::uint32_t>(hollow.triangles.size());
        const auto vertexOffset = static_cast<std::uint32_t>(hollow.vertices.size());
        for (auto v : outer.vertices)
            hollow.vertices.push_back({2 * v.x, 2 * v.y, 2 * v.z});
        for (auto triangle : outer.triangles) {
            for (auto &v : triangle.vertices)
                v += vertexOffset;
            triangle.faceId = 2;
            triangle.boundary = BoundaryCondition::wall;
            hollow.triangles.push_back(triangle);
        }
        BemSubdomain shell{1, "Schale", BemDomainKind::interior, {}};
        BemSubdomain exterior{2, "Aussen", BemDomainKind::exterior, {}};
        BemInterface sphereInterface{3, "Apertur", 1, 2, {}};
        for (std::uint32_t t = 0; t < innerCount; ++t)
            shell.panels.push_back({t, -1});
        for (auto t = innerCount; t < hollow.triangles.size(); ++t) {
            shell.panels.push_back({t, 1});
            exterior.panels.push_back({t, -1});
            sphereInterface.triangleIndices.push_back(t);
        }
        BemSubdomainModel shellModel{{shell, exterior}, {sphereInterface}};
        validateAcousticPartition(hollow, shellModel);
        auto inverted = shellModel;
        for (std::uint32_t t = 0; t < innerCount; ++t)
            inverted.domains[0].panels[t].normalSign = 1;
        rejects([&] { validateAcousticPartition(hollow, inverted); });
        const auto physical = testing::acousticCup();
        const auto sourceQuality = inspectMesh(physical);
        require(sourceQuality.watertight(), "Cup fixture is invalid.");
        require(sourceQuality.signedVolume > 0, "Cup fixture has reversed physical normals.");
        MeshingOptions options;
        options.targetEdgeLength = .02;
        const std::vector<AcousticRegionDefinition> regions{{7, "Hornluft", {2, 3}}};
        const auto split = partitionAcousticSurface(physical, regions, options);
        require(split.model.domains.size() == 2 && split.model.interfaces.size() == 1,
                "Cup did not form one interior, one exterior and one interface.");
        require(split.mesh.triangles.size() > physical.triangles.size(), "No interface generated.");
        require(std::equal(physical.triangles.begin(), physical.triangles.end(),
                           split.mesh.triangles.begin()),
                "Physical mesh tags or topology changed.");
        require(std::equal(physical.vertices.begin(), physical.vertices.end(),
                           split.mesh.vertices.begin()),
                "Physical coordinates were moved.");
        const auto &interface = split.model.interfaces.front();
        require(interface.firstDomainId == 7 && interface.secondDomainId != 7,
                "Wrong interface owners.");
        for (const auto index : interface.triangleIndices) {
            const auto &t = split.mesh.triangles.at(index);
            for (const auto vertex : t.vertices)
                require(std::abs(split.mesh.vertices.at(vertex).z - .08) < 1e-12,
                        "Interface left its explicit opening plane.");
            unsigned count = 0;
            int sum = 0;
            for (const auto &domain : split.model.domains)
                for (const auto &panel : domain.panels)
                    if (panel.triangleIndex == index) {
                        ++count;
                        sum += panel.normalSign;
                    }
            require(count == 2 && sum == 0, "Interface did not retain opposite domain sides.");
        }
        const auto interfaceFace =
            split.mesh.triangles.at(interface.triangleIndices.front()).faceId;
        const std::vector<AcousticRegionDefinition> explicitRegions{
            {7, "Hornluft", {2, 3}, {interfaceFace}}};
        const auto imported = partitionAcousticSurface(split.mesh, explicitRegions, options);
        require(imported.mesh == split.mesh && imported.model.interfaces.size() == 1,
                "An explicitly imported interface was remeshed or duplicated.");
        auto incomplete = split.mesh;
        incomplete.triangles.erase(incomplete.triangles.begin());
        rejects([&] { partitionAcousticSurface(incomplete, explicitRegions, options); });
        auto wrongRoles = explicitRegions;
        wrongRoles[0].faceIds.push_back(interfaceFace);
        rejects([&] { partitionAcousticSurface(split.mesh, wrongRoles, options); });
        auto three = physical;
        for (auto &t : three.triangles)
            if (t.faceId == 2) {
                double z = 0;
                for (auto vertex : t.vertices)
                    z += three.vertices[vertex].z / 3;
                if (z > .04)
                    t.faceId = 6;
            }
        const auto chain = partitionAcousticSurface(
            three,
            std::vector<AcousticRegionDefinition>{{7, "Hals", {2, 3}, {}}, {8, "Glocke", {6}, {}}},
            options);
        require(chain.model.domains.size() == 3 && chain.model.interfaces.size() == 2,
                "Interior/interior/exterior chain was not created.");
        bool innerLink = false, exteriorLink = false;
        for (const auto &link : chain.model.interfaces) {
            innerLink = innerLink || (link.firstDomainId == 7 && link.secondDomainId == 8);
            exteriorLink = exteriorLink || (link.firstDomainId == 8 && link.secondDomainId != 7);
        }
        require(innerLink && exteriorLink, "Wrong domain graph for the split horn.");
        rejects([&] { partitionAcousticSurface(physical, {}, options); });
        rejects([&] {
            partitionAcousticSurface(
                physical, std::vector<AcousticRegionDefinition>{{7, "bad", {999}}}, options);
        });
        rejects([&] {
            partitionAcousticSurface(
                physical, std::vector<AcousticRegionDefinition>{{7, "a", {2}}, {8, "b", {2}}},
                options);
        });
        rejects([&] {
            partitionAcousticSurface(
                physical, std::vector<AcousticRegionDefinition>{{7, "solid", {1, 5}}}, options);
        });
        auto broken = physical;
        broken.triangles.pop_back();
        rejects([&] { partitionAcousticSurface(broken, regions, options); });
        auto bent = physical;
        bent.vertices[4 * 16].z += .001;
        rejects([&] { partitionAcousticSurface(bent, regions, options); });
        auto strict = options;
        strict.quality.minimumAngleDegrees = 59;
        rejects([&] { partitionAcousticSurface(physical, regions, strict); });
        std::stop_source stop;
        stop.request_stop();
        bool cancelled = false;
        try {
            (void)partitionAcousticSurface(physical, regions, options, stop.get_token());
        } catch (const std::runtime_error &) {
            cancelled = true;
        }
        require(cancelled, "Partition ignored cancellation.");
        std::cout << "Acoustic partition: " << physical.triangles.size() << " physical, "
                  << interface.triangleIndices.size()
                  << " common interface triangles; quality, signs and cancellation PASS.\n";
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
