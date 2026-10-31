#include "acoustic_regions_dialog.hpp"
#include <QApplication>
#include <QLabel>
#include <QListWidget>
#include <QPushButton>
#include <QTreeWidget>
#include <iostream>
#include <stdexcept>
using namespace hornsculpt;
namespace {
void require(bool value, const char *message) {
    if (!value)
        throw std::runtime_error(message);
}
AcousticModelDefinition fixture() {
    return {{{42, "Innen", BemDomainKind::interior, {1}}, {77, "Außen", BemDomainKind::exterior}},
            {{3, "Öffnung", 42, 77, {3}}},
            std::nullopt};
}
std::map<std::uint32_t, QString> faces() {
    return {{1, "A"}, {2, "B"}, {3, "C"}, {4, "D"}};
}
void tests() {
    auto baffle = fixture();
    baffle.exteriorBaffle = ExteriorBaffle{{0, 0, .1}, {0, 0, 1}, 1e-6};
    AcousticRegionsDialog baffleDialog(baffle, faces());
    auto *baffleTree = baffleDialog.findChild<QTreeWidget *>("acousticModelTree");
    require(baffleTree && baffleTree->topLevelItem(1)->childCount() >= 3 &&
                baffleTree->topLevelItem(1)->child(2)->text(0).contains("Schallwand") &&
                baffleTree->topLevelItem(2)->text(1).contains("Innen") &&
                baffleTree->topLevelItem(2)->text(1).contains("Außen"),
            "Baffle or named interface domains are missing from the model tree.");
    auto grouped = fixture();
    grouped.domains.front().faceIds = {1, 2};
    grouped.domains.front().groups = {{1, "Driver", {1}}};
    AcousticRegionsDialog groupEdit(grouped, faces());
    auto *groupTree = groupEdit.findChild<QTreeWidget *>("acousticModelTree");
    auto *groupAssign = groupEdit.findChild<QPushButton *>("acousticAssignFaces");
    groupTree->setCurrentItem(groupTree->topLevelItem(0)->child(2));
    groupEdit.selectFaces({1, 2});
    groupAssign->click();
    groupEdit.accept();
    require(groupEdit.result() == QDialog::Accepted &&
                groupEdit.model().domains.front().groups.front().faceIds ==
                    std::vector<std::uint32_t>({1, 2}),
            "Group multi-face edit failed.");
    AcousticRegionsDialog canceled(fixture(), faces());
    auto *tree = canceled.findChild<QTreeWidget *>("acousticModelTree");
    auto *assign = canceled.findChild<QPushButton *>("acousticAssignFaces");
    require(tree && assign && tree->topLevelItemCount() == 3, "Model tree missing entries.");
    tree->setCurrentItem(tree->topLevelItem(2));
    canceled.selectFaces({2, 4});
    assign->click();
    canceled.reject();
    require(canceled.model().domains.empty() && canceled.result() == QDialog::Rejected,
            "Cancel committed edits.");
    AcousticRegionsDialog accepted(fixture(), faces());
    tree = accepted.findChild<QTreeWidget *>("acousticModelTree");
    assign = accepted.findChild<QPushButton *>("acousticAssignFaces");
    tree->setCurrentItem(tree->topLevelItem(2));
    accepted.selectFaces({2, 4});
    assign->click();
    accepted.accept();
    require(accepted.result() == QDialog::Accepted &&
                accepted.model().interfaces.front().faceIds == std::vector<std::uint32_t>({2, 4}),
            "Multi-face interface edit failed.");
    require(accepted.model().interfaces.front().firstDomainId == 42 &&
                accepted.model().interfaces.front().secondDomainId == 77,
            "Interface identities changed.");
    AcousticRegionsDialog duplicate(fixture(), faces());
    tree = duplicate.findChild<QTreeWidget *>("acousticModelTree");
    assign = duplicate.findChild<QPushButton *>("acousticAssignFaces");
    tree->setCurrentItem(tree->topLevelItem(2));
    duplicate.selectFaces({1});
    assign->click();
    duplicate.accept();
    require(duplicate.result() != QDialog::Accepted,
            "Overlapping physical and interface faces accepted.");
    auto broken = fixture();
    broken.interfaces[0].secondDomainId = 99;
    AcousticRegionsDialog invalid(broken, faces());
    invalid.accept();
    require(invalid.result() != QDialog::Accepted, "Missing domain reference accepted.");
}
} // namespace
int main(int argc, char **argv) {
    QApplication app(argc, argv);
    try {
        tests();
        std::cout << "Acoustic model dialog tests PASS.\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
