#include "acoustic_regions_dialog.hpp"
#include <QApplication>
#include <QColor>
#include <QLabel>
#include <QListWidget>
#include <QPushButton>
#include <QTableWidget>
#include <iostream>
#include <stdexcept>

using namespace hornsculpt;
namespace {
void require(bool value, const char *message) {
    if (!value)
        throw std::runtime_error(message);
}
std::map<std::uint32_t, QString> faces() {
    return {{1, "Throat"}, {2, "Wall"}, {3, "Interface"}, {4, "Other"}};
}
std::vector<AcousticRegionDefinition> initial() {
    return {{42, "Interior", {1}, {3}}, {77, "Other", {2}, {3}}};
}
QTableWidget *table(AcousticRegionsDialog &d) {
    auto *t = d.findChild<QTableWidget *>("acousticRegionsTable");
    require(t, "Missing table");
    return t;
}
void identityAndCancel() {
    AcousticRegionsDialog d(initial(), faces());
    auto *t = table(d);
    t->item(0, 0)->setText("Renamed");
    d.accept();
    require(d.result() == QDialog::Accepted && d.regions().size() == 2 && d.regions()[0].id == 42 &&
                d.regions()[1].id == 77 && d.regions()[0].name == "Renamed",
            "Rename changed ID");
    AcousticRegionsDialog canceled(initial(), faces());
    table(canceled)->item(0, 0)->setText("discarded");
    canceled.reject();
    require(canceled.regions().empty() && canceled.result() == QDialog::Rejected,
            "Cancel published edits");
}
void assignmentReview() {
    auto original = initial();
    original[0].needsAssignmentReview = true;
    AcousticRegionsDialog canceled(original, faces());
    auto *item = table(canceled)->item(0, 0);
    require(item->foreground().color() == QColor(Qt::red) &&
                item->toolTip().contains(QStringLiteral("Flächenzuordnung prüfen")) &&
                table(canceled)->item(0, 1)->text() == "1",
            "Review indication missing");
    canceled.reject();
    require(original[0].needsAssignmentReview && canceled.regions().empty(),
            "Cancel lost review state");
    AcousticRegionsDialog accepted(original, faces());
    accepted.accept();
    require(accepted.result() == QDialog::Accepted && accepted.regions()[0].id == 42 &&
                !accepted.regions()[0].needsAssignmentReview &&
                accepted.regions()[0].faceIds == original[0].faceIds,
            "Confirmed review did not preserve identity and clear flag");
    AcousticRegionsDialog missing(original, {{2, "Wall"}, {3, "Interface"}});
    missing.accept();
    require(missing.result() != QDialog::Accepted && missing.regions().empty(),
            "Unknown face allowed review confirmation");
}
void invalid() {
    auto reject = [](int row, int col, const QString &value) {
        AcousticRegionsDialog d(initial(), faces());
        table(d)->item(row, col)->setText(value);
        d.accept();
        auto *error = d.findChild<QLabel *>("acousticRegionsError");
        require(d.result() != QDialog::Accepted && error && !error->text().isEmpty(),
                "Invalid edit accepted");
    };
    reject(1, 1, "1");
    reject(0, 1, "999");
    reject(0, 1, "1,1");
    reject(0, 2, "1");
    reject(0, 0, " ");
    AcousticRegionsDialog d(initial(), faces());
    auto *t = table(d);
    d.findChild<QPushButton *>("acousticRegionAdd")->click();
    t->item(2, 0)->setText("Third");
    t->item(2, 1)->setText("4");
    t->item(2, 2)->setText("3");
    d.accept();
    require(d.result() != QDialog::Accepted, "Third interface owner accepted");
}
void selectionAndNewId() {
    AcousticRegionsDialog d(initial(), faces());
    auto *t = table(d);
    t->setCurrentCell(1, 0);
    d.findChild<QPushButton *>("acousticRegionRemove")->click();
    d.findChild<QPushButton *>("acousticRegionAdd")->click();
    t->item(1, 0)->setText("New");
    auto *list = d.findChild<QListWidget *>("acousticRegionFaces");
    require(list && list->count() == 4, "Missing face list");
    list->item(3)->setSelected(true);
    d.findChild<QPushButton *>("acousticRegionUsePhysical")->click();
    d.accept();
    require(d.result() == QDialog::Accepted && d.regions()[0].id == 42 && d.regions()[1].id > 77 &&
                d.regions()[1].faceIds == std::vector<std::uint32_t>{4},
            "Face picker or new identity failed");
}
} // namespace
int main(int argc, char **argv) {
    QApplication app(argc, argv);
    try {
        identityAndCancel();
        assignmentReview();
        invalid();
        selectionAndNewId();
        std::cout << "Acoustic regions dialog tests PASS.\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
