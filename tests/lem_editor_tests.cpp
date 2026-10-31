#include "lem_editor.hpp"

#include <QApplication>
#include <QCheckBox>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QDragEnterEvent>
#include <QDropEvent>
#include <QGraphicsScene>
#include <QGraphicsView>
#include <QLabel>
#include <QLineEdit>
#include <QListWidget>
#include <QMimeData>
#include <QMouseEvent>
#include <QPushButton>
#include <QTableWidget>
#include <QWheelEvent>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numbers>
#include <stdexcept>

using namespace hornsculpt;
namespace {
void require(bool value, const char *message) {
    if (!value)
        throw std::runtime_error(message);
}
template <class Control> Control *control(LemEditor &editor, const char *name) {
    auto *result = editor.findChild<Control *>(QString::fromLatin1(name));
    require(result != nullptr, "Missing LEM control.");
    return result;
}
void mouse(QGraphicsView &view, QEvent::Type type, QPointF position, Qt::MouseButton button,
           Qt::MouseButtons buttons, Qt::KeyboardModifiers modifiers = Qt::NoModifier) {
    QMouseEvent event(type, position, view.viewport()->mapToGlobal(position.toPoint()), button,
                      buttons, modifiers);
    QApplication::sendEvent(view.viewport(), &event);
    QApplication::processEvents();
}
void clickScene(QGraphicsView &view, QPointF position) {
    const auto screen = view.mapFromScene(position);
    mouse(view, QEvent::MouseButtonPress, screen, Qt::LeftButton, Qt::LeftButton);
    mouse(view, QEvent::MouseButtonRelease, screen, Qt::LeftButton, Qt::NoButton);
}
void selectElement(LemEditor &editor, std::size_t index) {
    auto *view = control<QGraphicsView>(editor, "lemCircuitView");
    const auto position = editor.network().elements.at(index).position;
    view->centerOn(position.x + 100, position.y + 35);
    clickScene(*view, {position.x + 100, position.y + 18});
}
template <class Component> std::size_t findComponent(const LemEditor &editor) {
    const auto &parts = editor.network().elements;
    const auto found = std::find_if(parts.begin(), parts.end(), [](const auto &value) {
        return std::holds_alternative<Component>(value.component);
    });
    require(found != parts.end(), "Expected circuit component is missing.");
    return static_cast<std::size_t>(found - parts.begin());
}
void chooseNode(LemEditor &editor, const char *name, lem::NodeId node) {
    auto *combo = control<QComboBox>(editor, name);
    const auto index = combo->findData(node.value);
    require(index >= 0, "Compatible terminal node is unavailable.");
    combo->setCurrentIndex(index);
    require(QMetaObject::invokeMethod(combo, "activated", Q_ARG(int, index)),
            "Cannot activate terminal connection.");
    QApplication::processEvents();
}

void interfaceWorkflow() {
    LemEditor editor;
    editor.resize(1300, 700);
    editor.setFaces({{7, "Driven", BoundaryCondition::driven},
                     {8, "Virtual gap", BoundaryCondition::wall},
                     {9, "Radiation aperture", BoundaryCondition::wall, false}});
    editor.setSelectedFaces({7, 8, 9});
    editor.show();
    QApplication::processEvents();
    control<QPushButton>(editor, "lemAddInterfaceCircuit")->click();
    require(editor.network().elements.size() == 2,
            "The interface action did not create a real chamber/port network");
    const auto portIndex = findComponent<lem::BemPort>(editor);
    const auto cavityIndex = findComponent<lem::Cavity>(editor);
    const auto port = [&]() -> const lem::BemPort & {
        return std::get<lem::BemPort>(editor.network().elements[portIndex].component);
    };
    require(port().kind == lem::BemPortKind::interfacePort &&
                port().faceIds == std::vector<std::uint32_t>{8} && !port().needsAssignmentReview,
            "The interface action assigned Driven/aperture faces or lost its role");
    require(std::get<lem::Cavity>(editor.network().elements[cavityIndex].component).volume == 0,
            "A new chamber invented a positive volume");
    bool rejected = false;
    try {
        lem::validateNetwork(editor.network());
    } catch (const std::invalid_argument &) {
        rejected = true;
    }
    require(rejected, "The unfinished chamber was accepted for simulation");
    selectElement(editor, cavityIndex);
    auto *volume = control<QDoubleSpinBox>(editor, "lemVolume");
    require(volume->suffix().contains("cm³") &&
                volume->specialValueText().contains("nicht gesetzt"),
            "Chamber volume units or the unset value are hidden");
    volume->setValue(5);
    require(
        std::abs(std::get<lem::Cavity>(editor.network().elements[cavityIndex].component).volume -
                 5e-6) < 1e-16,
        "cm³ chamber input did not convert to SI");
    lem::validateNetwork(editor.network());
    selectElement(editor, portIndex);
    auto *faceList = control<QListWidget>(editor, "lemPortFaces");
    require(faceList->count() == 1 && faceList->item(0)->data(Qt::UserRole).toUInt() == 8,
            "The interface checklist offers Driven or aperture surfaces");
    require(!editor.findChild<QLineEdit *>("lemPortVelocityWeights"),
            "A uniform interface offers Driven velocity weights");
    require(control<QComboBox>(editor, "lemPortPolarity")->currentText().contains("BEM-Frontluft"),
            "Interface volume-flow direction is not described as entering front air");
    auto review = editor.network();
    auto &lost = std::get<lem::BemPort>(review.elements[portIndex].component);
    lost.faceIds.clear();
    lost.needsAssignmentReview = true;
    editor.setNetwork(review);
    selectElement(editor, portIndex);
    require(control<QLabel>(editor, "lemPortAssignmentReview")->text().contains("prüfen"),
            "Lost virtual-interface assignment is not visible");
    editor.setSelectedFaces({7, 9});
    require(port().needsAssignmentReview,
            "Merely selecting surfaces cleared required assignment review");
    control<QPushButton>(editor, "lemTakeSelectedFaces")->click();
    require(port().needsAssignmentReview && port().faceIds.empty(),
            "Taking incompatible/aperture surfaces cleared assignment review");
    editor.setSelectedFaces({8});
    control<QPushButton>(editor, "lemTakeSelectedFaces")->click();
    require(!port().needsAssignmentReview && port().faceIds == std::vector<std::uint32_t>{8},
            "Explicit valid reassignment did not confirm the interface");
    const auto solved = lem::solve(editor.network(), 1000, {{{18000, 3000}}, {{2, -.5}}});
    editor.setSolution(solved);
    const auto *results = control<QTableWidget>(editor, "lemPortResults");
    require(
        results->rowCount() == 1 && results->item(0, 1) &&
            std::abs(results->item(0, 1)->text().toDouble() - std::abs(solved.ports[0].pressure)) <
                1e-5 * std::max(1., std::abs(solved.ports[0].pressure)) &&
            std::abs(results->item(0, 2)->text().toDouble() -
                     1e6 * std::abs(solved.ports[0].volumeVelocity)) <
                1e-5 * std::max(1., 1e6 * std::abs(solved.ports[0].volumeVelocity)),
        "Actual coupled pressure/volume-flow results are not displayed with correct units");
    editor.setFaces({{8, "Role changed", BoundaryCondition::driven}});
    selectElement(editor, portIndex);
    faceList = control<QListWidget>(editor, "lemPortFaces");
    require(faceList->count() == 1 && faceList->item(0)->text().contains("nicht verfügbar"),
            "A virtual interface with changed geometry role disappeared from assignment review");
}

void weightedDriveWorkflow() {
    LemEditor editor;
    editor.resize(1300, 700);
    editor.setFaces(
        {{7, "Dome"}, {8, "Surround"}, {9, "Other Driven"}, {10, "Gap", BoundaryCondition::wall}});
    lem::Network network;
    network.nodes = {{{1}, lem::Domain::acoustic, false, "Front", {}},
                     {{2}, lem::Domain::acoustic, true, "Reference", {}}};
    lem::BemPort port{{{1}, {2}}, {7, 8}, 1};
    port.velocityWeights = {{8, .5}};
    network.elements = {{{1}, "Weighted drive", {}, port},
                        {{2}, "Cavity", {330, 0}, lem::Cavity{{{1}, {2}}, 5e-6}}};
    editor.setNetwork(network);
    editor.show();
    QApplication::processEvents();
    selectElement(editor, 0);
    const auto currentPort = [&]() -> const lem::BemPort & {
        return std::get<lem::BemPort>(editor.network().elements[0].component);
    };
    auto *input = control<QLineEdit>(editor, "lemPortVelocityWeights");
    require(input->text() == "8:0.5", "Loaded sparse drive weights are not editable");
    input->setText("7:1.25; 8:0.5");
    control<QPushButton>(editor, "lemApplyVelocityWeights")->click();
    require(currentPort().velocityWeights == std::map<std::uint32_t, double>{{7, 1.25}, {8, .5}} &&
                !lem::samePhysics(network, editor.network()),
            "Editing drive weights did not change the physical network");
    const auto accepted = currentPort().velocityWeights;
    for (const auto *invalid : {"7:0", "8:nan", "99:1", "7:1;7:2", "7:0,5", "7:1;"}) {
        input->setText(invalid);
        control<QPushButton>(editor, "lemApplyVelocityWeights")->click();
        require(currentPort().velocityWeights == accepted &&
                    !control<QLabel>(editor, "lemPortWeightsError")->text().isEmpty(),
                "Invalid drive weights were silently accepted or changed the physical network");
    }
    editor.setSelectedFaces({8, 9});
    control<QPushButton>(editor, "lemTakeSelectedFaces")->click();
    require(currentPort().faceIds == std::vector<std::uint32_t>{8, 9} &&
                currentPort().velocityWeights == std::map<std::uint32_t, double>{{8, .5}},
            "3D selection takeover lost retained weights, retained removed weights, or invented "
            "weights for new faces");
    auto *faces = control<QListWidget>(editor, "lemPortFaces");
    for (int row = 0; row < faces->count(); ++row)
        if (faces->item(row)->data(Qt::UserRole).toUInt() == 8)
            faces->item(row)->setCheckState(Qt::Unchecked);
    require(currentPort().faceIds == std::vector<std::uint32_t>{9} &&
                currentPort().velocityWeights.empty(),
            "Removing an assigned face retained its velocity weight");
    control<QComboBox>(editor, "lemPortKind")->setCurrentIndex(1);
    require(currentPort().kind == lem::BemPortKind::interfacePort &&
                currentPort().velocityWeights.empty() &&
                !editor.findChild<QLineEdit *>("lemPortVelocityWeights"),
            "Switching to uniform interface retained Driven weights or its weight editor");
}

void acousticSlitWorkflow() {
    LemEditor editor;
    editor.resize(1300, 700);
    editor.show();
    QApplication::processEvents();
    auto *palette = control<QListWidget>(editor, "lemPalette");
    int slitRow = -1;
    for (int row = 0; row < palette->count(); ++row)
        if (palette->item(row)->data(Qt::UserRole).toByteArray() == "a-slit")
            slitRow = row;
    require(slitRow >= 0, "The component palette has no typed slit channel");
    palette->setCurrentRow(slitRow);
    control<QPushButton>(editor, "lemAddComponent")->click();
    const auto index = findComponent<lem::AcousticSlit>(editor);
    selectElement(editor, index);
    require(std::all_of(editor.network().nodes.begin(), editor.network().nodes.end(),
                        [](const auto &node) { return node.domain == lem::Domain::acoustic; }),
            "The slit action created non-acoustic terminals");
    const auto &unset = std::get<lem::AcousticSlit>(editor.network().elements[index].component);
    require(unset.length == 0 && unset.area == 0 && unset.gap == 0,
            "A generic slit invented geometry");
    control<QDoubleSpinBox>(editor, "lemSlitLength")->setValue(30);
    control<QDoubleSpinBox>(editor, "lemSlitArea")->setValue(150);
    control<QDoubleSpinBox>(editor, "lemSlitGap")->setValue(1.75);
    control<QDoubleSpinBox>(editor, "lemSlitViscosity")->setValue(18.4);
    const auto &slit = std::get<lem::AcousticSlit>(editor.network().elements[index].component);
    require(std::abs(slit.length - .03) < 1e-16 && std::abs(slit.area - 150e-6) < 1e-18 &&
                std::abs(slit.gap - .00175) < 1e-17 &&
                std::abs(slit.dynamicViscosity - 18.4e-6) < 1e-18,
            "Slit length, area, full gap or viscosity display units do not convert to SI");
    require(control<QDoubleSpinBox>(editor, "lemSlitGap")->toolTip().contains("keine Halbweite"),
            "The slit input hides its full-gap convention");
    auto saved = editor.network();
    saved.nodes.push_back({{100}, lem::Domain::electrical, true, "Electrical reference", {}});
    saved.nodes.push_back({{101}, lem::Domain::mechanical, true, "Mechanical reference", {}});
    editor.setNetwork(saved);
    selectElement(editor, index);
    auto *terminal = control<QComboBox>(editor, "lemTerminal0");
    require(terminal->findData(100) < 0 && terminal->findData(101) < 0 && editor.network() == saved,
            "Slit terminals permit cross-domain wiring or loading changes typed SI values");
}

void closedSlitCavityWorkflow() {
    LemEditor editor;
    editor.resize(1300, 700);
    editor.show();
    QApplication::processEvents();
    auto *palette = control<QListWidget>(editor, "lemPalette");
    int row = -1;
    for (int i = 0; i < palette->count(); ++i)
        if (palette->item(i)->data(Qt::UserRole).toByteArray() == "a-closed-slit")
            row = i;
    require(row >= 0, "Closed 1D annular chamber is absent from the component palette");
    palette->setCurrentRow(row);
    control<QPushButton>(editor, "lemAddComponent")->click();
    const auto index = findComponent<lem::ClosedSlitCavity>(editor);
    selectElement(editor, index);
    const auto &unset = std::get<lem::ClosedSlitCavity>(editor.network().elements[index].component);
    require(unset.length == 0 && unset.area == 0 && unset.gap == 0 &&
                std::all_of(editor.network().nodes.begin(), editor.network().nodes.end(),
                            [](const auto &node) { return node.domain == lem::Domain::acoustic; }),
            "Generic 1D chamber invented geometry or non-acoustic terminals");
    control<QDoubleSpinBox>(editor, "lemClosedSlitLength")->setValue(30);
    control<QDoubleSpinBox>(editor, "lemClosedSlitArea")->setValue(150);
    control<QDoubleSpinBox>(editor, "lemClosedSlitGap")->setValue(1.75);
    control<QDoubleSpinBox>(editor, "lemClosedSlitViscosity")->setValue(18.4);
    control<QDoubleSpinBox>(editor, "lemClosedSlitGamma")->setValue(1.4);
    control<QDoubleSpinBox>(editor, "lemClosedSlitPrandtl")->setValue(.71);
    const auto &chamber =
        std::get<lem::ClosedSlitCavity>(editor.network().elements[index].component);
    require(std::abs(chamber.length - .03) < 1e-16 && std::abs(chamber.area - 150e-6) < 1e-18 &&
                std::abs(chamber.gap - .00175) < 1e-17 &&
                std::abs(chamber.dynamicViscosity - 18.4e-6) < 1e-18 &&
                chamber.heatCapacityRatio == 1.4 && chamber.prandtlNumber == .71,
            "Closed 1D chamber parameters do not use the documented SI quantities");
    require(control<QLabel>(editor, "lemClosedSlitVolume")->text() == "4.5 cm³" &&
                control<QDoubleSpinBox>(editor, "lemClosedSlitGamma")->minimum() == 1,
            "Derived chamber volume or the gamma bound is incorrect");
    control<QDoubleSpinBox>(editor, "lemClosedSlitLength")->setValue(40);
    require(control<QLabel>(editor, "lemClosedSlitVolume")->text() == "6 cm³",
            "The derived volume did not update with the actual channel geometry");
    auto saved = editor.network();
    saved.nodes.push_back({{100}, lem::Domain::electrical, true, "Electrical", {}});
    saved.nodes.push_back({{101}, lem::Domain::mechanical, true, "Mechanical", {}});
    editor.setNetwork(saved);
    selectElement(editor, index);
    const auto *terminal = control<QComboBox>(editor, "lemTerminal1");
    require(terminal->findData(100) < 0 && terminal->findData(101) < 0 && editor.network() == saved,
            "Loading the closed chamber changes SI values or exposes cross-domain wiring");
}
} // namespace

int main(int argc, char **argv) {
    QApplication app(argc, argv);
    try {
        int changes = 0;
        int enableChanges = 0;
        LemEditor editor;
        editor.resize(1300, 700);
        editor.setFaces({{7, "Hals / Driven"}, {9, "Zweiter Hals / Driven"}});
        editor.onChanged = [&] { ++changes; };
        editor.onEnabledChanged = [&](bool) { ++enableChanges; };
        editor.show();
        QApplication::processEvents();
        editor.setCouplingEnabled(true);
        require(changes == 0 && enableChanges == 0 && editor.couplingEnabled(),
                "Loading the coupling toggle emits an edit.");
        control<QCheckBox>(editor, "lemEnabled")->click();
        require(enableChanges == 1 && !editor.couplingEnabled(),
                "The coupling toggle has no functional callback.");

        control<QPushButton>(editor, "lemAddDriverCircuit")->click();
        require(changes == 1 && editor.network().elements.size() == 4,
                "Adding a driver circuit did not produce a real source/driver/cavity/BEM network.");
        lem::validateNetwork(editor.network());
        const auto driverIndex = findComponent<lem::Driver>(editor);
        selectElement(editor, driverIndex);
        auto *mass = control<QDoubleSpinBox>(editor, "lemDriverMovingMass");
        require(std::abs(mass->value() - 4) < 1e-12 && mass->toolTip().contains("Mmd"),
                "Driver mass is not displayed in grams with the Mmd convention.");
        mass->setValue(3.5);
        control<QDoubleSpinBox>(editor, "lemDriverPistonArea")->setValue(53.2);
        const auto &parameters =
            std::get<lem::Driver>(editor.network().elements[driverIndex].component).parameters;
        require(std::abs(parameters.movingMass - .0035) < 1e-12 &&
                    std::abs(parameters.pistonArea - .00532) < 1e-12,
                "Driver display units do not convert to SI.");

        const auto sourceIndex = findComponent<lem::EffortSource>(editor);
        selectElement(editor, sourceIndex);
        control<QDoubleSpinBox>(editor, "lemAmplitude")->setValue(12);
        control<QDoubleSpinBox>(editor, "lemPhase")->setValue(90);
        const auto source =
            std::get<lem::EffortSource>(editor.network().elements[sourceIndex].component);
        require(std::abs(source.amplitude - lem::Complex{0, 12}) < 1e-10,
                "RMS amplitude and phase are not saved as a complex source.");

        const auto portIndex = findComponent<lem::BemPort>(editor);
        selectElement(editor, portIndex);
        std::vector<std::uint32_t> highlighted;
        editor.onFacesSelected = [&](auto faces) { highlighted = std::move(faces); };
        auto *faceList = control<QListWidget>(editor, "lemPortFaces");
        require(faceList->count() == 2 && faceList->item(0)->checkState() == Qt::Checked,
                "The preset does not bind an available BEM face.");
        faceList->item(1)->setCheckState(Qt::Checked);
        require(
            highlighted == std::vector<std::uint32_t>{7, 9} &&
                std::get<lem::BemPort>(editor.network().elements[portIndex].component).faceIds ==
                    highlighted,
            "BEM multi-face assignment or highlighting is disconnected.");
        editor.setFaces({{9, "Zweiter Hals / Driven"}});
        faceList = control<QListWidget>(editor, "lemPortFaces");
        require(faceList->count() == 2 && faceList->item(1)->text().contains("nicht verfügbar") &&
                    faceList->item(1)->checkState() == Qt::Checked,
                "Unavailable face assignments are silently lost.");
        faceList->item(1)->setCheckState(Qt::Unchecked);
        require(std::get<lem::BemPort>(editor.network().elements[portIndex].component).faceIds ==
                    std::vector<std::uint32_t>{9},
                "Unavailable BEM faces cannot be removed explicitly.");

        auto *view = control<QGraphicsView>(editor, "lemCircuitView");
        const auto original = editor.network().elements[driverIndex].position;
        view->centerOn(original.x + 100, original.y + 18);
        const QPoint start = view->mapFromScene({original.x + 100, original.y + 18});
        const QPoint end = start + QPoint(45, -25);
        mouse(*view, QEvent::MouseButtonPress, start, Qt::LeftButton, Qt::LeftButton);
        mouse(*view, QEvent::MouseMove, end, Qt::NoButton, Qt::LeftButton);
        mouse(*view, QEvent::MouseButtonRelease, end, Qt::LeftButton, Qt::NoButton);
        const auto moved = editor.network().elements[driverIndex].position;
        require(std::hypot(moved.x - original.x, moved.y - original.y) > 1,
                "Dragging a circuit component does not persist its position.");
        const auto sourceBeforeMove = editor.network().elements[sourceIndex].position;
        view->centerOn(sourceBeforeMove.x + 100, sourceBeforeMove.y + 18);
        const auto groupStart =
            view->mapFromScene({sourceBeforeMove.x + 100, sourceBeforeMove.y + 18});
        mouse(*view, QEvent::MouseButtonPress, groupStart, Qt::LeftButton, Qt::LeftButton,
              Qt::ControlModifier);
        mouse(*view, QEvent::MouseButtonRelease, groupStart, Qt::LeftButton, Qt::NoButton,
              Qt::ControlModifier);
        mouse(*view, QEvent::MouseButtonPress, groupStart, Qt::LeftButton, Qt::LeftButton);
        mouse(*view, QEvent::MouseMove, groupStart + QPoint(-20, 25), Qt::NoButton, Qt::LeftButton);
        mouse(*view, QEvent::MouseButtonRelease, groupStart + QPoint(-20, 25), Qt::LeftButton,
              Qt::NoButton);
        const auto sourceAfterMove = editor.network().elements[sourceIndex].position;
        const auto driverAfterMove = editor.network().elements[driverIndex].position;
        require(std::abs((sourceAfterMove.x - sourceBeforeMove.x) - (driverAfterMove.x - moved.x)) <
                        1e-9 &&
                    std::abs((sourceAfterMove.y - sourceBeforeMove.y) -
                             (driverAfterMove.y - moved.y)) < 1e-9 &&
                    std::hypot(sourceAfterMove.x - sourceBeforeMove.x,
                               sourceAfterMove.y - sourceBeforeMove.y) > 1,
                "Moving a selected group fails to persist every component position.");
        const auto loaded = editor.network();
        const auto countBeforeLoad = changes;
        editor.setNetwork(loaded);
        require(editor.network() == loaded && changes == countBeforeLoad,
                "Loading a saved schematic changes its topology, layout, or revision.");

        QMimeData mime;
        mime.setData("application/x-hornsculpt-lem-component", "e-r");
        const QPoint dropPosition = view->viewport()->rect().center();
        QDragEnterEvent enter(dropPosition, Qt::CopyAction, &mime, Qt::LeftButton, Qt::NoModifier);
        QApplication::sendEvent(view->viewport(), &enter);
        QDropEvent drop(dropPosition, Qt::CopyAction, &mime, Qt::LeftButton, Qt::NoModifier);
        QApplication::sendEvent(view->viewport(), &drop);
        QApplication::processEvents();
        require(enter.isAccepted() && drop.isAccepted() && editor.network().elements.size() == 5,
                "Palette drag/drop does not add an editable component.");
        const auto resistorIndex = findComponent<lem::Resistance>(editor);
        selectElement(editor, resistorIndex);
        control<QDoubleSpinBox>(editor, "lemResistance")->setValue(1000);
        auto *positive = control<QComboBox>(editor, "lemTerminal0");
        for (const auto &node : editor.network().nodes)
            if (node.domain == lem::Domain::acoustic)
                require(positive->findData(node.id.value) == -1,
                        "Terminal selection permits cross-domain wiring.");
        chooseNode(editor, "lemTerminal0", source.terminals.positive);
        chooseNode(editor, "lemTerminal1", source.terminals.negative);
        const auto resistor =
            std::get<lem::Resistance>(editor.network().elements[resistorIndex].component);
        require(resistor.terminals == source.terminals && resistor.resistance == 1000,
                "Editing terminal selectors does not wire the real circuit.");
        chooseNode(editor, "lemTerminal0", {0});
        require(std::get<lem::Resistance>(editor.network().elements[resistorIndex].component)
                        .terminals.positive != source.terminals.positive,
                "Creating a new terminal node does not disconnect the previous wire.");
        const auto resistorPosition = editor.network().elements[resistorIndex].position;
        view->centerOn(resistorPosition.x, resistorPosition.y + 66);
        clickScene(*view, {resistorPosition.x, resistorPosition.y + 66});
        const auto sourcePosition = editor.network().elements[sourceIndex].position;
        view->centerOn(sourcePosition.x, sourcePosition.y + 66);
        clickScene(*view, {sourcePosition.x, sourcePosition.y + 66});
        require(std::get<lem::Resistance>(editor.network().elements[resistorIndex].component)
                        .terminals.positive == source.terminals.positive,
                "Clicking two compatible schematic terminals does not join their nodes.");
        lem::validateNetwork(editor.network());
        selectElement(editor, resistorIndex);

        const auto solved = lem::solve(editor.network(), 300, {{lem::Complex{100000, 12000}}, {}});
        editor.setSolution(solved);
        auto *results = control<QTableWidget>(editor, "lemResults");
        require(
            results->rowCount() == 1 && results->item(0, 5) &&
                std::abs(results->item(0, 5)->text().toDouble() - solved.drivers[0].copperPower) <
                    1e-4 * std::max(1.0, solved.drivers[0].copperPower),
            "Calculated driver copper power is not displayed.");
        require(std::abs(results->item(0, 3)->text().toDouble() -
                         std::sqrt(2.0) * std::abs(solved.drivers[0].displacement) * 1000) < 1e-4,
                "Displayed peak displacement uses the wrong RMS or length convention.");
        control<QDoubleSpinBox>(editor, "lemResistance")->setValue(900);
        require(control<QLabel>(editor, "lemResultStatus")->text().contains("veraltet") &&
                    results->rowCount() == 1,
                "Editing the circuit fails to retain and mark stale solver results.");
        control<QPushButton>(editor, "lemDelete")->click();
        require(editor.network().elements.size() == 4,
                "Deleting a selected component does not remove it from the network.");
        lem::validateNetwork(editor.network());

        const auto originalDiagram = editor.network();
        auto largeDiagram = originalDiagram;
        largeDiagram.elements.front().position.x -= 10000;
        editor.setNetwork(largeDiagram);
        control<QPushButton>(editor, "lemFit")->click();
        require(view->transform().m11() < .25,
                "The large-graph fixture did not exercise fitting below the manual zoom range.");
        const QPointF zoomPosition = view->viewport()->rect().center();
        QWheelEvent wheel(zoomPosition, view->viewport()->mapToGlobal(zoomPosition.toPoint()),
                          QPoint{}, QPoint(0, 120), Qt::NoButton, Qt::ControlModifier,
                          Qt::NoScrollPhase, false);
        QApplication::sendEvent(view->viewport(), &wheel);
        require(view->transform().m11() >= .25 && view->transform().m11() <= 3,
                "Ctrl+wheel cannot zoom in after fitting a large circuit below 25 percent.");
        editor.setNetwork(originalDiagram);

        interfaceWorkflow();
        weightedDriveWorkflow();
        acousticSlitWorkflow();
        closedSlitCavityWorkflow();

        if (argc > 1) {
            selectElement(editor, driverIndex);
            control<QPushButton>(editor, "lemFit")->click();
            QApplication::processEvents();
            require(editor.grab().save(QString::fromLocal8Bit(argv[1])),
                    "Could not save the requested editor screenshot.");
        }
        std::cout
            << "LEM editor: source/driver SI editing, wiring, drag/drop, movement, BEM faces, "
               "deletion, persistence callbacks, and real solver result display passed.\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
