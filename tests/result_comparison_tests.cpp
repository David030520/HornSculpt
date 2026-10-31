#include "result_comparison.hpp"
#include <QApplication>
#include <QFile>
#include <QFileInfo>
#include <QListWidget>
#include <QTextEdit>
#include <QTemporaryDir>
#include <QTabWidget>
#include <hornsculpt/plot_widget.hpp>
#include <cmath>
#include <numbers>
#include <iostream>
#include <stdexcept>
using namespace hornsculpt;
int main(int argc, char **argv) {
    QApplication app(argc, argv);
    QCoreApplication::setApplicationVersion(QStringLiteral("QA-0.1.52"));
    try {
        QTemporaryDir directory;
        if (!directory.isValid()) throw std::runtime_error("Temporary directory unavailable");
        const bool external = argc == 4 && QString::fromLocal8Bit(argv[1]) == "--project";
        const auto path = external ? QString::fromLocal8Bit(argv[3]) : directory.filePath("report.pdf");
        const auto project = external ? readProject(QString::fromLocal8Bit(argv[2])) : Project{};
        writeCustomerReport(path, project, {}, external ? QStringLiteral("HF-Draft Referenz") : QString());
        writeCustomerReport(directory.filePath("long-name.pdf"), project, {}, QString(128, QChar('W')));
        QFile file(path);
        if (!file.open(QIODevice::ReadOnly) || !file.read(5).startsWith("%PDF-"))
            throw std::runtime_error("PDF report not written");
        ResultComparisonDialog dialog(Project{});
        if (dialog.windowTitle().isEmpty()) throw std::runtime_error("Comparison dialog unavailable");
        Project reference;
        reference.simulation.request = BemRequest{};
        reference.simulation.request->observers.push_back({0, 0, 1});
        for (int plane=0; plane<2; ++plane) for (int i=0; i<73; ++i) {
            const double angle = (-180. + 5.*i) * std::numbers::pi / 180.;
            reference.simulation.request->observers.push_back(plane == 0
                ? Vec3{std::sin(angle), 0, std::cos(angle)} : Vec3{0, std::sin(angle), std::cos(angle)});
        }
        BemExcitationResult excitation;
        for (int i=0; i<147; ++i) excitation.observerPressure.push_back({1, 0});
        BemFrequencyResult frequency; frequency.frequency = 1000; frequency.excitations.push_back(excitation);
        BemResult result; result.frequencies.push_back(frequency);
        reference.simulation.result = std::make_shared<const BemResult>(result);
        reference.simulation.mesh = std::make_shared<const SurfaceMesh>();
        reference.simulation.current = true;
        if (canCaptureComparison(reference)) throw std::runtime_error("Unknown completeness accepted");
        reference.simulation.resultComplete = false;
        if (canCaptureComparison(reference)) throw std::runtime_error("Partial result accepted");
        reference.simulation.resultComplete = true;
        if (!canCaptureComparison(reference)) throw std::runtime_error("Complete result rejected");
        Project variant = reference;
        variant.meshing.targetEdgeLength *= 1.1;
        variant.simulation.request->medium.density = 1.3;
        variant.simulation.request->farField = FarFieldObservation{};
        std::swap(variant.simulation.request->observers[2], variant.simulation.request->observers[3]);
        Project hostForTest;
        hostForTest.comparisonSnapshots = {{QStringLiteral("Referenz"), std::make_shared<const Project>(reference)},
                                           {QStringLiteral("Variante"), std::make_shared<const Project>(variant)}};
        ResultComparisonDialog comparison(std::move(hostForTest));
        auto *items = comparison.findChild<QListWidget *>();
        if (!items || items->count() != 2) throw std::runtime_error("Comparison items missing");
        items->item(0)->setSelected(true); items->item(1)->setSelected(true);
        auto *details = comparison.findChild<QTextEdit *>();
        if (!details || !details->toPlainText().contains(QStringLiteral("Referenz: Referenz")) ||
            !details->toPlainText().contains(QStringLiteral("Zielkante:")) ||
            !details->toPlainText().contains(QStringLiteral("density")) ||
            !details->toPlainText().contains(QStringLiteral("farField")))
            throw std::runtime_error("Comparison differences missing");
        const auto *tabs = comparison.findChild<QTabWidget *>();
        const auto *horizontal = tabs ? dynamic_cast<PlotWidget *>(tabs->widget(1)) : nullptr;
        if (!horizontal || horizontal->series().size() != 2 ||
            horizontal->series()[1].points.empty() || std::isfinite(horizontal->series()[1].points.front().y()))
            throw std::runtime_error("Invalid angular raster was used for beamwidth");
        Project reversed;
        reversed.comparisonSnapshots = {{QStringLiteral("Mit Fernfeld"), std::make_shared<const Project>(variant)},
                                         {QStringLiteral("Ohne Fernfeld"), std::make_shared<const Project>(reference)}};
        ResultComparisonDialog reverseComparison(std::move(reversed));
        auto *reverseItems = reverseComparison.findChild<QListWidget *>();
        reverseItems->item(0)->setSelected(true); reverseItems->item(1)->setSelected(true);
        if (!reverseComparison.findChild<QTextEdit *>()->toPlainText().contains(QStringLiteral("farField")))
            throw std::runtime_error("Removed request option missing from comparison");
        if (external) {
            auto first = project;
            auto second = project;
            first.comparisonSnapshots.clear(); second.comparisonSnapshots.clear();
            second.meshing.targetEdgeLength *= 1.1;
            second.simulation.current = false;
            auto host = project;
            host.comparisonSnapshots = {{QStringLiteral("Referenz"), std::make_shared<const Project>(first)},
                                        {QStringLiteral("Test: geänderter Entwurf, Ergebnis veraltet"), std::make_shared<const Project>(second)}};
            ResultComparisonDialog visual(std::move(host));
            visual.setAttribute(Qt::WA_DontShowOnScreen);
            visual.resize(853, 533);
            visual.show(); app.processEvents();
            auto *list = visual.findChild<QListWidget *>();
            if (!list || list->count() != 2) throw std::runtime_error("Comparison items missing");
            list->item(0)->setSelected(true); list->item(1)->setSelected(true);
            app.processEvents();
            if (!visual.grab().save(QFileInfo(path).absolutePath() + "/comparison-853.png"))
                throw std::runtime_error("Comparison capture failed");
        }
        std::cout << "PDF and comparison dialog smoke passed.\n";
        return 0;
    } catch (const std::exception &error) { std::cerr << error.what() << '\n'; return 1; }
}
