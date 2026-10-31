#include "result_comparison.hpp"
#include "project_json.hpp"
#include "result_export.hpp"
#include "simulation_request.hpp"
#include "sweep_model.hpp"
#include <QCoreApplication>
#include <QFileDialog>
#include <QGridLayout>
#include <QGuiApplication>
#include <QHBoxLayout>
#include <QInputDialog>
#include <QJsonDocument>
#include <QJsonObject>
#include <QLabel>
#include <QLineEdit>
#include <QListWidget>
#include <QMessageBox>
#include <QPageSize>
#include <QPainter>
#include <QPdfWriter>
#include <QPushButton>
#include <QSaveFile>
#include <QScreen>
#include <QTabWidget>
#include <QTextDocument>
#include <QTextEdit>
#include <QTextOption>
#include <QVBoxLayout>
#include <algorithm>
#include <array>
#include <cmath>
#include <hornsculpt/acoustic_analysis.hpp>
#include <hornsculpt/plot_widget.hpp>
#include <limits>
#include <stdexcept>

namespace hornsculpt {
namespace {
QString stateLabel(const Project &p) {
    const auto &a = p.simulation;
    if (!a.result)
        return QObject::tr("Kein Ergebnis");
    if (!a.current)
        return QObject::tr("Veraltet");
    if (!a.resultComplete)
        return QObject::tr("Aktuell · Vollständigkeit unbekannt");
    return *a.resultComplete ? QObject::tr("Zum gespeicherten Stand aktuell · vollständig")
                             : QObject::tr("Teilstand");
}
std::vector<QPointF> response(const Project &p) {
    std::vector<QPointF> points;
    if (!p.simulation.result)
        return points;
    const double missing = std::numeric_limits<double>::quiet_NaN();
    for (const auto &f : p.simulation.result->frequencies) {
        double level = missing;
        if (!f.excitations.empty() && !f.excitations.front().observerPressure.empty() &&
            p.simulation.request && observerValid(*p.simulation.request, 0))
            level = soundPressureLevelDecibels(f.excitations.front().observerPressure.front());
        points.emplace_back(f.frequency, level);
    }
    return points;
}
std::array<std::vector<QPointF>, 4> comparisonCurves(const Project &p) {
    std::array<std::vector<QPointF>, 4> curves;
    if (!p.simulation.result || !p.simulation.request)
        return curves;
    const double missing = std::numeric_limits<double>::quiet_NaN();
    for (const auto &f : p.simulation.result->frequencies) {
        double level = missing, horizontal = missing, vertical = missing, impedance = missing;
        if (!f.excitations.empty()) {
            const auto &e = f.excitations.front();
            if (!e.observerPressure.empty() && observerValid(*p.simulation.request, 0))
                level = soundPressureLevelDecibels(e.observerPressure.front());
            const auto count =
                e.observerPressure.size() > 1 ? (e.observerPressure.size() - 1) / 2 : 0;
            if ((count == 37 || count == 73) && e.observerPressure.size() == 1 + 2 * count &&
                angularSampleCount(*p.simulation.request) == count) {
                for (std::size_t cut = 0; cut < 2; ++cut) {
                    std::vector<AcousticsAngularSample> samples;
                    for (std::size_t a = 0; a < count; ++a)
                        samples.push_back({(count == 73 ? -180. : -90.) + 5. * a,
                                           observerValid(*p.simulation.request, 1 + cut * count + a)
                                               ? soundPressureLevelDecibels(
                                                     e.observerPressure[1 + cut * count + a])
                                               : missing});
                    if (const auto width = beamwidthDegrees(samples, 6.))
                        (cut == 0 ? horizontal : vertical) = width->widthDegrees;
                }
            }
            if (!e.impedances.empty()) {
                const auto &z = e.impedances.front();
                const auto medium = p.simulation.request->medium;
                impedance = std::abs(normalizedVolumeImpedance(z.volumeImpedance, z.area,
                                                               medium.density, medium.soundSpeed));
            }
        }
        for (std::size_t i = 0; i < curves.size(); ++i)
            curves[i].emplace_back(
                f.frequency, std::array<double, 4>{level, horizontal, vertical, impedance}[i]);
    }
    return curves;
}
void appendJsonDifferences(const QJsonValue &a, const QJsonValue &b, const QString &path,
                           QString &out, int &count) {
    if (a == b || count >= 40)
        return;
    if (a.isObject() && b.isObject()) {
        const auto ao = a.toObject(), bo = b.toObject();
        for (auto i = ao.begin(); i != ao.end(); ++i)
            appendJsonDifferences(i.value(), bo.value(i.key()),
                                  path.isEmpty() ? i.key() : path + "." + i.key(), out, count);
        for (auto i = bo.begin(); i != bo.end(); ++i)
            if (!ao.contains(i.key()))
                appendJsonDifferences(QJsonValue(QJsonValue::Undefined), i.value(),
                                      path.isEmpty() ? i.key() : path + "." + i.key(), out, count);
        return;
    }
    auto value = [](const QJsonValue &v) {
        if (v.isArray())
            return QString::fromUtf8(QJsonDocument(v.toArray()).toJson(QJsonDocument::Compact));
        if (v.isObject())
            return QString::fromUtf8(QJsonDocument(v.toObject()).toJson(QJsonDocument::Compact));
        if (v.isNull() || v.isUndefined())
            return QStringLiteral("—");
        return v.toVariant().toString();
    };
    out += QObject::tr("%1: %2 → %3\n").arg(path, value(a), value(b));
    ++count;
}
QString parameterDifferences(const Project &reference, const Project &variant) {
    QString output;
    if (reference.source != variant.source)
        output += QObject::tr("Geometriequelle unterschiedlich.\n");
    if (reference.source == ProjectSource::draft && variant.source == ProjectSource::draft) {
        for (const auto &parameter : sweepParameters()) {
            if (parameter.key == QStringLiteral("meshTargetEdgeLength"))
                continue;
            const double a = sweepParameterValue(reference.parameters, parameter.key);
            const double b = sweepParameterValue(variant.parameters, parameter.key);
            if (a != b)
                output += QObject::tr("%1: %2 → %3 %4\n")
                              .arg(parameter.label)
                              .arg(a * parameter.displayScale, 0, 'g', 7)
                              .arg(b * parameter.displayScale, 0, 'g', 7)
                              .arg(parameter.unit);
        }
        int count = 0;
        QString raw;
        appendJsonDifferences(project_json::encodeDraft(reference),
                              project_json::encodeDraft(variant), {}, raw, count);
        if (!raw.isEmpty())
            output += QObject::tr("Form-/Transformationsdaten (SI):\n") + raw;
        if (count >= 40)
            output += QObject::tr("Weitere Unterschiede nicht angezeigt.\n");
    } else
        output += QObject::tr("CAD-Geometrie und Flächenzuordnung im jeweiligen Stand prüfen.\n");
    if (reference.meshing.targetEdgeLength != variant.meshing.targetEdgeLength)
        output += QObject::tr("Zielkante: %1 → %2 mm\n")
                      .arg(reference.meshing.targetEdgeLength * 1000)
                      .arg(variant.meshing.targetEdgeLength * 1000);
    if (reference.simulation.request && variant.simulation.request) {
        SimulationArchive referenceRequest, variantRequest;
        referenceRequest.request = reference.simulation.request;
        variantRequest.request = variant.simulation.request;
        const QJsonValue a = project_json::encodeSimulation(referenceRequest).value("request");
        const QJsonValue b = project_json::encodeSimulation(variantRequest).value("request");
        int count = 0;
        QString raw;
        appendJsonDifferences(a, b, {}, raw, count);
        if (!raw.isEmpty())
            output += QObject::tr("Auftragsunterschiede (SI):\n") + raw;
        if (count >= 40)
            output += QObject::tr("Weitere Auftragsunterschiede nicht angezeigt.\n");
    } else
        output += QObject::tr("Anregungs-/Beobachterdaten fehlen in mindestens einem Stand.\n");
    return output;
}
QString summary(const Project &p) {
    QString result = QObject::tr("Quelle: %1 · Projekt-/Geometrie-/Netzrevision: %2 / %3 / %4\n")
                         .arg(p.source == ProjectSource::draft ? QStringLiteral("Draft")
                                                               : QStringLiteral("FreeCAD"))
                         .arg(p.revision)
                         .arg(p.geometryRevision)
                         .arg(p.meshRevision);
    if (!p.simulation.current && p.simulation.result)
        result += QObject::tr("Hinweis: Die folgenden Parameter beschreiben den aktuellen Entwurf, "
                              "nicht notwendigerweise das veraltete Ergebnis.\n");
    if (p.source == ProjectSource::draft)
        result += QObject::tr("Länge: %1 mm · Wand: %2 mm · H/V-Halbwinkel: %3° / %4°\n")
                      .arg(p.parameters.length * 1000, 0, 'g', 7)
                      .arg(p.parameters.wallThickness * 1000, 0, 'g', 7)
                      .arg(p.parameters.horizontalAngle * 180 / std::acos(-1.), 0, 'g', 6)
                      .arg(p.parameters.verticalAngle * 180 / std::acos(-1.), 0, 'g', 6);
    if (p.simulation.result && p.simulation.result->geometryMode == BemGeometryMode::axisymmetric) {
        result += QObject::tr(
                      "Achsensymmetrisches Integrationsmodell: %1 Druckfreiheitsgrade, %2 Dreiecke "
                      "im Integrations-/Ansichtssnapshot; kein isotropes 3D-Qualitätsnetz.\n")
                      .arg(p.simulation.result->pressureDofCount)
                      .arg(p.simulation.mesh ? p.simulation.mesh->triangles.size() : 0);
    } else if (p.mesh) {
        result += QObject::tr("Projekt-Netz: %1 Dreiecke · Zielkante: %2 mm\n")
                      .arg(p.mesh->mesh.triangles.size())
                      .arg(p.meshing.targetEdgeLength * 1000, 0, 'g', 6);
    } else
        result += QObject::tr("Kein Projekt-Netz vorhanden.\n");
    result += QObject::tr("Ergebnis: %1 · Revision: %2 · Frequenzen: %3\n")
                  .arg(stateLabel(p))
                  .arg(p.simulation.revision)
                  .arg(p.simulation.result ? p.simulation.result->frequencies.size() : 0);
    if (p.simulation.request) {
        result += QObject::tr("Anregungen: %1 · Beobachter: %2 · Medium: %3 kg/m³, %4 m/s\n")
                      .arg(p.simulation.request->excitations.size())
                      .arg(p.simulation.request->observers.size())
                      .arg(p.simulation.request->medium.density, 0, 'g', 6)
                      .arg(p.simulation.request->medium.soundSpeed, 0, 'g', 6);
        result += QObject::tr("Anregungsmodus: %1 · Pegelreferenz: 20 µPa RMS\n")
                      .arg(p.simulation.request->coupling
                               ? QObject::tr("gekoppeltes Netzwerk")
                               : QObject::tr("vorgegebene Flächengeschwindigkeit"));
        const auto &request = *p.simulation.request;
        if (!request.excitations.empty() && !request.excitations.front().faces.empty()) {
            const auto &face = request.excitations.front().faces.front();
            result +=
                QObject::tr(
                    "Erste Anregung, Fläche %1: Geschwindigkeit %2 + i·%3 m/s RMS (%4 Flächen)\n")
                    .arg(face.faceId)
                    .arg(face.normalVelocity.real(), 0, 'g', 6)
                    .arg(face.normalVelocity.imag(), 0, 'g', 6)
                    .arg(request.excitations.front().faces.size());
        }
        result += request.farField
                      ? QObject::tr("Fernfeld auf %1 m Bezugsabstand normiert.\n")
                            .arg(request.farField->referenceDistance, 0, 'g', 6)
                      : QObject::tr("Endliche Beobachterpositionen; kein Fernfeldabstand.\n");
    } else
        result += QObject::tr("Anregungs- und Beobachterdaten fehlen.\n");
    if (p.simulation.result) {
        double maximumResidual = 0.;
        bool hasResidual = false;
        for (const auto &frequency : p.simulation.result->frequencies)
            for (const auto &excitation : frequency.excitations)
                if (std::isfinite(excitation.boundaryResidual)) {
                    maximumResidual = std::max(maximumResidual, excitation.boundaryResidual);
                    hasResidual = true;
                }
        result += hasResidual ? QObject::tr("Größtes dokumentiertes Randresiduum: %1\n")
                                    .arg(maximumResidual, 0, 'g', 5)
                              : QObject::tr("Randresiduum nicht vorhanden.\n");
    }
    result += QObject::tr("Netz-/Quadraturkonvergenz: nicht aus diesem Bericht nachgewiesen.\n");
    return result;
}
} // namespace

bool canCaptureComparison(const Project &project) noexcept {
    const auto &a = project.simulation;
    return a.result && a.request && a.mesh && a.current && a.resultComplete == std::optional{true};
}

void writeCustomerReport(const QString &path, const Project &project, const QImage &currentView,
                         const QString &standName) {
    QSaveFile file(path);
    if (!file.open(QIODevice::WriteOnly))
        throw std::runtime_error(file.errorString().toStdString());
    {
        QPdfWriter pdf(&file);
        pdf.setPageSize(QPageSize(QPageSize::A4));
        pdf.setResolution(96);
        pdf.setTitle(QObject::tr("HornSculpt Ergebnisbericht"));
        QPainter painter(&pdf);
        if (!painter.isActive())
            throw std::runtime_error("PDF-Ausgabe konnte nicht gestartet werden.");
        const int width = pdf.width();
        QFont title = painter.font();
        title.setPointSize(17);
        title.setBold(true);
        painter.setFont(title);
        painter.drawText(QRect(0, 0, width, 40), Qt::AlignLeft,
                         QObject::tr("HornSculpt Ergebnisbericht"));
        QFont body = painter.font();
        body.setPointSize(10);
        body.setBold(false);
        painter.setFont(body);
        const QString intro =
            QObject::tr(
                "Version %1 · Stand: %3\n%2\nDie lineare Simulation belegt keine thermische "
                "Belastbarkeit oder Hochpegelfestigkeit. Fehlende Beobachterdaten bleiben Lücken; "
                "Netz- und Quadraturkonvergenz sind gesondert zu prüfen.")
                .arg(QCoreApplication::applicationVersion(), summary(project),
                     standName.isEmpty() ? QObject::tr("Aktueller Projektstand") : standName);
        QTextDocument introDocument;
        introDocument.setDocumentMargin(0);
        introDocument.setDefaultFont(body);
        QTextOption introOption;
        introOption.setWrapMode(QTextOption::WrapAtWordBoundaryOrAnywhere);
        introDocument.setDefaultTextOption(introOption);
        introDocument.setPlainText(intro);
        introDocument.setTextWidth(width);
        const QRect introBounds(0, 48, width,
                                static_cast<int>(std::ceil(introDocument.size().height())));
        if (introDocument.idealWidth() > width + 1. || introBounds.bottom() >= pdf.height())
            throw std::runtime_error("PDF-Einleitung passt nicht auf die Seite.");
        painter.save();
        painter.translate(0, introBounds.top());
        introDocument.drawContents(&painter);
        painter.restore();
        const int imageTop = introBounds.bottom() + 20;
        if (!currentView.isNull()) {
            const QRect imageRect(0, imageTop, width, 240);
            const QSize scaled = currentView.size().scaled(imageRect.size(), Qt::KeepAspectRatio);
            const QRect target(imageRect.center().x() - scaled.width() / 2,
                               imageRect.center().y() - scaled.height() / 2, scaled.width(),
                               scaled.height());
            painter.drawImage(target, currentView);
        } else
            painter.drawText(QRect(0, imageTop, width, 40),
                             QObject::tr("Keine Modellansicht übergeben."));
        PlotWidget plot(QObject::tr("Erster Beobachter · erste Anregung"),
                        QObject::tr("Frequenz [Hz]"), QObject::tr("Pegel [dB SPL]"));
        plot.resize(width, 280);
        plot.setFrequencyAxis();
        plot.setSeries(
            {PlotSeries{standName.isEmpty() ? QObject::tr("Aktueller Projektstand") : standName,
                        QColor(30, 100, 180), response(project)}});
        painter.save();
        painter.translate(0, imageTop + (currentView.isNull() ? 50 : 260));
        plot.render(&painter);
        painter.restore();
        pdf.newPage();
        painter.setFont(body);
        painter.drawText(QRect(0, 0, width, 48), Qt::TextWordWrap,
                         QObject::tr("Erste Anregung · H/V-Winkel bei −6 dB und normierte erste "
                                     "Flächenimpedanz. Fehlende Daten bleiben Lücken."));
        const auto values = comparisonCurves(project);
        const QString names[] = {QObject::tr("H-Abstrahlwinkel"), QObject::tr("V-Abstrahlwinkel"),
                                 QObject::tr("Erste Flächenimpedanz")};
        const QString units[] = {QObject::tr("Winkel [°]"), QObject::tr("Winkel [°]"),
                                 QObject::tr("Betrag normiert")};
        for (std::size_t i = 0; i < 3; ++i) {
            PlotWidget chart(names[i], QObject::tr("Frequenz [Hz]"), units[i]);
            chart.resize(width, 230);
            chart.setFrequencyAxis();
            chart.setSeries(
                {PlotSeries{QObject::tr("Projektstand"), QColor(30, 100, 180), values[i + 1]}});
            painter.save();
            painter.translate(0, 55 + static_cast<int>(i) * 240);
            chart.render(&painter);
            painter.restore();
        }
        QString cta = QObject::tr("CTA-2034: keine auswertbaren Fernfeld-Richtungsdaten.");
        if (project.simulation.request && project.simulation.result) {
            try {
                const auto rows =
                    ctaResults(*project.simulation.request, *project.simulation.result, 0);
                std::size_t valid = 0;
                for (const auto &row : rows)
                    if (row.curves.soundPowerDb)
                        ++valid;
                cta = QObject::tr("CTA-2034 Schallleistung: %1 auswertbare Frequenzen; "
                                  "70-Richtungs-Näherung, Raumkurve ist eine Schätzung.")
                          .arg(valid);
            } catch (const std::invalid_argument &) {
            }
        }
        painter.drawText(QRect(0, 800, width, 100), Qt::TextWordWrap, cta);
        pdf.newPage();
        painter.setFont(body);
        painter.drawText(QRect(0, 0, width, 70), Qt::TextWordWrap,
                         QObject::tr("CTA-2034 · 70-Richtungs-Näherung. Raumkurve ist eine "
                                     "Schätzung; fehlende Richtungen ergeben Lücken."));
        std::vector<QPointF> onAxis, listening, power;
        if (project.simulation.request && project.simulation.result) {
            try {
                const auto rows =
                    ctaResults(*project.simulation.request, *project.simulation.result, 0);
                const double missing = std::numeric_limits<double>::quiet_NaN();
                for (const auto &row : rows) {
                    onAxis.emplace_back(row.frequency, row.curves.onAxisDb.value_or(missing));
                    listening.emplace_back(row.frequency,
                                           row.curves.listeningWindowDb.value_or(missing));
                    power.emplace_back(row.frequency, row.curves.soundPowerDb.value_or(missing));
                }
            } catch (const std::invalid_argument &) {
            }
        }
        PlotWidget ctaPlot(QObject::tr("CTA-2034 · erste Anregung"), QObject::tr("Frequenz [Hz]"),
                           QObject::tr("Pegel [dB SPL]"));
        ctaPlot.resize(width, 400);
        ctaPlot.setFrequencyAxis();
        ctaPlot.setSeries(
            {{QObject::tr("On-Axis"), QColor(30, 100, 180), onAxis},
             {QObject::tr("Listening Window"), QColor(50, 145, 80), listening},
             {QObject::tr("Schallleistung (Näherung)"), QColor(210, 100, 30), power}});
        painter.save();
        painter.translate(0, 80);
        ctaPlot.render(&painter);
        painter.restore();
        painter.end();
    }
    if (!file.commit())
        throw std::runtime_error(file.errorString().toStdString());
}

struct ResultComparisonDialog::Implementation {
    Project current;
    QImage image;
    QListWidget *list{};
    std::array<PlotWidget *, 4> plots{};
    QTextEdit *details{};
};

ResultComparisonDialog::ResultComparisonDialog(Project current, QImage currentView, QWidget *parent)
    : QDialog(parent), implementation_(std::make_shared<Implementation>()) {
    auto &s = *implementation_;
    s.current = std::move(current);
    s.image = std::move(currentView);
    setWindowTitle(tr("Ergebnisstände vergleichen"));
    const auto *screen =
        parent && parent->screen() ? parent->screen() : QGuiApplication::primaryScreen();
    const QSize available = screen ? screen->availableGeometry().size() : QSize(853, 533);
    resize(std::min(980, available.width()), std::min(650, available.height()));
    auto *layout = new QVBoxLayout(this);
    auto *hint = new QLabel(tr("Stände enthalten ihre berechneten Daten. Veraltete Ergebnisse "
                               "bleiben ausdrücklich gekennzeichnet."),
                            this);
    hint->setWordWrap(true);
    layout->addWidget(hint);
    auto *row = new QHBoxLayout;
    layout->addLayout(row, 1);
    s.list = new QListWidget(this);
    s.list->setSelectionMode(QAbstractItemView::ExtendedSelection);
    row->addWidget(s.list, 1);
    auto *right = new QVBoxLayout;
    row->addLayout(right, 3);
    auto *tabs = new QTabWidget(this);
    tabs->setUsesScrollButtons(true);
    tabs->setElideMode(Qt::ElideRight);
    const QString titles[] = {tr("Erster Beobachter"), tr("H-Abstrahlwinkel −6 dB"),
                              tr("V-Abstrahlwinkel −6 dB"), tr("Erste Flächenimpedanz")};
    const QString units[] = {tr("Pegel [dB SPL]"), tr("Winkel [°]"), tr("Winkel [°]"),
                             tr("Betrag normiert")};
    for (std::size_t i = 0; i < s.plots.size(); ++i) {
        s.plots[i] = new PlotWidget(titles[i], tr("Frequenz [Hz]"), units[i], tabs);
        s.plots[i]->setFrequencyAxis();
        tabs->addTab(s.plots[i], titles[i]);
    }
    right->addWidget(tabs, 2);
    s.details = new QTextEdit(this);
    s.details->setReadOnly(true);
    right->addWidget(s.details, 1);
    auto *buttons = new QGridLayout;
    layout->addLayout(buttons);
    auto *add = new QPushButton(tr("Aktuellen Stand merken …"), this);
    auto *rename = new QPushButton(tr("Umbenennen …"), this);
    auto *remove = new QPushButton(tr("Entfernen"), this);
    auto *report = new QPushButton(tr("PDF-Bericht …"), this);
    auto *close = new QPushButton(tr("Schließen"), this);
    buttons->addWidget(add, 0, 0);
    buttons->addWidget(rename, 0, 1);
    buttons->addWidget(remove, 0, 2);
    buttons->addWidget(report, 1, 0);
    buttons->addWidget(close, 1, 2);
    auto refresh = [this] {
        auto &s = *implementation_;
        const auto selected = s.list->selectedItems();
        std::array<std::vector<PlotSeries>, 4> curves;
        QString description;
        const Project *reference = nullptr;
        const QColor colors[] = {QColor(30, 100, 180), QColor(210, 100, 30), QColor(50, 145, 80),
                                 QColor(155, 70, 155)};
        for (auto *item : selected) {
            const int i = s.list->row(item);
            if (i < 0 || static_cast<std::size_t>(i) >= s.current.comparisonSnapshots.size())
                continue;
            const auto &entry = s.current.comparisonSnapshots[static_cast<std::size_t>(i)];
            if (!entry.project)
                continue;
            const auto values = comparisonCurves(*entry.project);
            const auto color = colors[curves[0].size() % 4];
            for (std::size_t j = 0; j < curves.size(); ++j)
                curves[j].push_back({entry.name, color, values[j]});
            if (!reference) {
                reference = entry.project.get();
                description +=
                    tr("Referenz: %1 (erster ausgewählter Stand in Listenreihenfolge)\n\n")
                        .arg(entry.name);
            }
            description += entry.name + QStringLiteral(" — ") + stateLabel(*entry.project) + '\n' +
                           summary(*entry.project);
            if (reference != entry.project.get())
                description += parameterDifferences(*reference, *entry.project);
            description += '\n';
        }
        for (std::size_t i = 0; i < curves.size(); ++i)
            s.plots[i]->setSeries(std::move(curves[i]));
        s.details->setPlainText(description);
    };
    auto fill = [this, refresh] {
        auto &s = *implementation_;
        s.list->clear();
        for (const auto &entry : s.current.comparisonSnapshots)
            s.list->addItem(entry.name + QStringLiteral(" — ") +
                            (entry.project ? stateLabel(*entry.project) : tr("Ungültig")));
        refresh();
    };
    connect(s.list, &QListWidget::itemSelectionChanged, this, refresh);
    connect(add, &QPushButton::clicked, this, [this, fill] {
        auto &s = *implementation_;
        if (!canCaptureComparison(s.current)) {
            QMessageBox::warning(this, tr("Vergleich"),
                                 tr("Ein aktuelles Ergebnis mit Netz und Auftrag sowie bestätigter "
                                    "Vollständigkeit ist erforderlich."));
            return;
        }
        bool ok = false;
        const auto name = QInputDialog::getText(this, tr("Vergleichsstand"), tr("Name:"),
                                                QLineEdit::Normal, {}, &ok)
                              .trimmed();
        if (!ok)
            return;
        const auto duplicate =
            std::any_of(s.current.comparisonSnapshots.begin(), s.current.comparisonSnapshots.end(),
                        [&](const auto &e) { return e.name == name; });
        if (name.isEmpty() || name.size() > 128 || duplicate ||
            s.current.comparisonSnapshots.size() >= 64) {
            QMessageBox::warning(this, tr("Vergleich"),
                                 tr("Name ist leer, zu lang, bereits vergeben oder die maximale "
                                    "Anzahl ist erreicht."));
            return;
        }
        Project snapshot = s.current;
        snapshot.comparisonSnapshots.clear();
        s.current.comparisonSnapshots.push_back(
            {name, std::make_shared<const Project>(std::move(snapshot))});
        if (onSnapshotsChanged)
            onSnapshotsChanged(s.current.comparisonSnapshots);
        fill();
    });
    connect(rename, &QPushButton::clicked, this, [this, fill] {
        auto &s = *implementation_;
        const int i = s.list->currentRow();
        if (i < 0)
            return;
        bool ok = false;
        const auto name =
            QInputDialog::getText(this, tr("Vergleichsstand"), tr("Neuer Name:"), QLineEdit::Normal,
                                  s.current.comparisonSnapshots[i].name, &ok)
                .trimmed();
        if (!ok)
            return;
        for (int j = 0; j < static_cast<int>(s.current.comparisonSnapshots.size()); ++j)
            if (j != i && s.current.comparisonSnapshots[j].name == name) {
                QMessageBox::warning(this, tr("Vergleich"), tr("Name bereits vergeben."));
                return;
            }
        if (name.isEmpty() || name.size() > 128) {
            QMessageBox::warning(this, tr("Vergleich"), tr("Ungültiger Name."));
            return;
        }
        s.current.comparisonSnapshots[i].name = name;
        if (onSnapshotsChanged)
            onSnapshotsChanged(s.current.comparisonSnapshots);
        fill();
    });
    connect(remove, &QPushButton::clicked, this, [this, fill] {
        auto &s = *implementation_;
        const int i = s.list->currentRow();
        if (i < 0)
            return;
        s.current.comparisonSnapshots.erase(s.current.comparisonSnapshots.begin() + i);
        if (onSnapshotsChanged)
            onSnapshotsChanged(s.current.comparisonSnapshots);
        fill();
    });
    connect(report, &QPushButton::clicked, this, [this] {
        const auto &entries = implementation_->current.comparisonSnapshots;
        QStringList choices{tr("Aktueller Projektstand")};
        for (const auto &entry : entries)
            choices.push_back(tr("Vergleichsstand: %1").arg(entry.name));
        bool accepted = false;
        const auto choice = QInputDialog::getItem(this, tr("PDF-Bericht"),
                                                  tr("Welcher Stand soll berichtet werden?"),
                                                  choices, 0, false, &accepted);
        if (!accepted)
            return;
        const int selected = choices.indexOf(choice) - 1;
        const auto path =
            QFileDialog::getSaveFileName(this, tr("PDF-Bericht speichern"), {}, tr("PDF (*.pdf)"));
        if (path.isEmpty())
            return;
        const Project *target = selected >= 0 && static_cast<std::size_t>(selected) < entries.size()
                                    ? entries[static_cast<std::size_t>(selected)].project.get()
                                    : &implementation_->current;
        if (!target)
            return;
        const QImage image =
            target == &implementation_->current ? implementation_->image : QImage{};
        const QString standName =
            selected >= 0 && static_cast<std::size_t>(selected) < entries.size()
                ? entries[static_cast<std::size_t>(selected)].name
                : tr("Aktueller Projektstand");
        try {
            writeCustomerReport(path.endsWith(".pdf", Qt::CaseInsensitive) ? path : path + ".pdf",
                                *target, image, standName);
        } catch (const std::exception &e) {
            QMessageBox::critical(this, tr("PDF-Bericht"), QString::fromUtf8(e.what()));
        }
    });
    connect(close, &QPushButton::clicked, this, &QDialog::accept);
    fill();
}
} // namespace hornsculpt
