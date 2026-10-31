#include "step_export_dialog.hpp"
#include "numeric_input.hpp"

#include <QApplication>
#include <QComboBox>
#include <QFileDialog>
#include <QFileInfo>
#include <QFormLayout>
#include <QHBoxLayout>
#include <QLabel>
#include <QProgressBar>
#include <QPushButton>
#include <QSaveFile>
#include <QTimer>
#include <QVBoxLayout>
#include <atomic>
#include <hornsculpt/draft_step_export.hpp>
#include <limits>
#include <stdexcept>
#include <thread>
#include <utility>

namespace hornsculpt {
namespace {
struct ExportState {
    DraftStepExportResult result;
    QString error;
    std::atomic_bool finished{};
};
} // namespace

// A dialog can disappear immediately while Open CASCADE finishes its current
// operation. This controller owns the worker until completion, never captures
// GUI widgets in that worker and cannot publish without its live dialog.
class StepExportJob final : public QObject {
  public:
    StepExportJob(StepExportDialog *dialog, std::shared_ptr<const DraftGeometry> geometry,
                  DraftStepExportOptions options)
        : QObject(qApp), dialog_(dialog), state_(std::make_shared<ExportState>()) {
        worker_ = std::jthread(
            [geometry = std::move(geometry), options, state = state_](std::stop_token stop) {
                try {
                    state->result = exportDraftStep(*geometry, options, stop);
                } catch (const std::exception &error) {
                    state->error = QString::fromUtf8(error.what());
                }
                state->finished.store(true, std::memory_order_release);
            });
        connect(&poll_, &QTimer::timeout, this, [this] {
            if (!state_->finished.load(std::memory_order_acquire))
                return;
            poll_.stop();
            const bool cancelled = worker_.get_stop_token().stop_requested();
            worker_.join();
            if (dialog_ && dialog_->job_ == this)
                dialog_->finished(std::move(state_->result), state_->error, cancelled);
            deleteLater();
        });
        poll_.start(25);
    }

    void cancel() {
        worker_.request_stop();
    }

  private:
    QPointer<StepExportDialog> dialog_;
    std::shared_ptr<ExportState> state_;
    std::jthread worker_;
    QTimer poll_;
};

StepExportDialog::StepExportDialog(std::shared_ptr<const DraftGeometry> geometry, quint64 revision,
                                   QWidget *parent)
    : QDialog(parent), geometry_(std::move(geometry)), revision_(revision) {
    if (!geometry_)
        throw std::invalid_argument("STEP-Export benötigt eine gültige Draft-Geometrie.");
    setWindowTitle(tr("Draft als STEP exportieren"));
    setObjectName("stepExportDialog");
    setModal(false);
    resize(510, 290);
    auto *layout = new QVBoxLayout(this);
    auto *form = new QFormLayout;
    auto *revisionLabel = new QLabel(tr("Draft · Geometriestand %1").arg(revision_));
    revisionLabel->setObjectName("stepExportRevision");
    revisionLabel->setToolTip(
        tr("Unveränderlicher Geometriestand beim Öffnen dieses Dialogs. Spätere Draft-Änderungen "
           "gelten erst für einen neu geöffneten Export."));
    revisionLabel->setWhatsThis(revisionLabel->toolTip());
    layout->addWidget(revisionLabel);
    kind_ = new QComboBox;
    kind_->setObjectName("stepExportKind");
    kind_->addItem(tr("Hornkörper mit Wand"), static_cast<int>(DraftStepExportKind::solid));
    kind_->addItem(tr("Innere Hornfläche"), static_cast<int>(DraftStepExportKind::innerSurface));
    kind_->addItem(tr("H/V-Profile und Querschnitte"),
                   static_cast<int>(DraftStepExportKind::sections));
    kind_->setToolTip(
        tr("Hornkörper: Materialwand mit radialer Draft-Wandstärke und ringförmigen Hals-/"
           "Mundrändern; beide Luftöffnungen bleiben frei. Innenfläche: offene Horninnenwand. "
           "Konturen: vier Profile rechts/links/oben/unten und geschlossene axiale Querschnitte."));
    kind_->setWhatsThis(kind_->toolTip());
    form->addRow(tr("Umfang"), kind_);
    tolerance_ = new DecimalInput;
    tolerance_->setObjectName("stepExportTolerance");
    tolerance_->setDecimals(4);
    tolerance_->setRange(.001, 1);
    tolerance_->setSingleStep(.01);
    tolerance_->setValue(.05);
    tolerance_->setSuffix(tr(" mm"));
    tolerance_->setToolTip(
        tr("Zulässige Abweichung der CAD-Approximation von der finalen Draft-Geometrie an den "
           "adaptiven Prüfstellen in Millimetern. Kleinere Werte erhöhen Aufbauzeit und "
           "Dateigröße. Verletzungen an Prüfstellen werden abgewiesen; kein globaler "
           "Abstandsnachweis."));
    tolerance_->setWhatsThis(tolerance_->toolTip());
    form->addRow(tr("Formabweichung"), tolerance_);
    sections_ = new IntegerInput;
    sections_->setObjectName("stepExportSectionCount");
    sections_->setRange(2, 101);
    sections_->setValue(9);
    sections_->setToolTip(
        tr("Anzahl geschlossener Querschnittskurven, gleichmäßig entlang der Hornlänge verteilt, "
           "einschließlich Hals und Mund. Nur für den Konturenexport."));
    sections_->setWhatsThis(sections_->toolTip());
    sectionsLabel_ = new QLabel(tr("Querschnitte"));
    form->addRow(sectionsLabel_, sections_);
    layout->addLayout(form);
    destination_ = new QLabel;
    destination_->setObjectName("stepExportDestination");
    destination_->setTextFormat(Qt::PlainText);
    destination_->setWordWrap(true);
    destination_->setTextInteractionFlags(Qt::TextSelectableByMouse);
    layout->addWidget(destination_);
    progress_ = new QProgressBar;
    progress_->setObjectName("stepExportProgress");
    progress_->setTextVisible(false);
    layout->addWidget(progress_);
    status_ = new QLabel(tr("Bereit"));
    status_->setObjectName("stepExportStatus");
    status_->setTextFormat(Qt::PlainText);
    status_->setWordWrap(true);
    status_->setTextInteractionFlags(Qt::TextSelectableByMouse);
    layout->addWidget(status_);
    auto *buttons = new QHBoxLayout;
    buttons->addStretch();
    export_ = new QPushButton(tr("Exportieren …"));
    export_->setObjectName("stepExportStart");
    cancel_ = new QPushButton(tr("Abbrechen"));
    cancel_->setObjectName("stepExportCancel");
    cancel_->setToolTip(tr("Fordert den Abbruch an. Die aktuelle CAD-Einzeloperation darf noch "
                           "enden; die Zieldatei wird bei Abbruch nicht ersetzt."));
    cancel_->setWhatsThis(cancel_->toolTip());
    auto *close = new QPushButton(tr("Schließen"));
    close->setObjectName("stepExportClose");
    buttons->addWidget(export_);
    buttons->addWidget(cancel_);
    buttons->addWidget(close);
    layout->addLayout(buttons);
    connect(kind_, &QComboBox::currentIndexChanged, this, [this] { updateControls(); });
    connect(export_, &QPushButton::clicked, this, [this] { chooseDestination(); });
    connect(cancel_, &QPushButton::clicked, this, [this] { cancelExport(); });
    connect(close, &QPushButton::clicked, this, &StepExportDialog::reject);
    updateControls();
}

StepExportDialog::~StepExportDialog() {
    if (job_)
        job_->cancel();
}

void StepExportDialog::chooseDestination() {
    QFileDialog dialog(this, tr("Draft als STEP speichern"));
    dialog.setAcceptMode(QFileDialog::AcceptSave);
    dialog.setFileMode(QFileDialog::AnyFile);
    dialog.setNameFilter(tr("STEP-Geometrie (*.step *.stp)"));
    dialog.setDefaultSuffix("step");
    dialog.selectFile(targetPath_.isEmpty() ? tr("Horn-Draft-r%1.step").arg(revision_)
                                            : targetPath_);
    if (dialog.exec() == QDialog::Accepted && !dialog.selectedFiles().isEmpty())
        startExport(dialog.selectedFiles().front());
}

void StepExportDialog::startExport(const QString &confirmedPath) {
    if (busy())
        throw std::logic_error("STEP-Export läuft bereits.");
    if (confirmedPath.trimmed().isEmpty())
        throw std::invalid_argument("STEP-Export benötigt einen Zieldateinamen.");
    targetPath_ = QFileInfo(confirmedPath).absoluteFilePath();
    error_.clear();
    cancellationRequested_ = false;
    destination_->setText(targetPath_);
    const DraftStepExportOptions options{
        static_cast<DraftStepExportKind>(kind_->currentData().toInt()),
        tolerance_->value() / 1000.0, static_cast<std::size_t>(sections_->value())};
    try {
        job_ = new StepExportJob(this, geometry_, options);
        status_->setText(tr("STEP-Geometrie wird aufgebaut …"));
    } catch (const std::exception &error) {
        error_ = QString::fromUtf8(error.what());
        status_->setText(tr("Export konnte nicht starten: %1").arg(error_));
    }
    updateControls();
}

void StepExportDialog::cancelExport() {
    if (!job_)
        return;
    cancellationRequested_ = true;
    job_->cancel();
    status_->setText(tr("Abbruch angefordert …"));
    updateControls();
}

bool StepExportDialog::busy() const noexcept {
    return !job_.isNull();
}

QString StepExportDialog::lastError() const {
    return error_;
}

void StepExportDialog::reject() {
    cancelExport();
    QDialog::reject();
}

void StepExportDialog::updateControls() {
    const bool running = busy();
    kind_->setEnabled(!running);
    tolerance_->setEnabled(!running);
    const bool contours =
        kind_->currentData().toInt() == static_cast<int>(DraftStepExportKind::sections);
    sections_->setVisible(contours);
    sectionsLabel_->setVisible(contours);
    sections_->setEnabled(!running && contours);
    export_->setEnabled(!running);
    cancel_->setEnabled(running && !cancellationRequested_);
    progress_->setRange(0, running ? 0 : 1);
    progress_->setValue(0);
    progress_->setVisible(running);
}

void StepExportDialog::finished(DraftStepExportResult result, const QString &error,
                                bool cancelled) {
    job_.clear();
    if (cancelled || cancellationRequested_) {
        status_->setText(tr("Abgebrochen · Zieldatei unverändert"));
    } else if (!error.isEmpty()) {
        error_ = error;
        status_->setText(tr("Export fehlgeschlagen: %1").arg(error));
    } else {
        QSaveFile file(targetPath_);
        // Never fall back to truncating the destination if the sibling temporary
        // file cannot be created; cancellation/error must preserve existing data.
        file.setDirectWriteFallback(false);
        if (result.data.empty() ||
            result.data.size() > static_cast<std::size_t>(std::numeric_limits<qint64>::max())) {
            error_ = tr("STEP-Ausgabe ist leer oder zu groß.");
        } else if (!file.open(QIODevice::WriteOnly)) {
            error_ = file.errorString();
        } else if (file.write(result.data.data(), static_cast<qint64>(result.data.size())) !=
                   static_cast<qint64>(result.data.size())) {
            error_ = file.errorString();
            file.cancelWriting();
        } else if (cancellationRequested_) {
            file.cancelWriting();
        } else if (!file.commit()) {
            error_ = file.errorString();
        }
        if (!error_.isEmpty())
            status_->setText(tr("Export fehlgeschlagen: %1").arg(error_));
        else if (cancellationRequested_)
            status_->setText(tr("Abgebrochen · Zieldatei unverändert"));
        else
            status_->setText(tr("Gespeichert · Stand %1 · größte geprüfte Abweichung %2 mm")
                                 .arg(revision_)
                                 .arg(result.maximumDeviationMetres * 1000.0, 0, 'g', 5));
    }
    updateControls();
}
} // namespace hornsculpt
