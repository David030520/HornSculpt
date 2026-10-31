#include "chart_export.hpp"
#include <QContextMenuEvent>
#include <QFileDialog>
#include <QImage>
#include <QMenu>
#include <QMessageBox>
#include <QPainter>
#include <QWidget>

namespace hornsculpt {
namespace {
class ChartExport final : public QObject {
  public:
    ChartExport(QWidget &chart, std::function<void()> scaling)
        : QObject(&chart), chart_(chart), scaling_(std::move(scaling)) {
        chart.setFocusPolicy(Qt::StrongFocus);
        chart.installEventFilter(this);
    }

  protected:
    bool eventFilter(QObject *watched, QEvent *event) override {
        if (watched != &chart_ || event->type() != QEvent::ContextMenu)
            return false;
        QMenu menu(&chart_);
        auto *scale = menu.addAction(tr("Skalierung …"));
        scale->setEnabled(static_cast<bool>(scaling_));
        auto *save = menu.addAction(tr("PNG exportieren …"));
        const auto *context = static_cast<QContextMenuEvent *>(event);
        const auto *chosen = menu.exec(context->globalPos());
        if (chosen == scale) {
            if (scaling_)
                scaling_();
            return true;
        }
        if (chosen != save)
            return true;
        const auto path = QFileDialog::getSaveFileName(
            &chart_, tr("Diagramm exportieren"), "HornSculpt-Diagramm.png", tr("PNG (*.png)"));
        if (path.isEmpty())
            return true;
        // Render vectors/text at twice the widget resolution for a usable presentation image.
        QImage image(chart_.size() * 2, QImage::Format_ARGB32_Premultiplied);
        image.setDevicePixelRatio(2);
        image.fill(Qt::white);
        QPainter painter(&image);
        chart_.render(&painter);
        painter.end();
        if (!image.save(path, "PNG"))
            QMessageBox::warning(&chart_, tr("PNG-Export"),
                                 tr("Die Bilddatei konnte nicht gespeichert werden."));
        return true;
    }

  private:
    QWidget &chart_;
    std::function<void()> scaling_;
};
} // namespace
void enableChartExport(QWidget &chart, std::function<void()> scaling) {
    new ChartExport(chart, std::move(scaling));
}
} // namespace hornsculpt
