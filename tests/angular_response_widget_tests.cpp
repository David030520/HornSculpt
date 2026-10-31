#include "angular_response_widget.hpp"
#include <QApplication>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QLabel>
#include <QLineEdit>
#include <QSettings>
#include <QTemporaryDir>
#include <cmath>
#include <cstdio>
#include <hornsculpt/plot_widget.hpp>
#include <limits>
#include <stdexcept>

using namespace hornsculpt;
namespace {
void require(bool value, const char *message) {
    if (!value)
        throw std::runtime_error(message);
}
void close(double actual, double expected) {
    require(std::isfinite(actual) && std::abs(actual - expected) < 1e-10,
            "Angular frequency level differs from independent fixture.");
}
} // namespace
int main(int argc, char **argv) {
    QApplication app(argc, argv);
    QCoreApplication::setOrganizationName("HornSculptTests");
    QCoreApplication::setApplicationName("AngularResponses");
    QTemporaryDir directory;
    QSettings::setDefaultFormat(QSettings::IniFormat);
    QSettings::setPath(QSettings::IniFormat, QSettings::UserScope, directory.path());
    try {
        AngularResponseWidget widget;
        auto *h = dynamic_cast<PlotWidget *>(widget.findChild<QWidget *>("angularFrequencyH"));
        auto *v = dynamic_cast<PlotWidget *>(widget.findChild<QWidget *>("angularFrequencyV"));
        auto *angles = widget.findChild<QLineEdit *>("angularResponseAngles");
        auto *mode = widget.findChild<QComboBox *>("angularResponseNormalization");
        auto *reference = widget.findChild<QDoubleSpinBox *>("angularResponseReferenceAngle");
        auto *error = widget.findChild<QLabel *>("angularResponseError");
        require(h && v && angles && mode && reference && error, "Missing native angular controls.");
        require(widget.settings() == AngularResponseSettings{} && !reference->isEnabled(),
                "Absolute display or default zero-degree reference changed.");
        const auto zero = -std::numeric_limits<double>::infinity();
        widget.setRows({{100, {{0, 80}, {30, 74}, {60, 70}}, {{0, 80}, {30, 70}, {60, 60}}},
                        {400, {{0, 90}, {30, 84}, {60, 80}}, {{0, 90}, {30, 78}, {60, 65}}},
                        {1000, {{0, 100}, {30, zero}, {60, 90}}, {{0, 100}, {30, 95}, {60, 85}}}},
                       "Rev. 7");
        require(h->series().size() == 5 && h->series()[0].points.size() == 3,
                "Default curves or actual row count missing.");
        close(h->series()[0].points[1].x(), 400);
        close(h->series()[0].points[1].y(), 90);
        close(v->series()[2].points[1].y(), 78);
        require(std::isnan(h->series()[2].points[2].y()),
                "A true pressure zero became a fabricated level or removed frequency row.");
        int changes = 0;
        widget.onChanged = [&] { ++changes; };
        mode->setCurrentIndex(1);
        require(reference->isEnabled() && changes == 1,
                "Reference mode did not enable its editable zero-degree reference.");
        close(h->series()[0].points[0].y(), 0);
        close(h->series()[2].points[0].y(), -6);
        close(v->series()[2].points[1].y(), -12);
        reference->setValue(15);
        close(h->series()[0].points[0].y(), 3);
        close(v->series()[0].points[0].y(), 5);
        require(std::isnan(h->series()[0].points[2].y()) &&
                    std::isnan(h->series()[4].points[2].y()),
                "Missing reference pressure did not break every curve at that frequency.");
        angles->setText("22,5; -15");
        require(widget.settings().anglesDegrees == std::vector<double>({22.5, -15}),
                "Decimal comma or signed angular selection failed.");
        close(h->series()[0].points[0].y(), -1.5);
        close(v->series()[0].points[0].y(), -2.5);
        require(std::isnan(h->series()[1].points[0].y()), "An absent angle was extrapolated.");
        const auto valid = widget.settings();
        const int validChanges = changes;
        for (const auto &invalid : {"0;", "0; 0", "181", "nan", "0,1,2"}) {
            angles->setText(invalid);
            require(widget.settings() == valid && changes == validChanges && !error->isHidden(),
                    "Invalid angular input replaced the last valid selection or lacked feedback.");
        }
        AngularResponseSettings restored{{-30, 0, 22.5}, true, -22.5};
        widget.setSettings(restored);
        require(widget.settings() == restored && changes == validChanges && error->isHidden(),
                "Restoring display settings emitted an edit or lost its reference.");
        widget.setRows({{100, {}, {}}, {400, {{0, 90}}, {{0, 85}}}, {1000, {}, {}}});
        require(h->series()[1].points.size() == 3 && std::isnan(h->series()[1].points[0].y()),
                "Partial rows were removed, allowing a curve to bridge a gap.");
        widget.setSettings({{0}, false, 0});
        close(h->series()[0].points[1].y(), 90);
        close(v->series()[0].points[1].y(), 85);
        widget.resize(1200, 450);
        widget.show();
        app.processEvents();
        require(!widget.grab().isNull(), "Native angular charts did not render.");
        widget.setRows({});
        require(h->series().empty() && v->series().empty(),
                "Cleared result retained stale traces.");
        PlotWidget isolated("", "Hz", "dB");
        isolated.resize(600, 300);
        isolated.setFrequencyAxis();
        isolated.setAxisScales({100, 1000, {}}, {0, 100, 20});
        isolated.setSeries(
            {{"",
              Qt::blue,
              {{200, 20}, {400, std::numeric_limits<double>::quiet_NaN()}, {800, 80}}}});
        const auto pixels = isolated.grab().toImage();
        for (const auto &point : {QPointF(200, 20), QPointF(800, 80)}) {
            const int x = qRound(65 + (std::log10(point.x()) - 2) * (isolated.width() - 81));
            const int y = qRound(28 + (1 - point.y() / 100) * (isolated.height() - 67));
            bool visible = false;
            for (int dx = -3; dx <= 3; ++dx)
                for (int dy = -3; dy <= 3; ++dy) {
                    const auto color = pixels.pixelColor(x + dx, y + dy);
                    visible |= color.blue() > 150 && color.red() < 100 && color.green() < 100;
                }
            require(visible, "An isolated computed point next to a data gap disappeared.");
        }
        std::puts(
            "H/V angular responses, reference normalization, gaps, input and restore passed.");
        return 0;
    } catch (const std::exception &error) {
        std::fprintf(stderr, "%s\n", error.what());
        return 1;
    }
}
