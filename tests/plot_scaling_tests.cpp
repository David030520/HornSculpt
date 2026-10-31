#include "../modules/visualization/src/chart_scaling.hpp"
#include "../modules/visualization/src/plot_legend_layout.hpp"
#include <QApplication>
#include <QDialog>
#include <QDialogButtonBox>
#include <QImage>
#include <QLineEdit>
#include <QPainter>
#include <QPushButton>
#include <QSettings>
#include <QTemporaryDir>
#include <QTimer>
#include <cmath>
#include <hornsculpt/heatmap_widget.hpp>
#include <hornsculpt/level_palette.hpp>
#include <hornsculpt/plot_widget.hpp>
#include <hornsculpt/polar_plot_widget.hpp>
#include <hornsculpt/pressure_field_plane.hpp>
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace hornsculpt;
namespace {
void require(bool condition, const char *message) {
    if (!condition)
        throw std::runtime_error(message);
}
template <class F> void rejects(F f) {
    bool caught = false;
    try {
        f();
    } catch (const std::invalid_argument &) {
        caught = true;
    }
    require(caught, "Invalid scale accepted.");
}
QImage render(QWidget &widget) {
    widget.resize(720, 440);
    QImage image(widget.size(), QImage::Format_ARGB32);
    image.fill(Qt::transparent);
    QPainter painter(&image);
    widget.render(&painter);
    return image;
}
void validationAndTicks() {
    require(PlotDbScale{}.maximum == 3, "Default dB maximum must be +3.");
    const auto ticks = plotTicks(-30, 3, 3);
    require(ticks.front() == -30 && ticks.back() == 3 && ticks.size() == 12,
            "Default dB ticks incorrect.");
    require(plotTicks(0, 180, 30).size() == 7, "Beamwidth angle ticks incorrect.");
    validatePlotAxisScale({.001, 1000000., .125}, true);
    validatePlotDbScale({-45, 12.5, .125});
    validatePlotDbScale({0, 999, 1});
    rejects([] { validatePlotDbScale({0, 1000, 1}); });
    rejects([] { validatePlotDbScale({-30, 3, .001}); });
    rejects([] { validatePlotAxisScale({1., {}, {}}); });
    rejects([] { validatePlotAxisScale({0., 1., 0.}); });
    rejects([] { validatePlotAxisScale({0., 1., {}}, true); });
    rejects([] { validatePlotDbScale({3, -30, 3}); });
    rejects([] { validatePlotDbScale({-30, std::numeric_limits<double>::infinity(), 3}); });
    const auto dense = plotTicks(-10, 10, 1e-12, 300);
    require(!dense.empty() && dense.size() <= 1002, "Dense grid exhausted rendering resources.");
}
void pressureFieldBands() {
    require(levelPalette().front() == QColor(153, 0, 0) &&
                levelPalette()[3] == QColor(255, 211, 0) &&
                levelPalette()[6] == QColor(0, 216, 255) &&
                levelPalette().back() == QColor(0, 0, 132),
            "The original isobar palette anchors changed.");
    PressureFieldAppearance appearance;
    appearance.minimumDb = 40;
    appearance.maximumDb = 110;
    const PlotDbScale scale{40, 110, 3};
    require(!steppedLevelColor(0., {-1e308, 1e308, 1e308}).isValid(),
            "An overflowing scale span must not produce a palette colour.");
    require(!pressureFieldColor(std::numeric_limits<double>::quiet_NaN(), appearance).isValid() &&
                !steppedLevelColor(std::numeric_limits<double>::quiet_NaN(), scale).isValid() &&
                pressureFieldColor(-std::numeric_limits<double>::infinity(), appearance) ==
                    pressureFieldColor(scale.minimum, appearance) &&
                steppedLevelColor(-std::numeric_limits<double>::infinity(), scale) ==
                    steppedLevelColor(scale.minimum, scale),
            "Non-finite samples must not acquire a palette colour.");
    for (double value : {20., 40., 42., 51.5, 107.5, 110., 130.})
        require(pressureFieldColor(value, appearance) == steppedLevelColor(value, scale),
                "Field SPL bands differ from the isobar palette.");
    require(pressureFieldColor(108., appearance) == pressureFieldColor(109., appearance) &&
                pressureFieldColor(106., appearance) != pressureFieldColor(109., appearance),
            "SPL quantization does not follow 3 dB bands.");
    appearance.quantity = 1;
    const auto seam = pressureFieldColor(-180., appearance);
    require(seam == pressureFieldColor(180., appearance) &&
                seam == pressureFieldColor(-179., appearance) &&
                seam == pressureFieldColor(179., appearance) &&
                seam != pressureFieldColor(-161., appearance),
            "Phase wrap must stay inside one 36-degree band.");
    appearance.quantity = 2;
    const auto range = pressureFieldRange(appearance);
    const double delta = range[1] * .01;
    require(pressureFieldColor(0., appearance) == pressureFieldColor(-delta, appearance) &&
                pressureFieldColor(0., appearance) != pressureFieldColor(delta, appearance) &&
                pressureFieldColor(range[1] * 2, appearance) == pressureFieldColor(range[1], appearance),
            "Real-pressure bins must meet at zero and clamp at the displayed limit.");
}
void externalLegendLayout() {
    const QFontMetricsF font(QApplication::font());
    const std::vector<QString> labels{QStringLiteral("On-Axis"), QStringLiteral("Listening Window"),
        QStringLiteral("Early Reflections"), QStringLiteral("Sound Power · Näherung"),
        QStringLiteral("Raumprognose"), QString(100, QLatin1Char('W'))};
    for (const double width : {244., 510., 1000.}) {
        const auto layout = chart::legendLayout(labels, font, width);
        require(layout.items.size() == labels.size() && layout.height > 0,
                "Legend lost a labelled series.");
        for (std::size_t i = 0; i < layout.items.size(); ++i) {
            const auto &item = layout.items[i];
            require(item.seriesIndex == i && item.text.right() <= width + .001 &&
                        item.text.bottom() <= layout.height + .001 &&
                        item.swatch.bottom() <= layout.height + .001,
                    "A legend entry exceeds its reserved region.");
            for (std::size_t j = 0; j < i; ++j)
                require(!item.text.intersects(layout.items[j].text) &&
                            !item.text.intersects(layout.items[j].swatch) &&
                            !item.swatch.intersects(layout.items[j].text),
                        "Legend entries overlap.");
        }
    }
    PlotWidget narrow("CTA", "Hz", "dB");
    std::vector<PlotSeries> series;
    for (const auto &label : labels)
        series.push_back({label, Qt::blue, {{100, 0}, {1000, 1}}});
    narrow.resize(260, 155);
    narrow.setSeries(series);
    auto legendFont = narrow.font();
    legendFont.setPointSizeF(std::max(7.0, legendFont.pointSizeF() - 1));
    const auto layout = chart::legendLayout(labels, QFontMetricsF(legendFont), 244.);
    require(narrow.height() >= 156 + layout.height,
            "A narrow plot did not reserve room for its legend and data area.");
    const int narrowMinimum = narrow.minimumHeight();
    QImage image(narrow.size(), QImage::Format_ARGB32);
    image.fill(Qt::white);
    QPainter painter(&image);
    narrow.render(&painter);
    painter.end();
    require(!image.isNull(), "A narrow plot with a long legend did not render.");
    narrow.show();
    narrow.resize(1000, 500);
    QApplication::processEvents();
    require(narrow.minimumHeight() < narrowMinimum && narrow.height() >= narrow.minimumHeight(),
            "A wide plot retained the narrow legend's unnecessary minimum height.");
}
void persistenceAndDataUpdates() {
    const PlotAxisScale x{100, 20000, 1000}, y{0, 180, 30};
    {
        PlotWidget plot("Beamwidth", "Hz", "°");
        plot.setFrequencyAxis();
        plot.setObjectName("testBeamwidth");
        plot.setAxisScales(x, y);
        plot.setSeries({{"A", Qt::blue, {{100, 30}, {1000, 160}}}});
        plot.setSeries({{"B", Qt::red, {{200, 20}, {2000, 150}}}});
        require(plot.axisScales() == PlotAxisScales{x, y}, "setSeries reset fixed scales.");
        require(!render(plot).isNull(), "Fixed-axis plot did not render.");
    }
    PlotWidget restored("Beamwidth", "Hz", "°");
    restored.setFrequencyAxis();
    restored.setObjectName("testBeamwidth");
    require(restored.axisScales() == PlotAxisScales{x, y}, "Named chart scale did not restore.");
    rejects([&] { restored.setAxisScales({0, 100, 10}, y); });
    restored.setDecibelScale({-24, 3, .125});
    require(restored.decibelScale() == PlotDbScale{-24, 3, .125}, "Shared dB API changed values.");
}
void heatmapAndPolar() {
    HeatmapWidget map("H");
    map.setObjectName("testHeatmap");
    map.setAxisScales({100, 1000, 100}, {-90, 90, 30});
    map.setLevels({-90, 90}, {100, 1000}, {3, 3, 3, 3});
    const auto positive = render(map).pixelColor(250, 160);
    map.setLevels({-90, 90}, {100, 1000}, {0, 0, 0, 0});
    const auto zero = render(map).pixelColor(250, 160);
    require(positive != zero, "Heatmap clamped +3 dB to 0 dB.");
    map.setLevels({-90, 90}, {100, 1000}, {-30, -30, -30, -30});
    const auto low = render(map).pixelColor(250, 160);
    map.setLevels({-90, 90}, {100, 1000}, {-INFINITY, -INFINITY, -INFINITY, -INFINITY});
    require(render(map).pixelColor(250, 160) == low,
            "Zero-pressure -infinity was hidden instead of clipped.");
    map.setAxisScales({50, 2000, 250}, {-120, 120, 30});
    map.setLevels({-90, 90}, {100, 1000}, {1, 1, 1, 1});
    require(map.axisScales().y.minimum == -120, "Heatmap data reset angle scale.");
    map.setDecibelScale({-60, 6, 2});
    require(map.decibelScale() == PlotDbScale{-60, 6, 2}, "Heatmap shared dB scale failed.");
    require(!render(map).isNull(), "Expanded heatmap failed to render.");
    PolarPlotWidget polar("Polar", Qt::blue);
    polar.setObjectName("testPolar");
    polar.setAngularTickStep(15);
    polar.setDecibelScale({-30, 3, 3});
    polar.setLevels({{-1, 0}, {0, 3}, {1, -INFINITY}});
    require(!render(polar).isNull(), "Polar +3/-infinity failed to render.");
    require(polar.angularTickStep() == 15 && polar.decibelScale().maximum == 3,
            "Polar scales changed with data.");
    rejects([&] { polar.setAngularTickStep(0); });
}
void dialogAppliesExactValues() {
    QWidget owner;
    QTimer::singleShot(0, [] {
        auto *dialog = qobject_cast<QDialog *>(QApplication::activeModalWidget());
        require(dialog, "No scaling dialog.");
        dialog->findChild<QLineEdit *>("scaleDbMinimum")->setText(QLocale().toString(-42.5));
        dialog->findChild<QLineEdit *>("scaleDbMaximum")->setText(QLocale().toString(3.));
        dialog->findChild<QLineEdit *>("scaleDbStep")->setText(QLocale().toString(.125));
        dialog->findChild<QLineEdit *>("scaleAngularStep")->setText(QLocale().toString(12.5));
        dialog->findChild<QDialogButtonBox *>()->button(QDialogButtonBox::Ok)->click();
    });
    const auto result = editPolarScales(owner, {{}, 30});
    require(result && result->decibels == PlotDbScale{-42.5, 3, .125} &&
                result->angularStep == 12.5,
            "Scaling dialog rounded or discarded user values.");
}
void radialHeatmap() {
    HeatmapWidget map("Isobaren rund H · Testdaten");
    map.setRadialLayout(true);
    map.setLevels({-90, 0, 90}, {100, 1000, 10000}, {0, 0, 0, -12, -12, -12, -30, -30, -30});
    const auto image = render(map);
    if (const auto path = qEnvironmentVariable("HORNSCULPT_RADIAL_PREVIEW"); !path.isEmpty())
        require(image.save(path), "Could not save radial preview.");
    // At 720x440 the disk center is (325,205), radius 160. Avoid grid lines.
    require(image.pixelColor(345, 155) != image.pixelColor(370, 75),
            "Radial frequency bands did not change with distance.");
    require(image.pixelColor(340, 260) == image.pixelColor(170, 50),
            "Uncomputed rear angles were filled.");
    const auto radialLow = image.pixelColor(370, 75);
    map.setLevels({-90, 0, 90}, {100, 1000, 10000}, {0, 0, 0, 0, 0, 0, 0, 0, 0});
    require(render(map).pixelColor(370, 75) != radialLow,
            "Radial raster did not refresh after result change.");
    map.setLevels({}, {}, {});
    const auto cleared = render(map);
    require(cleared.pixelColor(370, 75) == cleared.pixelColor(10, 75),
            "Clearing results retained radial colors.");
}
} // namespace
int main(int argc, char **argv) {
    QApplication app(argc, argv);
    QTemporaryDir settings;
    QCoreApplication::setOrganizationName("HornSculptTests");
    QCoreApplication::setApplicationName("PlotScaling");
    QSettings::setDefaultFormat(QSettings::IniFormat);
    QSettings::setPath(QSettings::IniFormat, QSettings::UserScope, settings.path());
    try {
        validationAndTicks();
        pressureFieldBands();
        externalLegendLayout();
        persistenceAndDataUpdates();
        heatmapAndPolar();
        radialHeatmap();
        dialogAppliesExactValues();
        std::cout << "Plot scaling: validation, fixed axes, persistence, dialog, +3 dB and "
                     "null-pressure rendering PASS.\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
