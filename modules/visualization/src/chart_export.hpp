#pragma once
#include <functional>
class QWidget;
namespace hornsculpt {
void enableChartExport(QWidget &chart, std::function<void()> scaling);
}
