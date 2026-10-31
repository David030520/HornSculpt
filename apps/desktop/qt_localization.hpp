#pragma once

class QTranslator;

namespace hornsculpt {
// The translator must outlive the QApplication's widgets.
bool installGermanQtTranslations(QTranslator &translator);
}
