#include "qt_localization.hpp"

#include <QCoreApplication>
#include <QTranslator>

namespace hornsculpt {
bool installGermanQtTranslations(QTranslator &translator) {
    return translator.load(QStringLiteral(":/hornsculpt/i18n/qtbase_de.qm")) &&
           QCoreApplication::installTranslator(&translator);
}
}
