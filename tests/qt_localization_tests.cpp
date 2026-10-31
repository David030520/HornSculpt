#include "qt_localization.hpp"

#include <QAction>
#include <QApplication>
#include <QDialog>
#include <QDialogButtonBox>
#include <QFileDialog>
#include <QLabel>
#include <QLineEdit>
#include <QLocale>
#include <QMenu>
#include <QMessageBox>
#include <QPushButton>
#include <QTranslator>
#include <algorithm>
#include <cstdio>
#include <memory>
#include <stdexcept>

namespace {
void require(bool condition, const char *message) {
    if (!condition)
        throw std::runtime_error(message);
}
QString clean(QString value) {
    return value.remove('&').section('\t', 0, 0);
}
} // namespace

int main(int argc, char **argv) {
    QApplication app(argc, argv);
    QLocale::setDefault(QLocale(QLocale::German, QLocale::Germany));
    try {
        QTranslator translator;
        require(hornsculpt::installGermanQtTranslations(translator),
                "Bundled Qt German translation did not load");
        QMessageBox message(QMessageBox::Question, QStringLiteral("Projekt speichern"),
                            QStringLiteral("Änderungen speichern?"),
                            QMessageBox::Save | QMessageBox::Discard | QMessageBox::Cancel);
        message.show();
        app.processEvents();
        require(clean(message.button(QMessageBox::Save)->text()) == QStringLiteral("Speichern") &&
                    clean(message.button(QMessageBox::Discard)->text()) == QStringLiteral("Verwerfen") &&
                    clean(message.button(QMessageBox::Cancel)->text()) == QStringLiteral("Abbrechen"),
                "Standard save confirmation buttons are not German");
        QDialog dialog;
        QDialogButtonBox buttons(QDialogButtonBox::Save | QDialogButtonBox::Cancel |
                                     QDialogButtonBox::Help, &dialog);
        dialog.show();
        buttons.show();
        app.processEvents();
        require(clean(buttons.button(QDialogButtonBox::Save)->text()) == QStringLiteral("Speichern") &&
                    clean(buttons.button(QDialogButtonBox::Cancel)->text()) == QStringLiteral("Abbrechen") &&
                    clean(buttons.button(QDialogButtonBox::Help)->text()) == QStringLiteral("Hilfe"),
                "Dialog button box labels are not German");
        QFileDialog file;
        file.setOption(QFileDialog::DontUseNativeDialog, true);
        file.setAcceptMode(QFileDialog::AcceptOpen);
        file.show();
        app.processEvents();
        const auto labels = file.findChildren<QLabel *>();
        const auto labelPresent = [&](const QString &expected) {
            return std::any_of(labels.begin(), labels.end(),
                               [&](const QLabel *label) { return clean(label->text()) == expected; });
        };
        require(labelPresent(QStringLiteral("Suchen in:")) &&
                    labelPresent(QStringLiteral("Dateiname:")) &&
                    labelPresent(QStringLiteral("Dateien des Typs:")),
                "Qt file dialog fallback contains untranslated labels");
        QLineEdit input;
        input.setText(QStringLiteral("test"));
        const auto menu = std::unique_ptr<QMenu>(input.createStandardContextMenu());
        QStringList actions;
        for (const auto *action : menu->actions())
            actions.push_back(clean(action->text()));
        for (const auto &expected : {QStringLiteral("Rückgängig"), QStringLiteral("Wiederherstellen"),
                                     QStringLiteral("Ausschneiden"), QStringLiteral("Kopieren"),
                                     QStringLiteral("Einfügen"), QStringLiteral("Löschen"),
                                     QStringLiteral("Alles auswählen")})
            require(actions.contains(expected), "Text context menu contains an untranslated action");
        std::puts("Bundled Qt translations: message, dialog, file fallback and text menu passed.");
        return 0;
    } catch (const std::exception &error) {
        std::fprintf(stderr, "%s\n", error.what());
        return 1;
    }
}
