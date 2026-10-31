#include "branding.hpp"

#include <QAction>
#include <QApplication>
#include <QCryptographicHash>
#include <QDialog>
#include <QEvent>
#include <QFile>
#include <QIcon>
#include <QLabel>
#include <QMainWindow>
#include <QPixmap>
#include <QPushButton>
#include <QTextBrowser>
#include <QUrl>
#include <cstdio>
#include <stdexcept>

namespace {
void require(bool condition, const char *message) {
    if (!condition)
        throw std::runtime_error(message);
}
} // namespace

int main(int argc, char **argv) {
    QApplication application(argc, argv);
    application.setApplicationVersion(QStringLiteral(HORNSCULPT_VERSION));
    try {
        require(argc <= 2, "Use branding_tests [screenshot.png]");
        QFile original(QStringLiteral(":/hornsculpt/branding/logo.png"));
        require(original.open(QIODevice::ReadOnly), "Original logo resource is missing");
        require(QCryptographicHash::hash(original.readAll(), QCryptographicHash::Sha256).toHex() ==
                    "988fd96e69c8496fd8f80035039617259d2da47b2e358bdcf3227cf6a5acf9f7",
                "Embedded logo differs from the supplied original");
        QMainWindow window;
        hornsculpt::installBranding(window);
        window.show();
        require(!window.windowIcon().isNull() && !QApplication::windowIcon().isNull(),
                "Application or window icon is missing");
        auto *action = window.findChild<QAction *>(QStringLiteral("aboutHornSculpt"));
        require(action != nullptr, "About menu action is missing");
        action->trigger();
        application.processEvents();
        auto *dialog = window.findChild<QDialog *>(QStringLiteral("aboutHornSculptDialog"));
        require(dialog && dialog->isVisible(), "About dialog did not open");
        const auto *image = dialog->findChild<QLabel *>(QStringLiteral("brandLogo"));
        const auto *version = dialog->findChild<QLabel *>(QStringLiteral("brandVersion"));
        require(image && !image->pixmap().isNull() && version &&
                    version->text().contains(QCoreApplication::applicationVersion()),
                "About dialog does not display the real logo and application version");
        action->trigger();
        require(window.findChildren<QDialog *>(QStringLiteral("aboutHornSculptDialog")).size() == 1,
                "Repeated About action created duplicate dialogs");
        if (argc == 2)
            require(dialog->grab().save(application.arguments()[1]), "About screenshot failed");
        dialog->close();
        QCoreApplication::sendPostedEvents(nullptr, QEvent::DeferredDelete);
        auto *guide = window.findChild<QAction *>(QStringLiteral("helpGuide"));
        auto *knowledge = window.findChild<QAction *>(QStringLiteral("helpKnowledge"));
        require(guide && knowledge, "Offline help actions are missing");
        guide->trigger();
        application.processEvents();
        auto *help = window.findChild<QDialog *>(QStringLiteral("helpDialog"));
        auto *browser = help ? help->findChild<QTextBrowser *>(QStringLiteral("helpBrowser")) : nullptr;
        auto *back = help ? help->findChild<QPushButton *>(QStringLiteral("helpBack")) : nullptr;
        require(help && help->isVisible() && browser && back &&
                    browser->toPlainText().contains(QStringLiteral("Kurzanleitung")),
                "German quick guide did not open from embedded resources");
        QMetaObject::invokeMethod(browser, "anchorClicked", Qt::DirectConnection,
                                  Q_ARG(QUrl, QUrl(QStringLiteral("FreeCAD-LINK.md"))));
        require(browser->source().path().endsWith(QStringLiteral("/FreeCAD-LINK.md")) &&
                    browser->toPlainText().contains(QStringLiteral("FreeCAD")) && back->isEnabled(),
                "Relative guide link or navigation history failed");
        back->click();
        require(browser->source().path().endsWith(QStringLiteral("/CUSTOMER-GUIDE.md")),
                "Back navigation did not restore the guide");
        knowledge->trigger();
        require(browser->source().path().endsWith(QStringLiteral("/knowledge/README.md")) &&
                    browser->toPlainText().contains(QStringLiteral("Wissen")),
                "Knowledge index is not available offline");
        QMetaObject::invokeMethod(browser, "anchorClicked", Qt::DirectConnection,
                                  Q_ARG(QUrl, QUrl(QStringLiteral("validation.md"))));
        QMetaObject::invokeMethod(browser, "anchorClicked", Qt::DirectConnection,
                                  Q_ARG(QUrl, QUrl(QStringLiteral("../CTA-2034.md"))));
        require(browser->source().path().endsWith(QStringLiteral("/CTA-2034.md")) &&
                    browser->toPlainText().contains(QStringLiteral("CTA")),
                "Knowledge links do not resolve across resource directories");
        help->close();
        QCoreApplication::sendPostedEvents(nullptr, QEvent::DeferredDelete);
        std::puts("Original branding resource, icons and native About dialog passed.");
        return 0;
    } catch (const std::exception &error) {
        std::fprintf(stderr, "%s\n", error.what());
        return 1;
    }
}
