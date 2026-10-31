#include "branding.hpp"

#include <QAction>
#include <QApplication>
#include <QDialog>
#include <QDialogButtonBox>
#include <QDesktopServices>
#include <QHBoxLayout>
#include <QIcon>
#include <QLabel>
#include <QMainWindow>
#include <QMenuBar>
#include <QPainter>
#include <QPixmap>
#include <QPushButton>
#include <QScreen>
#include <QTextBrowser>
#include <QUrl>
#include <QVBoxLayout>
#include <stdexcept>

namespace hornsculpt {
namespace {
QPixmap logo() {
    const QPixmap original(QStringLiteral(":/hornsculpt/branding/logo.png"));
    if (original.isNull())
        throw std::runtime_error("Eingebettetes HornSculpt-Logo fehlt.");
    // The original asset includes a large black canvas. Crop only the display,
    // retaining a 20-pixel black margin around the complete, unchanged wordmark.
    const auto cropped = original.copy(QRect(91, 725, 1855, 550));
    QPixmap rendered(cropped.size());
    rendered.fill(Qt::black);
    {
        QPainter painter(&rendered);
        painter.drawPixmap(0, 0, cropped);
    }
    return rendered;
}
void showAbout(QMainWindow &window) {
    if (auto *existing = window.findChild<QDialog *>(QStringLiteral("aboutHornSculptDialog"))) {
        existing->raise();
        existing->activateWindow();
        return;
    }
    auto *dialog = new QDialog(&window);
    dialog->setObjectName(QStringLiteral("aboutHornSculptDialog"));
    dialog->setWindowTitle(QObject::tr("Über HornSculpt"));
    dialog->setAttribute(Qt::WA_DeleteOnClose);
    auto *layout = new QVBoxLayout(dialog);
    layout->setContentsMargins(16, 16, 16, 12);
    layout->setSpacing(12);
    auto *image = new QLabel;
    image->setObjectName(QStringLiteral("brandLogo"));
    image->setAlignment(Qt::AlignCenter);
    image->setAccessibleName(QObject::tr("HornSculpt-Logo"));
    const double ratio = window.devicePixelRatioF();
    auto pixmap = logo().scaled(QSize(qRound(440 * ratio), qRound(131 * ratio)),
                                Qt::KeepAspectRatio, Qt::SmoothTransformation);
    pixmap.setDevicePixelRatio(ratio);
    image->setPixmap(pixmap);
    layout->addWidget(image);
    auto *version =
        new QLabel(QObject::tr("Version %1").arg(QCoreApplication::applicationVersion()));
    version->setObjectName(QStringLiteral("brandVersion"));
    version->setTextInteractionFlags(Qt::TextSelectableByMouse);
    layout->addWidget(version);
    auto *buttons = new QDialogButtonBox(QDialogButtonBox::Close);
    buttons->button(QDialogButtonBox::Close)->setText(QObject::tr("Schließen"));
    QObject::connect(buttons, &QDialogButtonBox::rejected, dialog, &QDialog::close);
    layout->addWidget(buttons);
    dialog->open();
}
void showHelp(QMainWindow &window, const QString &page) {
    auto *dialog = window.findChild<QDialog *>(QStringLiteral("helpDialog"));
    QTextBrowser *browser = nullptr;
    if (!dialog) {
        dialog = new QDialog(&window);
        dialog->setObjectName(QStringLiteral("helpDialog"));
        dialog->setWindowTitle(QObject::tr("HornSculpt-Hilfe"));
        dialog->setAttribute(Qt::WA_DeleteOnClose);
        if (const auto *screen = window.screen()) {
            const auto available = screen->availableGeometry().size();
            dialog->resize(qMin(850, qMax(1, available.width() - 40)),
                           qMin(650, qMax(1, available.height() - 40)));
        } else {
            dialog->resize(850, 650);
        }
        auto *layout = new QVBoxLayout(dialog);
        auto *navigation = new QHBoxLayout;
        auto *back = new QPushButton(QObject::tr("Zurück"));
        back->setObjectName(QStringLiteral("helpBack"));
        back->setEnabled(false);
        navigation->addWidget(back);
        browser = new QTextBrowser;
        browser->setObjectName(QStringLiteral("helpBrowser"));
        browser->setOpenLinks(false);
        browser->setOpenExternalLinks(false);
        QObject::connect(back, &QPushButton::clicked, browser, &QTextBrowser::backward);
        QObject::connect(browser, &QTextBrowser::backwardAvailable, back, &QPushButton::setEnabled);
        const auto addPage = [&](const QString &label, const QString &path) {
            auto *button = new QPushButton(label);
            navigation->addWidget(button);
            QObject::connect(button, &QPushButton::clicked, browser,
                             [browser, path] { browser->setSource(QUrl(path)); });
        };
        addPage(QObject::tr("Kurzanleitung"), QStringLiteral("qrc:/hornsculpt/docs/CUSTOMER-GUIDE.md"));
        addPage(QObject::tr("Wissensbasis"), QStringLiteral("qrc:/hornsculpt/docs/knowledge/README.md"));
        addPage(QObject::tr("FreeCAD"), QStringLiteral("qrc:/hornsculpt/docs/FreeCAD-LINK.md"));
        addPage(QObject::tr("CTA-2034"), QStringLiteral("qrc:/hornsculpt/docs/CTA-2034.md"));
        addPage(QObject::tr("Schallfeld"), QStringLiteral("qrc:/hornsculpt/docs/SCHALLFELD-PHASE.md"));
        navigation->addStretch();
        layout->addLayout(navigation);
        layout->addWidget(browser);
        QObject::connect(browser, &QTextBrowser::anchorClicked, browser, [browser](const QUrl &link) {
            const auto target = browser->source().resolved(link);
            if (target.scheme() == QStringLiteral("qrc") &&
                target.path().startsWith(QStringLiteral("/hornsculpt/")))
                browser->setSource(target);
            else if (target.scheme() == QStringLiteral("http") || target.scheme() == QStringLiteral("https"))
                QDesktopServices::openUrl(target);
        });
    } else {
        browser = dialog->findChild<QTextBrowser *>(QStringLiteral("helpBrowser"));
    }
    browser->setSource(QUrl(page));
    dialog->show();
    dialog->raise();
    dialog->activateWindow();
}
} // namespace

void installBranding(QMainWindow &window) {
    const QIcon icon(logo());
    QApplication::setWindowIcon(icon);
    window.setWindowIcon(icon);
    auto *help = window.menuBar()->addMenu(QObject::tr("&Hilfe"));
    auto *guide = help->addAction(QObject::tr("Kurzanleitung …"));
    guide->setObjectName(QStringLiteral("helpGuide"));
    QObject::connect(guide, &QAction::triggered, &window, [&window] {
        showHelp(window, QStringLiteral("qrc:/hornsculpt/docs/CUSTOMER-GUIDE.md"));
    });
    auto *knowledge = help->addAction(QObject::tr("Wissensbasis …"));
    knowledge->setObjectName(QStringLiteral("helpKnowledge"));
    QObject::connect(knowledge, &QAction::triggered, &window, [&window] {
        showHelp(window, QStringLiteral("qrc:/hornsculpt/docs/knowledge/README.md"));
    });
    help->addSeparator();
    auto *about = help->addAction(QObject::tr("Über HornSculpt …"));
    about->setObjectName(QStringLiteral("aboutHornSculpt"));
    QObject::connect(about, &QAction::triggered, &window, [&window] { showAbout(window); });
}
} // namespace hornsculpt
