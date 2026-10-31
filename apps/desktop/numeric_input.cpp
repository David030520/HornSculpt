#include "numeric_input.hpp"

#include <QApplication>
#include <QEvent>
#include <QLineEdit>
#include <QLocale>
#include <QMouseEvent>
#include <QTimer>

namespace hornsculpt {
namespace {
class SelectOnEntry final : public QObject {
  public:
    explicit SelectOnEntry(QLineEdit *editor) : QObject(editor), editor_(editor) {
        editor_->installEventFilter(this);
        editor_->parentWidget()->installEventFilter(this);
    }

  protected:
    bool eventFilter(QObject *watched, QEvent *event) override {
        if (watched == editor_ || watched == editor_->parentWidget()) {
            switch (event->type()) {
            case QEvent::FocusIn:
                firstClick_ = true;
                break;
            case QEvent::FocusOut:
                armed_ = false;
                firstClick_ = true;
                break;
            case QEvent::MouseButtonPress:
                if (watched == editor_) {
                    armed_ = firstClick_ &&
                             static_cast<QMouseEvent *>(event)->button() == Qt::LeftButton;
                    press_ = static_cast<QMouseEvent *>(event)->position();
                }
                break;
            case QEvent::MouseButtonRelease:
                if (watched != editor_)
                    break;
                if (armed_ &&
                    (static_cast<QMouseEvent *>(event)->position() - press_).manhattanLength() <
                        QApplication::startDragDistance()) {
                    // Select after QLineEdit has positioned its cursor for the click.
                    QTimer::singleShot(0, editor_, [editor = editor_] {
                        if (editor->hasFocus())
                            editor->selectAll();
                    });
                }
                armed_ = false;
                firstClick_ = false;
                break;
            case QEvent::KeyPress:
                firstClick_ = false;
                break;
            default:
                break;
            }
        }
        return QObject::eventFilter(watched, event);
    }

  private:
    QLineEdit *editor_;
    QPointF press_;
    bool firstClick_{true};
    bool armed_{};
};
} // namespace

DecimalInput::DecimalInput(QWidget *parent) : QDoubleSpinBox(parent) {
    setKeyboardTracking(false);
    setGroupSeparatorShown(false);
    new SelectOnEntry(lineEdit());
}

QString DecimalInput::normalized(QString text) const {
    // Both decimal keys mean a decimal separator; grouping is not used in these fields.
    const auto separator = locale().decimalPoint();
    text.replace(QStringLiteral("."), separator);
    text.replace(QStringLiteral(","), separator);
    return text;
}

QValidator::State DecimalInput::validate(QString &text, int &position) const {
    text = normalized(std::move(text));
    return QDoubleSpinBox::validate(text, position);
}

double DecimalInput::valueFromText(const QString &text) const {
    return QDoubleSpinBox::valueFromText(normalized(text));
}

IntegerInput::IntegerInput(QWidget *parent) : QSpinBox(parent) {
    setKeyboardTracking(false);
    setGroupSeparatorShown(false);
    new SelectOnEntry(lineEdit());
}

} // namespace hornsculpt
