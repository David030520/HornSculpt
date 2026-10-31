#pragma once

#include <QDoubleSpinBox>
#include <QSpinBox>

namespace hornsculpt {

class DecimalInput final : public QDoubleSpinBox {
  public:
    explicit DecimalInput(QWidget *parent = nullptr);

  protected:
    QValidator::State validate(QString &text, int &position) const override;
    double valueFromText(const QString &text) const override;

  private:
    QString normalized(QString text) const;
};

class IntegerInput final : public QSpinBox {
  public:
    explicit IntegerInput(QWidget *parent = nullptr);
};

} // namespace hornsculpt
