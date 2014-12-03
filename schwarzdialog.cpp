#include "schwarzdialog.h"
#include "ui_schwarzdialog.h"

SchwarzDialog::SchwarzDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::SchwarzDialog)
{
    ui->setupUi(this);
}

SchwarzDialog::~SchwarzDialog()
{
    delete ui;
}

void SchwarzDialog::setDefaultParameters(double r, double h, int n, int m, double angle)
{
    ui->radiusLineEdit->setText(QString::number(r));
    ui->heightLineEdit->setText(QString::number(h));
    ui->nLineEdit->setText(QString::number(n));
    ui->mLineEdit->setText(QString::number(m));
    ui->angleLineEdit->setText(QString::number(angle));
}

void SchwarzDialog::getChosenParameters(double &r, double &h, int &n, int &m, double &angle, bool& open, bool& springs)
{
    r = ui->radiusLineEdit->text().toDouble();
    h = ui->heightLineEdit->text().toDouble();
    n = ui->nLineEdit->text().toInt();
    m = ui->mLineEdit->text().toInt();
    angle = ui->angleLineEdit->text().toDouble();
    open = ui->checkbox_openness->isChecked();
    springs = ui->checkbox_springs->isChecked();
}
