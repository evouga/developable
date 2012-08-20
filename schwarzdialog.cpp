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

void SchwarzDialog::setDefaultParameters(double r, double h, int n, int m)
{
    ui->radiusLineEdit->setText(QString::number(r));
    ui->heightLineEdit->setText(QString::number(h));
    ui->nLineEdit->setText(QString::number(n));
    ui->mLineEdit->setText(QString::number(m));
}

void SchwarzDialog::getChosenParameters(double &r, double &h, int &n, int &m)
{
    r = ui->radiusLineEdit->text().toDouble();
    h = ui->heightLineEdit->text().toDouble();
    n = ui->nLineEdit->text().toInt();
    m = ui->mLineEdit->text().toInt();
}
