#ifndef SCHWARZDIALOG_H
#define SCHWARZDIALOG_H

#include <QDialog>

namespace Ui {
    class SchwarzDialog;
}

class SchwarzDialog : public QDialog
{
    Q_OBJECT

public:
    explicit SchwarzDialog(QWidget *parent = 0);
    ~SchwarzDialog();

    void setDefaultParameters(double r, double h, int n, int m);
    void getChosenParameters(double &r, double &h, int &n, int &m);

private:
    Ui::SchwarzDialog *ui;
};

#endif // SCHWARZDIALOG_H
