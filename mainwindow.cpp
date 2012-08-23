#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "controller.h"
#include <QFileDialog>
#include <QMessageBox>
#include <QDir>
#include <QDateTime>
#include "schwarzdialog.h"

using namespace std;

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow),
    cont_(NULL)
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::setController(Controller &cont)
{
    cont_ = &cont;
    ui->GLwidget->setController(cont);
}

string MainWindow::launchMeshOpenDialog()
{
    string filename = string(QFileDialog::getOpenFileName(this,
        tr("Open Mesh"), "", tr("Mesh Files (*.obj *.ply)")).toAscii());

    return filename;
}

void MainWindow::launchSchwarzLanternDialog(double &r, double &h, int &n, int &m)
{
    SchwarzDialog sd(this);
    sd.setDefaultParameters(r, h, n, m);
    if(sd.exec() == QDialog::Accepted)
        sd.getChosenParameters(r,h,n,m);
}

void MainWindow::showError(const string &error)
{
    QMessageBox::warning(this, tr("VViewer"),
                                    QString(error.c_str()),
                                    QMessageBox::Ok, QMessageBox::NoButton);
}

void MainWindow::centerCamera()
{
    ui->GLwidget->centerCamera();
}


void MainWindow::saveScreenshot()
{
    QString curdir = QDir::currentPath();

    QDateTime dateTime = QDateTime::currentDateTime();
    QString dateTimeString = dateTime.toString("dd_MM_yy_hh_mm_ss_zzz");
    string filename = string(curdir.toAscii()) + "/screen_" + string(dateTimeString.toAscii()) + ".png";
    cout << "path " << filename << endl;
    ui->GLwidget->saveScreenshot(filename);
}

bool MainWindow::showWireframe() const
{
    return ui->wireframeCheckBox->isChecked();
}

bool MainWindow::smoothShade() const
{
    return ui->smoothShadeCheckBox->isChecked();
}

void MainWindow::updateGL()
{
    ui->GLwidget->updateGL();
}

void MainWindow::setCylinderHeight(double height)
{
    ui->heightSlider->setMaximum(100*height);
    ui->heightSlider->setValue(100*height);
    repaint();
}


void MainWindow::on_actionExit_triggered()
{
    assert(cont_);
    cont_->quit();
}

void MainWindow::on_actionLoad_OBJ_triggered()
{
    assert(cont_);
    cont_->loadOBJ();
}

void MainWindow::on_actionReset_Camera_triggered()
{
    centerCamera();
    updateGL();
}

void MainWindow::on_actionTake_Screenshot_triggered()
{
    saveScreenshot();
    updateGL();
}

void MainWindow::on_wireframeCheckBox_clicked()
{
    updateGL();
}

void MainWindow::on_smoothShadeCheckBox_clicked()
{
    updateGL();
}

void MainWindow::on_actionSchwarz_Lantern_triggered()
{
    cont_->newSchwarzLantern();
    updateGL();
}

void MainWindow::on_optimizeButton_clicked()
{
    cont_->deformLantern();
    updateGL();
}

void MainWindow::on_heightSlider_actionTriggered(int action)
{
    int position = ui->heightSlider->value();
    double height = position/100.0;
    cont_->updateLanternHeight(height);
    updateGL();
}
