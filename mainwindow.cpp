#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "controller.h"
#include <QFileDialog>
#include <QMessageBox>
#include <QDir>
#include <QDateTime>

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

bool MainWindow::showWireframe()
{
    return ui->wireframeCheckBox->isChecked();
}

bool MainWindow::smoothShade()
{
    return ui->smoothShadeCheckBox->isChecked();
}

Mesh::HeatMap MainWindow::getHeatMapType()
{
    if(ui->noneCurvatureButton->isChecked())
        return Mesh::HM_NONE;
    else if(ui->meanCurvatureButton->isChecked())
        return Mesh::HM_MEAN;
    else if(ui->gaussianCurvatureButton->isChecked())
        return Mesh::HM_GAUSSIAN;
    return Mesh::HM_NONE;
}

void MainWindow::updateGL()
{
    ui->GLwidget->updateGL();
}

double MainWindow::curvatureCutoff()
{
    return ui->cutoffSlider->value()*0.001;
}

bool MainWindow::showRulings()
{
    return ui->rulingsBox->isChecked();
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

void MainWindow::on_gaussianCurvatureButton_clicked()
{
    updateGL();
}

void MainWindow::on_meanCurvatureButton_clicked()
{
    updateGL();
}

void MainWindow::on_noneCurvatureButton_clicked()
{
    updateGL();
}

void MainWindow::on_cutoffSlider_actionTriggered(int )
{
    ui->curvatureCutoff->setText(QString::number(ui->cutoffSlider->value()*0.001));
    updateGL();
}

void MainWindow::on_rulingsBox_clicked()
{
    updateGL();
}
