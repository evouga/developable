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
    ui->GLwidget->updateGL();
}

void MainWindow::on_actionTake_Screenshot_triggered()
{
    saveScreenshot();
    ui->GLwidget->updateGL();
}

void MainWindow::on_wireframeCheckBox_clicked()
{
    ui->GLwidget->updateGL();
}

void MainWindow::on_smoothShadeCheckBox_clicked()
{
    ui->GLwidget->updateGL();
}
