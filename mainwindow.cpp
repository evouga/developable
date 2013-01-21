#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "controller.h"
#include <QFileDialog>
#include <QMessageBox>
#include <QDir>
#include <QDateTime>
#include "schwarzdialog.h"
#include <QPixmap>

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
    ui->GLwidget2D->setController(cont);
}

string MainWindow::launchMeshOpenDialog()
{
    string filename = string(QFileDialog::getOpenFileName(this,
        tr("Open Mesh"), "", tr("Mesh Files (*.obj *.ply)")).toAscii());

    return filename;
}

void MainWindow::launchSchwarzLanternDialog(double &r, double &h, int &n, int &m, double &angle)
{
    SchwarzDialog sd(this);
    sd.setDefaultParameters(r, h, n, m, angle);
    if(sd.exec() == QDialog::Accepted)
        sd.getChosenParameters(r,h,n,m, angle);
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
    ui->GLwidget2D->centerCamera();
}


void MainWindow::saveScreenshot()
{
    QDateTime dateTime = QDateTime::currentDateTime();
    QString dateTimeString = dateTime.toString("dd_MM_yy_hh_mm_ss_zzz");
    string filename = "screen_" + string(dateTimeString.toAscii()) + ".png";
    saveScreenshot(filename);
}

void MainWindow::saveScreenshot(const string &filename)
{
    updateGL();
    QPixmap p = QPixmap::grabWidget(this);
    QString curdir = QDir::currentPath();
    string fullname = string(curdir.toAscii()) + "/output/" + filename;
    p.save(QString::fromUtf8(fullname.c_str()));
    //ui->GLwidget->saveScreenshot(fullname);
}

bool MainWindow::showWireframe() const
{
    return ui->wireframeCheckBox->isChecked();
}

bool MainWindow::smoothShade() const
{
    return ui->smoothShadeCheckBox->isChecked();
}

void MainWindow::repaintMesh()
{
    updateGL();
}

void MainWindow::updateGL()
{
    ui->GLwidget->updateGL();
    ui->GLwidget2D->updateGL();
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

void MainWindow::on_actionExport_OBJ_triggered()
{
    QFileDialog savedialog(this, "Export 3D Geometry", ".", "Mesh Files (*.obj)");
    savedialog.setFileMode(QFileDialog::AnyFile);
    savedialog.setDefaultSuffix("obj");
    savedialog.setViewMode(QFileDialog::List);
    savedialog.setAcceptMode(QFileDialog::AcceptSave);
    if(savedialog.exec())
    {
        QStringList filenames = savedialog.selectedFiles();
        if(filenames.size() > 0)
        {
            QString filename = filenames[0];
            cont_->exportOBJ(filename.toStdString().c_str());
        }
    }
}
