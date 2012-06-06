#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <string>
#include "mesh.h"

class Controller;

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT
    
public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

    void setController(Controller &cont);

    std::string launchMeshOpenDialog();
    void saveScreenshot();
    void showError(const std::string &error);
    void centerCamera();
    const bool showWireframe();
    const bool smoothShade();
    const Mesh::HeatMap getHeatMapType();
    const double curvatureCutoff();
    const bool showRulings();
    const bool showContours();
    
private slots:
    void on_actionExit_triggered();

    void on_actionLoad_OBJ_triggered();

    void on_actionReset_Camera_triggered();

    void on_actionTake_Screenshot_triggered();

    void on_wireframeCheckBox_clicked();

    void on_smoothShadeCheckBox_clicked();

    void on_gaussianCurvatureButton_clicked();

    void on_meanCurvatureButton_clicked();

    void on_noneCurvatureButton_clicked();

    void on_cutoffSlider_actionTriggered(int action);

    void on_rulingsBox_clicked();

    void on_contoursBox_clicked();

    void on_contoursSlider_valueChanged(int value);

private:
    void updateGL();

    Ui::MainWindow *ui;
    Controller *cont_;
};

#endif // MAINWINDOW_H
