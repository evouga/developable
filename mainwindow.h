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
    void launchSchwarzLanternDialog(double &r, double &h, int &n, int &m);
    void saveScreenshot();
    void saveScreenshot(const std::string &filename);
    void showError(const std::string &error);
    void centerCamera();
    bool showWireframe() const;
    bool smoothShade() const;
    void setCylinderHeight(double height);
\
private slots:
    void on_actionExit_triggered();

    void on_actionLoad_OBJ_triggered();

    void on_actionReset_Camera_triggered();

    void on_actionTake_Screenshot_triggered();

    void on_wireframeCheckBox_clicked();

    void on_smoothShadeCheckBox_clicked();

    void on_actionSchwarz_Lantern_triggered();

    void on_optimizeButton_clicked();

    void on_heightSlider_actionTriggered(int action);

private:
    void updateGL();

    Ui::MainWindow *ui;
    Controller *cont_;
};

#endif // MAINWINDOW_H
