#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <string>

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
    bool showWireframe();
    bool smoothShade();
    
private slots:
    void on_actionExit_triggered();

    void on_actionLoad_OBJ_triggered();

    void on_actionReset_Camera_triggered();

    void on_actionTake_Screenshot_triggered();

    void on_wireframeCheckBox_clicked();

    void on_smoothShadeCheckBox_clicked();

private:
    Ui::MainWindow *ui;
    Controller *cont_;
};

#endif // MAINWINDOW_H
