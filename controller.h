#ifndef CONTROLLER_H
#define CONTROLLER_H

#include "developablemesh.h"
#include <Eigen/Core>

class MainWindow;

class Controller
{
public:
    Controller(MainWindow &mw);

    void renderMesh();
    void quit();
    void loadOBJ();
    void getSceneBounds(Eigen::Vector3d &center, double &radius);
    void newSchwarzLantern();

private:
    MainWindow &mw_;
    DevelopableMesh m_;
};

#endif // CONTROLLER_H
