#ifndef CONTROLLER_H
#define CONTROLLER_H

#include "mesh.h"
#include <Eigen/Core>
#include "meshcurvature.h"

class MainWindow;

class Controller
{
public:
    Controller(MainWindow &mw);

    void renderMesh();
    void quit();
    void loadOBJ();
    void getSceneBounds(Eigen::Vector3d &center, double &radius);

private:
    MainWindow &mw_;
    Mesh m_;
    MeshCurvature mc_;
};

#endif // CONTROLLER_H
