#ifndef CONTROLLER_H
#define CONTROLLER_H

#include "developablemesh.h"
#include <Eigen/Core>

class MainWindow;

class Controller : public DeformCallback
{
public:
    Controller(MainWindow &mw);

    void renderMesh();
    void renderMaterial();
    void quit();
    void loadSimulation();
    void saveSimulation();
    void getSceneBounds(Eigen::Vector3d &center, double &radius);
    void getMaterialBounds(Eigen::Vector2d &center, double &radius);
    void newSchwarzLantern();
    void deformLantern();
    void updateLanternHeight(double newheight);
    void exportOBJ(const char *filename);
    void jitterMesh();

    // Callback during solve

    virtual void repaintCallback();

private:
    MainWindow &mw_;
    DevelopableMesh m_;
};

#endif // CONTROLLER_H
