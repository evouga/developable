#include "controller.h"
#include "mainwindow.h"
#include "mesh.h"
#include <string>
#include <cassert>

using namespace std;
using namespace Eigen;

Controller::Controller(MainWindow &mw) : mw_(mw), m_()
{
    m_.buildSchwarzLantern(1.0, 1.0, 4, 4);
}

void Controller::quit()
{
    mw_.close();
}

void Controller::renderMesh()
{
    bool showWireframe = mw_.showWireframe();
    bool smoothShade = mw_.smoothShade();
    m_.render(showWireframe, smoothShade);
}

void Controller::loadOBJ()
{
    string filename = mw_.launchMeshOpenDialog();
    if(!m_.loadMesh(filename))
        mw_.showError(string("Couldn't open mesh file ") + filename);
    mw_.centerCamera();

}

void Controller::getSceneBounds(Eigen::Vector3d &center, double &radius)
{
    center = m_.centroid();
    radius = m_.radius();
}

void Controller::newSchwarzLantern()
{
    double r = 1.0;
    double h = 1.0;
    int n = 4;
    int m = 4;
    mw_.launchSchwarzLanternDialog(r, h, n, m);
    m_.buildSchwarzLantern(r, h, n, m);
}
