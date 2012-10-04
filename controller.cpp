#include "controller.h"
#include "mainwindow.h"
#include "mesh.h"
#include <string>
#include <cassert>
#include <iomanip>
#include <sstream>

using namespace std;
using namespace Eigen;

Controller::Controller(MainWindow &mw) : mw_(mw), m_()
{
    m_.buildSchwarzLantern(1.0, 3.0, 4, 4, 3.14159/4);
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

void Controller::renderMaterial()
{
    m_.renderMaterial();
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

void Controller::getMaterialBounds(Eigen::Vector2d &center, double &radius)
{
    center = m_.materialCenter();
    radius = m_.materialRadius();
}

void Controller::newSchwarzLantern()
{
 //   double r = 1.0;
    double r = 0.53033;
    double h = 2.0;
    int n = 4;
    int m = 3;
    //double angle = 3.14159/n;
    double angle = -6.28319;
    mw_.launchSchwarzLanternDialog(r, h, n, m, angle);
    m_.buildSchwarzLantern(r, h, n, m, angle);
    mw_.centerCamera();
}

void Controller::deformLantern()
{
    for(int i=0; i<1; i++)
    {
        m_.deformLantern(50);
        stringstream ss;
        ss << "frame_" << setfill('0') << setw(6) << i << ".png";
        mw_.saveScreenshot(ss.str());
    }
}
