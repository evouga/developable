#include "controller.h"
#include "mainwindow.h"
#include "mesh.h"
#include <string>

using namespace std;
using namespace Eigen;

Controller::Controller(MainWindow &mw) : mw_(mw), m_(), mc_(m_)
{
}

void Controller::quit()
{
    mw_.close();
}

void Controller::renderMesh()
{
    bool showWireframe = mw_.showWireframe();
    bool smoothShade = mw_.smoothShade();
    Mesh::HeatMap type = mw_.getHeatMapType();
    double cutoff = mw_.curvatureCutoff();
    m_.render(mc_,showWireframe, smoothShade, type, cutoff);
    if(mw_.showRulings())
        mc_.renderCurvatureDirs(m_);
}

void Controller::loadOBJ()
{
    string filename = mw_.launchMeshOpenDialog();
    if(!m_.loadMesh(filename))
        mw_.showError(string("Couldn't open mesh file ") + filename);
    mc_ = MeshCurvature(m_);
    mw_.centerCamera();

}

void Controller::getSceneBounds(Eigen::Vector3d &center, double &radius)
{
    center = m_.centroid();
    radius = m_.radius();
}
