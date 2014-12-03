#include "controller.h"
#include "mainwindow.h"
#include "mesh.h"
#include <string>
#include <cassert>
#include <iomanip>
#include <sstream>
#include <QMessageBox>
#include <fstream>

using namespace std;
using namespace Eigen;

Controller::Controller(MainWindow &mw) : mw_(mw), m_()
{
    m_.buildOpenSchwarzLantern(1.0, 3.0, 4, 4, 3.14159/4.0);
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

void Controller::loadSimulation()
{
    string filename = mw_.launchSimulationOpenDialog();
    ifstream ifs(filename.c_str());
    if(!ifs || !m_.loadFromStream(ifs))
        mw_.showError(string("Couldn't open simulation file ") + filename);
    mw_.centerCamera();
}

void Controller::saveSimulation()
{
    string filename = mw_.launchSimulationSaveDialog();
    ofstream ofs(filename.c_str());
    if(!ofs || !m_.saveToStream(ofs))
        mw_.showError(string("Couldn't save simulation file ") + filename);
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
    bool open = FALSE;
    bool springs = FALSE;
    mw_.launchSchwarzLanternDialog(r, h, n, m, angle,open,springs);
//    if(open){
//        m_.buildOpenSchwarzLantern(r, h, n, m, angle,springs);
//    }
//    else{
    m_.buildSchwarzLantern(r,h,n,m,angle,open,springs);
//    }
    mw_.centerCamera();
}

void Controller::deformLantern()
{
    m_.crushLantern(*this, 1e-3);
}

void Controller::repaintCallback()
{
    static int frame=1;
    stringstream fname;
    fname << "frame" << setw(6) << setfill('0') << frame++ << ".png";
    mw_.repaintMesh();
    mw_.saveScreenshot(fname.str());
}

void Controller::exportOBJ(const char *filename)
{
    if(!m_.exportOBJ(filename))
    {
        QString msg = "Couldn't write file " + QString(filename) + ". Save failed.";
        QMessageBox::warning(&mw_, "Couldn't Write File", msg, QMessageBox::Ok);
        return;
    }
}
void Controller::importOBJ(const char *filename)
{
    if(!m_.importOBJ(filename))
    {
        QString msg = "Couldn't read file " + QString(filename) + ". Import failed.";
        QMessageBox::warning(&mw_, "Couldn't read File", msg, QMessageBox::Ok);
        return;
    }
}
