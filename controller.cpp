#include "controller.h"
#include "mainwindow.h"
#include "mesh.h"
#include <string>
#include <cassert>

using namespace std;
using namespace Eigen;

Controller::Controller(MainWindow &mw) : mw_(mw), m_(), mc_(m_), mcontours_()
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
    
    if(mw_.showContours())
    {
        mcontours_.renderContours();
    }
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

void Controller::setNumContours( int num )
{
    if( num < 2 ) num = 2;
    
    std::vector< double > isovals;
    isovals.reserve( num );
    const double cutoff = mw_.curvatureCutoff();
    for( int i = 0; i < num; ++i )
    {
        const double isoval = -cutoff + ( 2*cutoff )*double(i)/(num-1);
        isovals.push_back( isoval );
    }
    
    std::vector< double > vertex_vals;
    vertex_vals.reserve( m_.getMesh().n_vertices() );
    const Mesh::HeatMap type = mw_.getHeatMapType();
    if( Mesh::HM_NONE == type )
    {
        mcontours_.clearPrecomputedContours();
        return;
    }
    
    for( unsigned int i = 0; i < m_.getMesh().n_vertices(); i++ )
    {
        switch( type )
        {
            case Mesh::HM_MEAN:
                vertex_vals.push_back( mc_.meanCurvature( i ) );
                break;
            
            case Mesh::HM_GAUSSIAN:
                vertex_vals.push_back( mc_.gaussianCurvature( i ) );
                break;
            
            default:
                assert(!"Shouldn't be here");
                vertex_vals.push_back( -31337 );
        }
    }
    assert( vertex_vals.size() == m_.getMesh().n_vertices() );
    
    mcontours_.precomputeContoursAtValuesForMeshUsingVertexValues( isovals, m_.getMesh(), vertex_vals );
}
