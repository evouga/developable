#include <OpenMesh/Core/IO/MeshIO.hh>
#include "mesh.h"
#include <GL/glu.h>

using namespace std;
using namespace Eigen;

Mesh::Mesh()
{
    quadric_ = gluNewQuadric();
}

Mesh::~Mesh()
{
    gluDeleteQuadric(quadric_);
}

bool Mesh::loadMesh(const string &filename)
{
    return OpenMesh::IO::read_mesh(mesh_, filename);
}

void Mesh::render(bool showWireframe, bool smoothShade)
{
    glEnable(GL_LIGHTING);
    glEnable(GL_DITHER);

    glPolygonOffset(1.0, 1.0);
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    if(smoothShade)
    {
        glShadeModel(GL_SMOOTH);
    }
    else
    {
        glShadeModel(GL_FLAT);
    }

    OMMesh::ConstFaceIter f, fEnd = mesh_.faces_end();
    int i=0;
    for (f = mesh_.faces_begin(); f != fEnd; ++f,i++) {

        glBegin(GL_POLYGON);
        for (OMMesh::ConstFaceVertexIter v = mesh_.cfv_iter(f); v; ++v) {
            glColor3f(156/255., 186/255., 214/255.);

            OMMesh::Point pt = mesh_.point(v);
            OMMesh::Point n;
            mesh_.calc_vertex_normal_correct(v, n);
            n.normalize();
            glNormal3d(n[0], -fabs(n[1]), n[2]);
            double offset = 0.0;
            glVertex3d(pt[0]+offset*n[0],pt[1]+offset*n[1],pt[2]+offset*n[2]);

        }
        glEnd();
    }

    glDisable(GL_POLYGON_OFFSET_FILL);

    if(showWireframe)
    {
        glLineWidth(1.0);
        glBegin(GL_LINES);
        for(OMMesh::ConstEdgeIter ei = mesh_.edges_begin(); ei != mesh_.edges_end(); ++ei)
        {
            glColor3f(0.0, 0.0, 0.0);
            OMMesh::Point pt1, pt2;
            edgeEndpoints(ei.handle(), pt1, pt2);
            glVertex3d(pt1[0], pt1[1], pt1[2]);
            glVertex3d(pt2[0], pt2[1], pt2[2]);
        }
        glEnd();
    }

    /*glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    for(OMMesh::ConstVertexIter v = mesh_.vertices_begin(); v != mesh_.vertices_end(); ++v)
    {
        glColor3f(0,0,0);
        drawSphere(v.handle().idx());
    }*/
}

void Mesh::edgeEndpoints(OMMesh::EdgeHandle eh, OMMesh::Point &pt1, OMMesh::Point &pt2)
{
    OMMesh::HalfedgeHandle heh1 = mesh_.halfedge_handle(eh, 0);
    pt1 = mesh_.point(mesh_.from_vertex_handle(heh1));
    pt2 = mesh_.point(mesh_.to_vertex_handle(heh1));
}

void Mesh::drawSphere(int vertex)
{
    OMMesh::VertexHandle v = mesh_.vertex_handle(vertex);
    OMMesh::Point pt = mesh_.point(v);
    double radius = std::numeric_limits<double>::infinity();
    for(OMMesh::ConstVertexEdgeIter vei = mesh_.cve_iter(v); vei; ++vei)
    {
        double len = 0.3*mesh_.calc_edge_length(vei.handle());
        if(len < radius)
            radius = len;
    }
    glPushMatrix();
    glTranslatef(pt[0],pt[1],pt[2]);
    gluSphere(quadric_, radius, 10, 10);
    glPopMatrix();
}

Vector3d Mesh::centroid()
{
    Vector3d centroid(0,0,0);
    int numpts = mesh_.n_vertices();

    for(int i=0; i<numpts; i++)
    {
        OMMesh::Point pt = mesh_.point(mesh_.vertex_handle(i));
        for(int j=0; j<3; j++)
            centroid[j] += pt[j];
    }
    centroid /= numpts;
    return centroid;
}

double Mesh::radius()
{
    Vector3d cent = centroid();
    int numpts = mesh_.n_vertices();
    double maxradius = 0;
    for(int i=0; i<numpts; i++)
    {
        OMMesh::Point pt = mesh_.point(mesh_.vertex_handle(i));
        Vector3d ept(pt[0],pt[1],pt[2]);
        double radius = (ept-cent).squaredNorm();
        if(radius > maxradius)
            maxradius = radius;
    }
    return sqrt(maxradius);
}
