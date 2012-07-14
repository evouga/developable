#include <OpenMesh/Core/IO/MeshIO.hh>
#include "mesh.h"
#include <GL/glu.h>
#include <Eigen/Dense>
#include "meshcurvature.h"

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

void Mesh::render(MeshCurvature &mc, bool showWireframe, bool smoothShade, HeatMap type, double cutoff)
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

    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);
    glEnableClientState(GL_COLOR_ARRAY);

    static vector<GLfloat> colors;
    static vector<int> indices;
    static vector<GLfloat> pos;
    static vector<GLfloat> normal;

    colors.clear();
    indices.clear();
    pos.clear();
    normal.clear();

    for(int i=0; i<(int)mesh_.n_vertices(); i++)
    {
        OMMesh::VertexHandle v = mesh_.vertex_handle(i);
        Vector3d color(0.0, 186/255., 0.0);
        if(type == HM_GAUSSIAN)
        {
            double curvature = mc.gaussianCurvature(v.idx());
            color = heatmap(curvature, cutoff);
        }
        else if(type == HM_MEAN)
        {
            double curvature = mc.meanCurvature(v.idx());
            color = heatmap(curvature, cutoff);
        }
        else if(type == HM_SPREAD)
        {
            double curvature = mc.curvatureSpread(v.idx());
            color = heatmap(curvature, cutoff);
        }

        OMMesh::Point pt = mesh_.point(v);
        OMMesh::Point n;
        mesh_.calc_vertex_normal_correct(v, n);
        n.normalize();
        for(int j=0; j<3; j++)
        {
            pos.push_back(pt[j]);
            normal.push_back(n[j]);
            colors.push_back(color[j]);
        }
    }

    glVertexPointer(3, GL_FLOAT, 0, &pos[0]);
    glNormalPointer(GL_FLOAT, 0, &normal[0]);
    glColorPointer(3, GL_FLOAT, 0, &colors[0]);

    OMMesh::ConstFaceIter f, fEnd = mesh_.faces_end();
    int i=0;
    for (f = mesh_.faces_begin(); f != fEnd; ++f,i++) {
        for (OMMesh::ConstFaceVertexIter v = mesh_.cfv_iter(f); v; ++v) {
            indices.push_back(v.handle().idx());
        }
    }
    glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, &indices[0]);

    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_NORMAL_ARRAY);
    glDisableClientState(GL_COLOR_ARRAY);
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

Vector3d Mesh::vertexNormal(int vidx) const
{
    assert(0 <= vidx && vidx < (int)mesh_.n_vertices());
    OMMesh::VertexHandle vh = mesh_.vertex_handle(vidx);
    OMMesh::Normal normal;
    mesh_.calc_vertex_normal_correct(vh, normal);
    Vector3d result(normal[0], normal[1], normal[2]);
    result.normalize();
    return result;
}

double Mesh::shortestAdjacentEdge(int vidx) const
{
    assert(0 <= vidx && vidx < (int)mesh_.n_vertices());
    OMMesh::VertexHandle vh = mesh_.vertex_handle(vidx);
    double mindist = std::numeric_limits<double>::infinity();
    for(OMMesh::ConstVertexEdgeIter vei = mesh_.cve_iter(vh); vei; ++vei)
    {
        double dist = mesh_.calc_edge_sqr_length(vei.handle());
        if(dist < mindist)
            mindist = dist;
    }
    return sqrt(mindist);
}

Vector3d Mesh::heatmap(double val, double max)
{
    if(val < 0)
    {
        val = 1.0 - std::min(fabs(val), max)/max;
        return Vector3d(val, val, 1.0);
    }
    else
    {
        val = 1.0 - std::min(val, max)/max;
        return Vector3d(1.0, val, val);
    }
}

double Mesh::areaOfInfluence(int vidx) const
{
    assert(0 <= vidx && vidx < (int)mesh_.n_vertices());
    double result = 0;
    OMMesh::VertexHandle vh = mesh_.vertex_handle(vidx);
    for(OMMesh::ConstVertexFaceIter vfi = mesh_.cvf_iter(vh); vfi; ++vfi)
    {
        double area = faceArea(vfi.handle().idx());
        result += 1.0/3.0 * area;
    }
    return result;
}

double Mesh::faceArea(int fidx) const
{
    assert(0 <= fidx && fidx < (int)mesh_.n_faces());
    OMMesh::FaceHandle fh = mesh_.face_handle(fidx);
    Vector3d points[3];
    int i=0;
    for(OMMesh::ConstFaceVertexIter fvi = mesh_.cfv_iter(fh); fvi; ++fvi,++i)
    {
        assert(i < 3);
        OMMesh::Point pt = mesh_.point(fvi.handle());
        for(int j=0; j<3; j++)
            points[i][j] = pt[j];
    }
    Vector3d e1 = points[1]-points[0];
    Vector3d e2 = points[2]-points[0];
    return 0.5 * e1.cross(e2).norm();
}
