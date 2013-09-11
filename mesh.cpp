#include <OpenMesh/Core/IO/MeshIO.hh>
#include "mesh.h"
#include <GL/glu.h>
#include <Eigen/Dense>
#include <vector>
\
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

bool Mesh::saveToStream(std::ostream &os)
{
    int nverts = mesh_.n_vertices();
    int nfaces = mesh_.n_faces();
    writeInt(os, nverts);
    for(int i=0; i<nverts; i++)
    {
        OMMesh::Point pt = mesh_.point(mesh_.vertex_handle(i));
        writeDouble(os, pt[0]);
        writeDouble(os, pt[1]);
        writeDouble(os, pt[2]);

        OMMesh::Point vel = mesh_.data(mesh_.vertex_handle(i)).vel();
        writeDouble(os, vel[0]);
        writeDouble(os, vel[1]);
        writeDouble(os, vel[2]);
    }
    writeInt(os, nfaces);
    for(int i=0; i<nfaces; i++)
    {
        OMMesh::FaceHandle fh = mesh_.face_handle(i);
        vector<int> nfverts;
        for(OMMesh::FaceVertexIter fvi = mesh_.fv_iter(fh); fvi; ++fvi)
            nfverts.push_back(fvi.handle().idx());
        writeInt(os, (int)nfverts.size());
        for(int j=0; j<(int)nfverts.size(); j++)
            writeInt(os, nfverts[j]);
    }
    return os;
}

bool Mesh::loadFromStream(std::istream &is)
{
    mesh_ = OMMesh();
    int nverts = readInt(is);
    if(!is)
    {
        assert(false);
        return false;
    }
    for(int i=0; i<nverts; i++)
    {
        OMMesh::Point newpt;
        newpt[0] = readDouble(is);
        newpt[1] = readDouble(is);
        newpt[2] = readDouble(is);
        if(!is)
        {
            assert(false);
            return false;
        }

        OMMesh::Point newvel;
        newvel[0] = readDouble(is);
        newvel[1] = readDouble(is);
        newvel[2] = readDouble(is);
        if(!is)
        {
            assert(false);
            return false;
        }

        OMMesh::VertexHandle newvh = mesh_.add_vertex(newpt);
        mesh_.data(newvh).set_vel(newvel);
    }
    int nfaces = readInt(is);
    if(!is)
    {
        assert(false);
        return false;
    }
    for(int i=0; i<nfaces; i++)
    {
        vector<OMMesh::VertexHandle> faceverts;
        int nfverts = readInt(is);
        if(!is)
        {
            assert(false);
            return false;
        }
        for(int j=0; j<nfverts; j++)
        {
            int vidx = readInt(is);
            if(!is)
            {
                assert(false);
                return false;
            }
            faceverts.push_back(mesh_.vertex_handle(vidx));
        }
        mesh_.add_face(faceverts);
    }
    return true;
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

    for(OMMesh::FaceIter fi = mesh_.faces_begin(); fi != mesh_.faces_end(); ++fi)
    {
        double straind = mesh_.data(fi.handle()).strainDensity();
        for(OMMesh::FaceVertexIter fvi = mesh_.fv_iter(fi.handle()); fvi; ++fvi)
        {
            Vector3d color;
            OMMesh::VertexHandle v = fvi.handle();
            color = Vector3d(straind, (1.0 - 0.5*straind)*186/255., 0.0);

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
    }

    glVertexPointer(3, GL_FLOAT, 0, &pos[0]);
    glNormalPointer(GL_FLOAT, 0, &normal[0]);
    glColorPointer(3, GL_FLOAT, 0, &colors[0]);

    int idx=0;
    for (int i=0; i<(int)mesh_.n_faces(); i++)
    {
        for(OMMesh::FaceVertexIter fvi = mesh_.fv_iter(mesh_.face_handle(i)); fvi; ++fvi)
        {
            indices.push_back(idx++);
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

int Mesh::findEdge(int vid1, int vid2)
{
    OMMesh::VertexHandle vhstart = mesh_.vertex_handle(vid1);
    for(OMMesh::VertexOHalfedgeIter voh = mesh_.voh_iter(vhstart); voh; ++voh)
    {
        int tov = mesh_.to_vertex_handle(voh.handle()).idx();
        if(tov == vid2)
            return mesh_.edge_handle(voh.handle()).idx();
    }
    return -1;
}

int Mesh::findHalfedge(int vid1, int vid2)
{
    OMMesh::VertexHandle vhstart = mesh_.vertex_handle(vid1);
    for(OMMesh::VertexOHalfedgeIter voh = mesh_.voh_iter(vhstart); voh; ++voh)
    {
        int tov = mesh_.to_vertex_handle(voh.handle()).idx();
        if(tov == vid2)
            return voh.handle().idx();
    }
    return -1;
}

double Mesh::edgeLength(int vid1, int vid2)
{
    int eidx = findEdge(vid1, vid2);
    return mesh_.calc_edge_length(mesh_.edge_handle(eidx));
}

double Mesh::edgeLength(int eid)
{
    OMMesh::EdgeHandle eh = mesh_.edge_handle(eid);
    OMMesh::HalfedgeHandle heh = mesh_.halfedge_handle(eh, 0);
    return edgeLength(mesh_.from_vertex_handle(heh).idx(), mesh_.to_vertex_handle(heh).idx());
}

bool Mesh::exportOBJ(const char *filename)
{
    OpenMesh::IO::Options opt;
    mesh_.request_face_normals();
    mesh_.request_vertex_normals();
    mesh_.update_normals();
    opt.set(OpenMesh::IO::Options::VertexNormal);
    return OpenMesh::IO::write_mesh(mesh_, filename, opt);
}

void Mesh::writeInt(std::ostream &os, int i)
{
    os.write((const char *)&i, sizeof(int));
}

void Mesh::writeDouble(std::ostream &os, double d)
{
    os.write((const char *)&d, sizeof(double));
}

void Mesh::writeBool(std::ostream &os, bool b)
{
    writeInt(os, b ? 1 : 0);
}

int Mesh::readInt(std::istream &is)
{
    int result;
    is.read((char *)&result, sizeof(int));
    return result;
}

double Mesh::readDouble(std::istream &is)
{
    double result;
    is.read((char *)&result, sizeof(double));
    return result;
}

bool Mesh::readBool(std::istream &is)
{
    return readInt(is) != 0;
}
