#ifndef MESH_H
#define MESH_H

#include <string>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <Eigen/Core>

typedef OpenMesh::TriMesh_ArrayKernelT<OpenMesh::DefaultTraits>  OMMesh;

class GLUquadric;

class Mesh
{
public:
    Mesh();
    virtual ~Mesh();

    bool exportOBJ(const char *filename);

    virtual bool loadMesh(const std::string &filename);
    OMMesh &getMesh() {return mesh_;}
    const OMMesh &getMesh() const {return mesh_;}
    void render(bool showWireframe, bool smoothShade);

    Eigen::Vector3d centroid();
    double radius();

    Eigen::Vector3d vertexNormal(int vidx) const;
    double shortestAdjacentEdge(int vidx) const;

    double areaOfInfluence(int vidx) const;
    double faceArea(int fidx) const;

    int findEdge(int vid1, int vid2);
    int findHalfedge(int vid1, int vid2);
    virtual double edgeLength(int vid1, int vid2);
    double edgeLength(int eid);


protected:
    OMMesh mesh_;
    void edgeEndpoints(OMMesh::EdgeHandle eh, OMMesh::Point &pt1, OMMesh::Point &pt2);


private:
    Mesh(const Mesh &other);
    Mesh &operator=(const Mesh &other);
    GLUquadric *quadric_;

    void drawSphere(int vertex);
};

#endif // MESH_H
