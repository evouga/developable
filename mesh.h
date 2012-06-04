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
    ~Mesh();

    bool loadMesh(const std::string &filename);
    OMMesh &getMesh() {return mesh_;}
    void render(bool showWireframe, bool smoothShade);

    Eigen::Vector3d centroid();
    double radius();

private:
    OMMesh mesh_;
    GLUquadric *quadric_;

    void edgeEndpoints(OMMesh::EdgeHandle eh, OMMesh::Point &pt1, OMMesh::Point &pt2);
    void drawSphere(int vertex);
};

#endif // MESH_H
