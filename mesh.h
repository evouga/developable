#ifndef MESH_H
#define MESH_H

#include <string>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <Eigen/Core>

typedef OpenMesh::TriMesh_ArrayKernelT<OpenMesh::DefaultTraits>  OMMesh;

class GLUquadric;
class MeshCurvature;

class Mesh
{
public:
    enum HeatMap {HM_NONE, HM_MEAN, HM_GAUSSIAN};

    Mesh();
    ~Mesh();

    bool loadMesh(const std::string &filename);
    OMMesh &getMesh() {return mesh_;}
    const OMMesh &getMesh() const {return mesh_;}
    void render(MeshCurvature &mc, bool showWireframe, bool smoothShade, HeatMap type, double cutoff);

    Eigen::Vector3d centroid();
    double radius();

    Eigen::Vector3d vertexNormal(int vidx) const;
    double shortestAdjacentEdge(int vidx) const;

    double areaOfInfluence(int vidx) const;
    double faceArea(int fidx) const;

private:
    OMMesh mesh_;
    GLUquadric *quadric_;

    void edgeEndpoints(OMMesh::EdgeHandle eh, OMMesh::Point &pt1, OMMesh::Point &pt2);
    void drawSphere(int vertex);
    Eigen::Vector3d heatmap(double val, double max);
};

#endif // MESH_H
