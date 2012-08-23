#ifndef DEVELOPABLEMESH_H
#define DEVELOPABLEMESH_H

#include <vector>
#include <string>
#include "mesh.h"

struct BoundaryCurve
{
    std::vector<int> edges;
    double arclength;
    double height;
};

class DevelopableMesh : public Mesh
{
public:
    DevelopableMesh();

    void buildSchwarzLantern(double r, double h, int n, int m);
    virtual bool loadMesh(const std::string &filename);

    void getBoundaryHeights(std::vector<double> &heights);

    void deformLantern(const std::vector<double> &newheights);

private:
    std::vector<BoundaryCurve> boundaries_;
    double surfacearea_;

    void identifyBoundaries();
    void calculateSurfaceArea();
    Eigen::Vector3d point2Vector(OMMesh::Point pt);
};

#endif // DEVELOPABLEMESH_H
