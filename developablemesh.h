#ifndef DEVELOPABLEMESH_H
#define DEVELOPABLEMESH_H

#include <vector>
#include <string>
#include "mesh.h"

struct BoundaryCurve
{
    std::vector<int> edges;
    double arclength;
};

class DevelopableMesh : public Mesh
{
public:
    DevelopableMesh();

    void buildSchwarzLantern(double r, double h, int n, int m);
    virtual bool loadMesh(const std::string &filename);

private:
    std::vector<BoundaryCurve> boundaries_;

    void identifyBoundaries();
};

#endif // DEVELOPABLEMESH_H
