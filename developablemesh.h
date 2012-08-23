#ifndef DEVELOPABLEMESH_H
#define DEVELOPABLEMESH_H

#include <vector>
#include <string>
#include "mesh.h"
#include <Eigen/Sparse>
#include <fadiff.h>

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

    void buildObjective(double &f, Eigen::VectorXd &Df, Eigen::SparseMatrix<double> &Hf);
    void buildConstraints(Eigen::VectorXd &g, Eigen::SparseMatrix<double> &Dg, const std::vector<double> &targetHeights);
    void projectOntoConstraints(const Eigen::VectorXd &v, const Eigen::SparseMatrix<double> &Dg, Eigen::VectorXd &result);

    fadbad::F<double> norm(fadbad::F<double> *v);
    void cross(fadbad::F<double> *v1, fadbad::F<double> *v2, fadbad::F<double> *result);
    void normalize(fadbad::F<double> *v);
    fadbad::F<double> dot(fadbad::F<double> *v1, fadbad::F<double> *v2);
};

#endif // DEVELOPABLEMESH_H
