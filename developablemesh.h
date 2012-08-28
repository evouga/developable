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

    void deformLantern(const std::vector<double> &newheights, int maxiters);

private:
    std::vector<BoundaryCurve> boundaries_;
    double surfacearea_;

    void identifyBoundaries();
    void calculateSurfaceArea();
    Eigen::Vector3d point2Vector(OMMesh::Point pt);

    void buildObjective(const Eigen::VectorXd &q, double &f, Eigen::VectorXd &Df, Eigen::SparseMatrix<double> &Hf);
    void buildConstraints(const Eigen::VectorXd &q, Eigen::VectorXd &g, Eigen::SparseMatrix<double> &Dg, const std::vector<double> &targetHeights);
    void projectOntoConstraints(const Eigen::VectorXd &v, const Eigen::SparseMatrix<double> &Dg, Eigen::VectorXd &result);
    void projectPositionsOntoConstraints(const Eigen::VectorXd &q, Eigen::VectorXd &result, const std::vector<double> &targetHeights);

    double lineSearch(const Eigen::VectorXd &q, const Eigen::VectorXd &dq, const Eigen::VectorXd &lambda, const std::vector<double> &targetHeights);

    template<class T> static T norm(T *v);
    template<class T> static void cross(T *v1, T *v2, T *result);
    template<class T> static void normalize(T *v);
    template<class T> static T dot(T *v1, T *v2);
};




template<class T>
void DevelopableMesh::cross(T *v1, T *v2, T *result)
{
    result[0] = v1[1]*v2[2] - v1[2]*v2[1];
    result[1] = v1[2]*v2[0] - v1[0]*v2[2];
    result[2] = v1[0]*v2[1] - v1[1]*v2[0];
}

template<class T> T DevelopableMesh::norm(T *v)
{
    T result = 0;
    for(int i=0; i<3; i++)
        result += v[i]*v[i];
    return sqrt(result);
}

template<class T> void DevelopableMesh::normalize(T *v)
{
    T n = norm(v);
    for(int i=0; i<3; i++)
        v[i] /= n;
}

template<class T> T DevelopableMesh::dot(T *v1, T *v2)
{
    T result = 0;
    for(int i=0; i<3; i++)
        result += v1[i]*v2[i];
    return result;
}


#endif // DEVELOPABLEMESH_H
