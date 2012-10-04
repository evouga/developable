#ifndef DEVELOPABLEMESH_H
#define DEVELOPABLEMESH_H

#include <vector>
#include <string>
#include <map>
#include "mesh.h"
#include <Eigen/Sparse>
#include <fadiff.h>
#include "autodiffTemplates.h"

typedef Eigen::Triplet<double> T;

class MaterialMesh : public Mesh
{
public:
    MaterialMesh(double W, double H) : W_(W), H_(H) {}
    double getW() {return W_;}
    double getH() {return H_;}

    std::vector<std::pair<int, int> > &getIdentifiedVerts() {return identifiedverts_;}
    std::vector<std::pair<int, double> > &getBoundaryVerts() {return bdryverts_;}
    int materialEdge(int embeddedEdge);

    void setMaterialEdge(int embeddedEdge, int materialEdge);

private:
    double W_;
    double H_;
    std::vector<std::pair<int, int> > identifiedverts_;
    std::vector<std::pair<int, double> > bdryverts_;
    std::map<int, int> embedge2matedge_;
};

struct Boundary
{
    std::vector<int> bdryVerts;
    std::vector<Eigen::Vector3d> bdryPos;
};

class DevelopableMesh : public Mesh
{
public:
    DevelopableMesh();
    ~DevelopableMesh();

    void buildSchwarzLantern(double r, double h, int n, int m, double angle);
    virtual bool loadMesh(const std::string &filename);

    void deformLantern(int maxiters);
    void buildObjective(const Eigen::VectorXd &q, double &f, Eigen::VectorXd &Df, std::vector<T> &Hf);
    void buildConstraints(const Eigen::VectorXd &q, Eigen::VectorXd &g, std::vector<T> &Dg, std::vector<std::vector<T> > &Hg);
    void buildInversionConstraints(const Eigen::VectorXd &q, Eigen::VectorXd &h, std::vector<T> &Dh, std::vector<std::vector<T> > &Hh);

    Eigen::Vector2d materialCenter();
    double materialRadius();
    void renderMaterial();

private:
    DevelopableMesh(const DevelopableMesh &other);
    DevelopableMesh &operator=(const DevelopableMesh &other);
    std::vector<Boundary> boundaries_;
    MaterialMesh *material_;

    Eigen::Vector3d point2Vector(OMMesh::Point pt);

    void centerCylinder();
    void checkInversion(const Eigen::VectorXd &q);

};

#endif // DEVELOPABLEMESH_H
