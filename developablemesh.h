#ifndef DEVELOPABLEMESH_H
#define DEVELOPABLEMESH_H

#include <vector>
#include <string>
#include <map>
#include "mesh.h"
#include <Eigen/Sparse>
#include <fadiff.h>
#include "autodiffTemplates.h"
#include "periodicmesh.h"

typedef Eigen::Triplet<double> T;

class DeformCallback
{
public:
    virtual void repaintCallback()=0;
};

class MaterialMesh : public PeriodicMesh
{
public:
    MaterialMesh(double H) : H_(H) {}
    double getH() {return H_;}

    std::vector<std::pair<int, double> > &getBoundaryVerts() {return bdryverts_;}
    int materialEdge(int embeddedEdge);

    void setMaterialEdge(int embeddedEdge, int materialEdge);

private:
    double H_;
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

    void deformLantern(DeformCallback &dc);
    void buildObjective(const Eigen::VectorXd &q, double &f, Eigen::VectorXd &Df, std::vector<T> &Hf);
    void buildConstraints(const Eigen::VectorXd &q, Eigen::VectorXd &g, std::vector<T> &Dg, std::vector<std::vector<T> > &Hg);
    void buildInversionConstraints(const Eigen::VectorXd &q, Eigen::VectorXd &h, std::vector<T> &Dh, std::vector<std::vector<T> > &Hh);
    void radiusOverlapConstraint(const Eigen::VectorXd &q, OMMesh::HalfedgeHandle heh, double &g, std::vector<std::pair<int, double> > &Dg, std::vector<T> &Hg);
    bool shouldCollapseEdges(const Eigen::VectorXd &q);
    int findCollapsibleEdge(const Eigen::VectorXd &q);

    Eigen::Vector2d materialCenter();
    double materialRadius();
    void renderMaterial();

    void repopulateDOFs(const Eigen::VectorXd &q);

private:
    DevelopableMesh(const DevelopableMesh &other);
    DevelopableMesh &operator=(const DevelopableMesh &other);

    std::vector<Boundary> boundaries_;
    MaterialMesh *material_;

    Eigen::Vector3d point2Vector(OMMesh::Point pt);

    void centerCylinder();

    void collapseEdges();

    bool canCollapseEdge(int eid);
    void collapseEdge(int eid);

    double equalityConstraintViolation(const Eigen::VectorXd &q);
    int activeInequalityConstraints(const Eigen::VectorXd &q);
    void buildConstraintBasis(const Eigen::VectorXd &q, std::vector<Eigen::VectorXd> &normalspace, std::vector<Eigen::VectorXd> &tangentspace);
    void checkConstrainedHessian(const Eigen::VectorXd &q);

    void flushOutNANs(const Eigen::VectorXd &q);
    bool hasNANs(const Eigen::VectorXd &v);
    bool hasNANs(const std::vector<T> &M);
    bool hasNANs(const std::vector<std::vector<T> > &H);
};

#endif // DEVELOPABLEMESH_H
