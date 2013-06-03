#ifndef DEVELOPABLEMESH_H
#define DEVELOPABLEMESH_H

#include <vector>
#include <string>
#include <map>
#include "mesh.h"
#include <Eigen/Sparse>
#include <fadiff.h>
#include "autodiffTemplates.h"
#include "materialmesh.h"

typedef Eigen::Triplet<double> T;

class DeformCallback
{
public:
    virtual void repaintCallback()=0;
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

    enum CriticalPointType {CPT_INDETERMINATE, CPT_SADDLE, CPT_MAXIMUM, CPT_MINIMUM};

    void buildSchwarzLantern(double r, double h, int n, int m, double angle);
    virtual bool loadFromStream(std::istream &is);
    virtual bool saveToStream(std::ostream &os);

    void projectOntoConstraintManifold(DeformCallback &dc);
    void crushLantern(DeformCallback &dc, double dt);
    void buildObjective(const Eigen::VectorXd &q, double &f, Eigen::VectorXd &Df, std::vector<T> &Hf);
    void buildConstraints(const Eigen::VectorXd &q, Eigen::VectorXd &g, std::vector<T> &Dg, std::vector<std::vector<T> > &Hg);
    void buildInversionConstraints(const Eigen::VectorXd &q, Eigen::VectorXd &h, std::vector<T> &Dh, std::vector<std::vector<T> > &Hh);
    void radiusOverlapConstraint(const Eigen::VectorXd &q, OMMesh::HalfedgeHandle heh, double &g, std::vector<std::pair<int, double> > &Dg, std::vector<T> &Hg);
    bool shouldCollapseEdges(const Eigen::VectorXd &q);
    int findCollapsibleEdge(const Eigen::VectorXd &q);

    Eigen::Vector2d materialCenter();
    double materialRadius();
    void renderMaterial();

    void gatherDOFs(Eigen::VectorXd &q, Eigen::VectorXd &v);
    void repopulateDOFs(const Eigen::VectorXd &q, const Eigen::VectorXd &v);

    double equalityConstraintViolation(const Eigen::VectorXd &q);

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

    int activeInequalityConstraints(const Eigen::VectorXd &q);
    void buildConstraintBasis(const Eigen::VectorXd &q, std::vector<Eigen::VectorXd> &normalspace, std::vector<Eigen::VectorXd> &tangentspace);
    CriticalPointType checkConstrainedHessian(const Eigen::VectorXd &q, std::vector<Eigen::VectorXd> &negdirs);
    void gatherInfo(const Eigen::VectorXd &q);

    void flushOutNANs(const Eigen::VectorXd &q);
    void perturbConfiguration(Eigen::VectorXd &q, const std::vector<Eigen::VectorXd> &dirs, double mag);
};

#endif // DEVELOPABLEMESH_H
