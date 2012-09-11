#ifndef DEVELOPABLEMESH_H
#define DEVELOPABLEMESH_H

#include <vector>
#include <string>
#include "mesh.h"
#include <Eigen/Sparse>
#include <fadiff.h>
#define HAVE_CSTDDEF
#include <coin/IpIpoptApplication.hpp>

typedef Eigen::Triplet<double> T;

struct BoundaryCurve
{
    std::vector<int> edges;
    double arclength;
    double height;
    double targetheight;
};

class DevelopableMesh : public Mesh
{
public:
    DevelopableMesh();

    void buildSchwarzLantern(double r, double h, int n, int m, double angle);
    virtual bool loadMesh(const std::string &filename);

    void getBoundaryHeights(std::vector<double> &heights);

    void deformLantern(int maxiters);

private:
    std::vector<BoundaryCurve> boundaries_;
    double surfacearea_;

    void identifyBoundaries();
    void calculateSurfaceArea();
    Eigen::Vector3d point2Vector(OMMesh::Point pt);

    void buildObjective(const Eigen::VectorXd &q, double &f, Eigen::VectorXd &Df, std::vector<T> &Hf);
    void buildConstraints(const Eigen::VectorXd &q, Eigen::VectorXd &g, std::vector<T> &Dg, std::vector<std::vector<T> > &Hg);

    bool collapseShortEdges();
    void computeCreaseWidths(std::vector<double> &widths);
    double turningAngle(int edge);

    void collapseEdge(int edgeidx);
    bool canCollapseEdge(int edgeidx);
    void centerCylinder();

    template<class T> static T norm(T *v);
    template<class T> static void cross(T *v1, T *v2, T *result);
    template<class T> static void normalize(T *v);
    template<class T> static T dot(T *v1, T *v2);

    friend class IpoptSolver;
};

class IpoptSolver : public Ipopt::TNLP
{
public:
    IpoptSolver(int n, int m, int nnz_j, int nnz_h, Eigen::VectorXd &initq, DevelopableMesh &mesh);

    virtual bool get_nlp_info(Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index& nnz_jac_g,
                              Ipopt::Index& nnz_h_lag, IndexStyleEnum& index_style);

    virtual bool get_bounds_info(Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u,
                                 Ipopt::Index m, Ipopt::Number* g_l, Ipopt::Number* g_u);

    virtual bool get_starting_point(Ipopt::Index n, bool init_x, Ipopt::Number* x,
                                    bool init_z, Ipopt::Number* z_L, Ipopt::Number* z_U,
                                    Ipopt::Index m, bool init_lambda, Ipopt::Number* lambda);

    virtual bool eval_f(Ipopt::Index n, const Ipopt::Number* x,
                        bool new_x, Ipopt::Number& obj_value);

    virtual bool eval_grad_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                             Ipopt::Number* grad_f);

    virtual bool eval_g(Ipopt::Index n, const Ipopt::Number* x,
                        bool new_x, Ipopt::Index m, Ipopt::Number* g);

    virtual bool eval_jac_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                            Ipopt::Index m, Ipopt::Index nele_jac, Ipopt::Index* iRow,
                            Ipopt::Index *jCol, Ipopt::Number* values);

    virtual bool eval_h(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                        Ipopt::Number obj_factor, Ipopt::Index m, const Ipopt::Number* lambda,
                        bool new_lambda, Ipopt::Index nele_hess, Ipopt::Index* iRow,
                        Ipopt::Index* jCol, Ipopt::Number* values);

    virtual void finalize_solution(Ipopt::SolverReturn status, Ipopt::Index n,
                                   const Ipopt::Number* x, const Ipopt::Number* z_L,
                                   const Ipopt::Number* z_U, Ipopt::Index m, const Ipopt::Number* g,
                                   const Ipopt::Number* lambda, Ipopt::Number obj_value,
                                   const Ipopt::IpoptData* ip_data,
                                   Ipopt::IpoptCalculatedQuantities* ip_cq);

    Eigen::VectorXd &getQ() {return initq_;}
private:
    int n_, m_;
    int nnz_j_, nnz_h_;
    Eigen::VectorXd initq_;
    DevelopableMesh &mesh_;
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
