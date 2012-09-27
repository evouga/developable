#ifndef DEVELOPABLEMESH_H
#define DEVELOPABLEMESH_H

#include <vector>
#include <string>
#include <map>
#include "mesh.h"
#include <Eigen/Sparse>
#include <fadiff.h>
#include "coin/IpTNLP.hpp"

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

private:
    DevelopableMesh(const DevelopableMesh &other);
    DevelopableMesh &operator=(const DevelopableMesh &other);
    std::vector<Boundary> boundaries_;
    MaterialMesh *material_;
    Eigen::Vector3d point2Vector(OMMesh::Point pt);

    void centerCylinder();

};

class DevTLNP : public Ipopt::TNLP
{
public:
    DevTLNP(DevelopableMesh &dm, Eigen::VectorXd startq) : dm_(dm), startq_(startq) {}

    virtual bool get_nlp_info(Ipopt::Index &n, Ipopt::Index &m, Ipopt::Index &nnz_jac_g, Ipopt::Index &nnz_h_lag, IndexStyleEnum &index_style);
    virtual bool get_bounds_info(Ipopt::Index n, Ipopt::Number *x_l, Ipopt::Number *x_u, Ipopt::Index m, Ipopt::Number *g_l, Ipopt::Number *g_u);
    virtual bool get_starting_point(Ipopt::Index n, bool init_x, Ipopt::Number *x, bool init_z, Ipopt::Number *z_L, Ipopt::Number *z_U, Ipopt::Index m, bool init_lambda, Ipopt::Number *lambda);
    virtual bool eval_f(Ipopt::Index n, const Ipopt::Number *x, bool new_x, Ipopt::Number &obj_value);
    virtual bool eval_grad_f(Ipopt::Index n, const Ipopt::Number *x, bool new_x, Ipopt::Number *grad_f);
    virtual bool eval_g(Ipopt::Index n, const Ipopt::Number *x, bool new_x, Ipopt::Index m, Ipopt::Number *g);
    virtual bool eval_jac_g(Ipopt::Index n, const Ipopt::Number *x, bool new_x, Ipopt::Index m, Ipopt::Index nele_jac, Ipopt::Index *iRow, Ipopt::Index *jCol, Ipopt::Number *values);
    virtual bool eval_h(Ipopt::Index n, const Ipopt::Number *x, bool new_x, Ipopt::Number obj_factor, Ipopt::Index m, const Ipopt::Number *lambda, bool new_lambda, Ipopt::Index nele_hess, Ipopt::Index *iRow, Ipopt::Index *jCol, Ipopt::Number *values);
    virtual void finalize_solution(Ipopt::SolverReturn status, Ipopt::Index n, const Ipopt::Number *x, const Ipopt::Number *z_L, const Ipopt::Number *z_U, Ipopt::Index m, const Ipopt::Number *g, const Ipopt::Number *lambda, Ipopt::Number obj_value, const Ipopt::IpoptData *ip_data, Ipopt::IpoptCalculatedQuantities *ip_cq);

    Eigen::VectorXd getFinalQ() {return startq_;}

private:
    DevelopableMesh &dm_;
    Eigen::VectorXd startq_;
};

template<typename T> void diff(const T *vec1, const T *vec2, T *difference)
{
    for(int i=0; i<3; i++)
        difference[i] = vec1[i]-vec2[i];
}

template<typename T> T norm(const T *vec)
{
    T result = 0;
    for(int i=0; i<3; i++)
        result += vec[i]*vec[i];
    return sqrt(result);
}

template<typename T> void cross(const T *vec1, const T *vec2, T *result)
{
    result[0] = vec1[1]*vec2[2] - vec1[2]*vec2[1];
    result[1] = vec1[2]*vec2[0] - vec1[0]*vec2[2];
    result[2] = vec1[0]*vec2[1] - vec1[1]*vec2[0];
}

template<typename T> void normalize(T *vec)
{
    T n = norm(vec);
    for(int i=0; i<3; i++)
        vec[i] /= n;
}

template<typename T> T dot(const T *vec1, const T *vec2)
{
    T result = 0;
    for(int i=0; i<3; i++)
        result += vec1[i]*vec2[i];
    return result;
}

#endif // DEVELOPABLEMESH_H
