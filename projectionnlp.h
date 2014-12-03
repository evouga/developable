#ifndef PROJECTIONNLP_H
#define PROJECTIONNLP_H

#include "coin/IpTNLP.hpp"
#include <Eigen/Core>
#include <Eigen/Sparse>

typedef Eigen::Triplet<double> T;

class DevelopableMesh;
class DeformCallback;

class ProjectionNLP : public Ipopt::TNLP
{
public:
    ProjectionNLP(DevelopableMesh &dm, DeformCallback &dc, Eigen::VectorXd startq) : dm_(dm), dc_(dc), startq_(startq) {}

    virtual bool get_nlp_info(Ipopt::Index &n, Ipopt::Index &m, Ipopt::Index &nnz_jac_g, Ipopt::Index &nnz_h_lag, IndexStyleEnum &index_style);
    virtual bool get_bounds_info(Ipopt::Index n, Ipopt::Number *x_l, Ipopt::Number *x_u, Ipopt::Index m, Ipopt::Number *g_l, Ipopt::Number *g_u);
    virtual bool get_starting_point(Ipopt::Index n, bool init_x, Ipopt::Number *x, bool init_z, Ipopt::Number *z_L, Ipopt::Number *z_U, Ipopt::Index m, bool init_lambda, Ipopt::Number *lambda);
    virtual bool eval_f(Ipopt::Index n, const Ipopt::Number *x, bool new_x, Ipopt::Number &obj_value);
    virtual bool eval_grad_f(Ipopt::Index n, const Ipopt::Number *x, bool new_x, Ipopt::Number *grad_f);
    virtual bool eval_g(Ipopt::Index n, const Ipopt::Number *x, bool new_x, Ipopt::Index m, Ipopt::Number *g);
    virtual bool eval_jac_g(Ipopt::Index n, const Ipopt::Number *x, bool new_x, Ipopt::Index m, Ipopt::Index nele_jac, Ipopt::Index *iRow, Ipopt::Index *jCol, Ipopt::Number *values);
    virtual bool eval_h(Ipopt::Index n, const Ipopt::Number *x, bool new_x, Ipopt::Number obj_factor, Ipopt::Index m, const Ipopt::Number *lambda, bool new_lambda, Ipopt::Index nele_hess, Ipopt::Index *iRow, Ipopt::Index *jCol, Ipopt::Number *values);
    virtual void finalize_solution(Ipopt::SolverReturn status, Ipopt::Index n, const Ipopt::Number *x, const Ipopt::Number *z_L, const Ipopt::Number *z_U, Ipopt::Index m, const Ipopt::Number *g, const Ipopt::Number *lambda, Ipopt::Number obj_value, const Ipopt::IpoptData *ip_data, Ipopt::IpoptCalculatedQuantities *ip_cq);

    virtual bool intermediate_callback(Ipopt::AlgorithmMode mode,
                                       Ipopt::Index iter,
                                       Ipopt::Number obj_value,
                                       Ipopt::Number inf_pr,
                                       Ipopt::Number inf_du,
                                       Ipopt::Number mu,
                                       Ipopt::Number d_norm,
                                       Ipopt::Number regularization_size,
                                       Ipopt::Number alpha_du,
                                       Ipopt::Number alpha_pr,
                                       Ipopt::Index ls_trials,
                                       const Ipopt::IpoptData *ip_data,
                                       Ipopt::IpoptCalculatedQuantities *ip_cq);

    Eigen::VectorXd getFinalQ() {return startq_;}

private:
    DevelopableMesh &dm_;
    DeformCallback &dc_;
    Eigen::VectorXd startq_;
};


#endif // PROJECTIONNLP_H
