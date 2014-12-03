#include "projectionnlp.h"
#include "developablemesh.h"
#include "coin/IpOrigIpoptNLP.hpp"
#include "coin/IpIpoptCalculatedQuantities.hpp"
#include "coin/IpTNLPAdapter.hpp"
#include "coin/IpIpoptData.hpp"

using namespace Eigen;
using namespace std;

bool ProjectionNLP::get_nlp_info(Ipopt::Index &n, Ipopt::Index &m, Ipopt::Index &nnz_jac_g, Ipopt::Index &nnz_h_lag, IndexStyleEnum &index_style)
{
    n = startq_.size();
    nnz_h_lag = n;

    VectorXd g;
    vector<T> Dg;
    vector<vector<T> > Hg;
    dm_.buildConstraints(startq_, g, Dg, Hg);

//    VectorXd h;
//    vector<T> Dh;
//    vector<vector<T> > Hh;
//    dm_.buildInversionConstraints(startq_, h, Dh, Hh);

//    m = g.size() + h.size();
    m = g.size();
//    nnz_jac_g = Dg.size() + Dh.size();
    nnz_jac_g = Dg.size();

    for(int i=0; i<(int)Hg.size(); i++)
        nnz_h_lag += Hg[i].size();
//    for(int i=0; i<(int)Hh.size(); i++)
//        nnz_h_lag += Hh[i].size();

    index_style = C_STYLE;

    return true;
}

bool ProjectionNLP::get_bounds_info(Ipopt::Index n, Ipopt::Number *x_l, Ipopt::Number *x_u, Ipopt::Index m, Ipopt::Number *g_l, Ipopt::Number *g_u)
{
    for(int i=0; i<n; i++)
    {
        x_l[i] = -std::numeric_limits<double>::infinity();
        x_u[i] = std::numeric_limits<double>::infinity();
    }

    VectorXd g;
    vector<T> Dg;
    vector<vector<T> > Hg;
    dm_.buildConstraints(startq_, g, Dg, Hg);

    for(int i=0; i<(int)g.size(); i++)
    {
        g_l[i] = 0;
        g_u[i] = 0;
    }

    return true;
}

bool ProjectionNLP::get_starting_point(Ipopt::Index n, bool init_x, Ipopt::Number *x, bool init_z, Ipopt::Number *, Ipopt::Number *, Ipopt::Index , bool init_lambda, Ipopt::Number *)
{
    assert(init_x == true);
    assert(init_z == false);
    assert(init_lambda == false);

    for(int i=0; i<n; i++)
        x[i] = startq_[i];

    return true;
}

bool ProjectionNLP::eval_f(Ipopt::Index n, const Ipopt::Number *x, bool , Ipopt::Number &obj_value)
{
    VectorXd q(n);
    for(int i=0; i<n; i++)
    {
        q[i] = x[i];
    }
    obj_value = 0.5*(q-startq_).squaredNorm();

    return true;
}

bool ProjectionNLP::eval_grad_f(Ipopt::Index n, const Ipopt::Number *x, bool , Ipopt::Number *grad_f)
{
    VectorXd q(n);
    for(int i=0; i<n; i++)
    {
        q[i] = x[i];
    }
    double val = 0.0;
    VectorXd Df(n);
    vector<T> Hf;
    dm_.elasticEnergy(q,val,Df,Hf);
    obj_value = val;

    for(int i=0; i <n; i++)
        grad_f[i] = Df[i];
    return true;
}

bool ProjectionNLP::eval_g(Ipopt::Index n, const Ipopt::Number *x, bool , Ipopt::Index , Ipopt::Number *g)
{
    VectorXd q(n);
    for(int i=0; i<n; i++)
    {
        q[i] = x[i];
    }
    VectorXd Vg;
    vector<T> dummyM;
    vector<vector<T> > dummyH;
    dm_.buildConstraints(q, Vg, dummyM, dummyH);
    for(int i=0; i<(int)Vg.size(); i++)
        g[i] = Vg[i];
//    VectorXd h;
//    dm_.buildInversionConstraints(q, h, dummyM, dummyH);
//    for(int i=0; i<(int)h.size(); i++)
//    {
//        g[i+Vg.size()] = h[i];
//    }

    return true;
}

bool ProjectionNLP::eval_jac_g(Ipopt::Index n, const Ipopt::Number *x, bool , Ipopt::Index, Ipopt::Index , Ipopt::Index *iRow, Ipopt::Index *jCol, Ipopt::Number *values)
{
    if(iRow != NULL && jCol != NULL)
    {
        VectorXd g;
        vector<T> Dg;
        vector<vector<T> > Hg;
        dm_.buildConstraints(startq_, g, Dg, Hg);

        for(int i=0; i<(int)Dg.size(); i++)
        {
            iRow[i] = Dg[i].row();
            jCol[i] = Dg[i].col();
        }
//        VectorXd h;
//        vector<T> Dh;
//        vector<vector<T> > Hh;
//        dm_.buildInversionConstraints(startq_, h, Dh, Hh);

//        for(int i=0; i<(int)Dh.size(); i++)
//        {
//            iRow[i+Dg.size()] = g.size() + Dh[i].row();
//            jCol[i+Dg.size()] = Dh[i].col();
//        }
    }
    else
    {
        VectorXd q(n);
        for(int i=0; i<n; i++)
            q[i] = x[i];

        VectorXd g;
        vector<T> Dg;
        vector<vector<T> > Hg;
        dm_.buildConstraints(q, g, Dg, Hg);
        for(int i=0; i<(int)Dg.size(); i++)
            values[i] = Dg[i].value();
//        VectorXd h;
//        vector<T> Dh;
//        vector<vector<T> > Hh;
//        dm_.buildInversionConstraints(q, h, Dh, Hh);
//        for(int i=0; i<(int)Dh.size(); i++)
//            values[i+Dg.size()] = Dh[i].value();
    }

    return true;
}

bool ProjectionNLP::eval_h(Ipopt::Index n, const Ipopt::Number *x, bool , Ipopt::Number obj_factor, Ipopt::Index m, const Ipopt::Number *lambda, bool , Ipopt::Index nele_hess, Ipopt::Index *iRow, Ipopt::Index *jCol, Ipopt::Number *values)
{
    if(iRow != NULL && jCol != NULL)
    {
        VectorXd g;
        vector<T> Dg;
        vector<vector<T> > Hg;
        dm_.buildConstraints(startq_, g, Dg, Hg);
//        VectorXd h;
//        vector<T> Dh;
//        vector<vector<T> > Hh;
//        dm_.buildInversionConstraints(startq_, h, Dh, Hh);
//        assert((int)Hg.size() + (int)Hh.size() == m);
        assert((int)Hg.size() == m);

        int row=0;
        for(int i=0; i<n; i++)
        {
            iRow[row] = i;
            jCol[row] = i;
            row++;
        }

        for(int i=0; i<(int)Hg.size(); i++)
        {
            for(int j=0; j<(int)Hg[i].size(); j++)
            {
                iRow[row] = Hg[i][j].row();
                jCol[row] = Hg[i][j].col();
                row++;
            }
        }
//        for(int i=0; i<(int)Hh.size(); i++)
//        {
//            for(int j=0; j<(int)Hh[i].size(); j++)
//            {
//                iRow[row] = Hh[i][j].row();
//                jCol[row] = Hh[i][j].col();
//                row++;
//            }
//        }
        assert(row == nele_hess);
    }
    else
    {
        VectorXd q(n);
        for(int i=0; i<n; i++)
            q[i] = x[i];

        int row=0;

        VectorXd g;
        vector<T> Dg;
        vector<vector<T> > Hg;
        dm_.buildConstraints(q, g, Dg, Hg);
//        VectorXd h;
//        vector<T> Dh;
//        vector<vector<T> > Hh;
//        dm_.buildInversionConstraints(q, h, Dh, Hh);

        for(int i=0; i<n; i++)
        {
            values[row] = obj_factor;
            row++;
        }

        for(int i=0; i<(int)Hg.size(); i++)
        {
            for(int j=0; j<(int)Hg[i].size(); j++)
            {
                values[row] = lambda[i]*Hg[i][j].value();
                row++;
            }
        }
//        for(int i=0; i<(int)Hh.size(); i++)
//        {
//            for(int j=0; j<(int)Hh[i].size(); j++)
//            {
//                values[row] = lambda[i+Hg.size()]*Hh[i][j].value();
//                row++;
//            }
//        }
        assert(row == nele_hess);
    }

    return true;
}

void ProjectionNLP::finalize_solution(Ipopt::SolverReturn status, Ipopt::Index n, const Ipopt::Number *x, const Ipopt::Number *, const Ipopt::Number *, Ipopt::Index, const Ipopt::Number *, const Ipopt::Number *, Ipopt::Number, const Ipopt::IpoptData *, Ipopt::IpoptCalculatedQuantities *)
{
    cout << "Status: " << status << endl;

    for(int i=0; i<n; i++)
        startq_[i] = x[i];
}

bool ProjectionNLP::intermediate_callback(Ipopt::AlgorithmMode , Ipopt::Index , Ipopt::Number , Ipopt::Number, Ipopt::Number, Ipopt::Number, Ipopt::Number, Ipopt::Number, Ipopt::Number, Ipopt::Number, Ipopt::Index, const Ipopt::IpoptData *ip_data, Ipopt::IpoptCalculatedQuantities *ip_cq)
{
    return true;
    /*
    assert(ip_cq);
    Ipopt::OrigIpoptNLP *orignlp = dynamic_cast<Ipopt::OrigIpoptNLP*>(GetRawPtr(ip_cq->GetIpoptNLP()));
    if(!orignlp)
    {
        return true;
    }
    Ipopt::TNLPAdapter *tnlp_adapter = dynamic_cast<Ipopt::TNLPAdapter*>(GetRawPtr(orignlp->nlp()));
    assert(tnlp_adapter);

    double *x = new double[startq_.size()];
    tnlp_adapter->ResortX(*ip_data->curr()->x(), x);
    Eigen::VectorXd q;
    Eigen::VectorXd v;
    dm_.gatherDOFs(q, v);
    for(int i=0; i<(int)startq_.size(); i++)
        q[i] = x[i];
    delete[] x;

    dm_.repopulateDOFs(q,v);
    dc_.repaintCallback();

    return true;*/
}
