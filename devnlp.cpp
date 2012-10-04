#include "devnlp.h"
#include "developablemesh.h"

using namespace Eigen;
using namespace std;



bool DevTLNP::get_nlp_info(Ipopt::Index &n, Ipopt::Index &m, Ipopt::Index &nnz_jac_g, Ipopt::Index &nnz_h_lag, IndexStyleEnum &index_style)
{
    double f;
    VectorXd Df;
    vector<T> Hf;
    dm_.buildObjective(startq_, f, Df, Hf);
    n = startq_.size();
    nnz_h_lag = Hf.size();

    VectorXd g;
    vector<T> Dg;
    vector<vector<T> > Hg;
    dm_.buildConstraints(startq_, g, Dg, Hg);

    VectorXd h;
    vector<T> Dh;
    vector<vector<T> > Hh;
    dm_.buildInversionConstraints(startq_, h, Dh, Hh);

    m = g.size() + h.size();
    nnz_jac_g = Dg.size() + Dh.size();

    for(int i=0; i<(int)Hg.size(); i++)
        nnz_h_lag += Hg[i].size();
    for(int i=0; i<(int)Hh.size(); i++)
        nnz_h_lag += Hh[i].size();

    index_style = C_STYLE;

    return true;
}

bool DevTLNP::get_bounds_info(Ipopt::Index n, Ipopt::Number *x_l, Ipopt::Number *x_u, Ipopt::Index m, Ipopt::Number *g_l, Ipopt::Number *g_u)
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
    for(int i=(int)g.size(); i < m; i++)
    {
        g_l[i] = 0;
        g_u[i] = std::numeric_limits<double>::infinity();
    }

    return true;
}

bool DevTLNP::get_starting_point(Ipopt::Index n, bool init_x, Ipopt::Number *x, bool init_z, Ipopt::Number *, Ipopt::Number *, Ipopt::Index , bool init_lambda, Ipopt::Number *)
{
    assert(init_x == true);
    assert(init_z == false);
    assert(init_lambda == false);

    for(int i=0; i<n; i++)
        x[i] = startq_[i];

    return true;
}

bool DevTLNP::eval_f(Ipopt::Index n, const Ipopt::Number *x, bool , Ipopt::Number &obj_value)
{
    VectorXd q(n);
    for(int i=0; i<n; i++)
    {
        q[i] = x[i];
    }
    double f;
    VectorXd dummy;
    vector<T> dummyM;
    dm_.buildObjective(q, f, dummy, dummyM);
    obj_value = f;

    return true;
}

bool DevTLNP::eval_grad_f(Ipopt::Index n, const Ipopt::Number *x, bool , Ipopt::Number *grad_f)
{
    VectorXd q(n);
    for(int i=0; i<n; i++)
    {
        q[i] = x[i];
    }
    double f;
    VectorXd Df;
    vector<T> dummyM;
    dm_.buildObjective(q, f, Df, dummyM);
    for(int i=0; i<n; i++)
        grad_f[i] = Df[i];

    return true;
}

bool DevTLNP::eval_g(Ipopt::Index n, const Ipopt::Number *x, bool , Ipopt::Index , Ipopt::Number *g)
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
    VectorXd h;
    dm_.buildInversionConstraints(q, h, dummyM, dummyH);
    for(int i=0; i<(int)h.size(); i++)
    {
        g[i+Vg.size()] = h[i];
    }

    return true;
}

bool DevTLNP::eval_jac_g(Ipopt::Index n, const Ipopt::Number *x, bool , Ipopt::Index, Ipopt::Index , Ipopt::Index *iRow, Ipopt::Index *jCol, Ipopt::Number *values)
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
        VectorXd h;
        vector<T> Dh;
        vector<vector<T> > Hh;
        dm_.buildInversionConstraints(startq_, h, Dh, Hh);

        for(int i=0; i<(int)Dh.size(); i++)
        {
            iRow[i+Dg.size()] = g.size() + Dh[i].row();
            jCol[i+Dg.size()] = Dh[i].col();
        }
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
        VectorXd h;
        vector<T> Dh;
        vector<vector<T> > Hh;
        dm_.buildInversionConstraints(q, h, Dh, Hh);
        for(int i=0; i<(int)Dh.size(); i++)
            values[i+Dg.size()] = Dh[i].value();
    }

    return true;
}

bool DevTLNP::eval_h(Ipopt::Index n, const Ipopt::Number *x, bool , Ipopt::Number obj_factor, Ipopt::Index m, const Ipopt::Number *lambda, bool , Ipopt::Index nele_hess, Ipopt::Index *iRow, Ipopt::Index *jCol, Ipopt::Number *values)
{
    if(iRow != NULL && jCol != NULL)
    {
        double f;
        VectorXd Df;
        vector<T> Hf;
        dm_.buildObjective(startq_, f, Df, Hf);
        VectorXd g;
        vector<T> Dg;
        vector<vector<T> > Hg;
        dm_.buildConstraints(startq_, g, Dg, Hg);
        VectorXd h;
        vector<T> Dh;
        vector<vector<T> > Hh;
        dm_.buildInversionConstraints(startq_, h, Dh, Hh);
        assert((int)Hg.size() + (int)Hh.size() == m);

        int row=0;
        for(int i=0; i<(int)Hf.size(); i++)
        {
            iRow[row] = Hf[i].row();
            jCol[row] = Hf[i].col();
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
        for(int i=0; i<(int)Hh.size(); i++)
        {
            for(int j=0; j<(int)Hh[i].size(); j++)
            {
                iRow[row] = Hh[i][j].row();
                jCol[row] = Hh[i][j].col();
                row++;
            }
        }
        assert(row == nele_hess);
    }
    else
    {
        VectorXd q(n);
        for(int i=0; i<n; i++)
            q[i] = x[i];
        double f;
        VectorXd Df;
        vector<T> Hf;
        dm_.buildObjective(q, f, Df, Hf);
        int row=0;
        for(int i=0; i<(int)Hf.size(); i++)
        {
            values[row] = obj_factor * Hf[i].value();
            row++;
        }

        VectorXd g;
        vector<T> Dg;
        vector<vector<T> > Hg;
        dm_.buildConstraints(q, g, Dg, Hg);
        VectorXd h;
        vector<T> Dh;
        vector<vector<T> > Hh;
        dm_.buildInversionConstraints(q, h, Dh, Hh);
        for(int i=0; i<(int)Hg.size(); i++)
        {
            for(int j=0; j<(int)Hg[i].size(); j++)
            {
                values[row] = lambda[i]*Hg[i][j].value();
                row++;
            }
        }
        for(int i=0; i<(int)Hh.size(); i++)
        {
            for(int j=0; j<(int)Hh[i].size(); j++)
            {
                values[row] = lambda[i+Hg.size()]*Hh[i][j].value();
                row++;
            }
        }
        assert(row == nele_hess);
    }

    return true;
}

void DevTLNP::finalize_solution(Ipopt::SolverReturn status, Ipopt::Index n, const Ipopt::Number *x, const Ipopt::Number *, const Ipopt::Number *, Ipopt::Index m, const Ipopt::Number *g, const Ipopt::Number *, Ipopt::Number obj_value, const Ipopt::IpoptData *, Ipopt::IpoptCalculatedQuantities *)
{
    cout << "Status: " << status << endl;

    for(int i=0; i<n; i++)
        startq_[i] = x[i];

    double viol = 0.0;
    for(int i=0; i<m; i++)
        viol += g[i]*g[i];
    cout << "Constraint violation " << sqrt(viol) << endl;

    cout << "Objective value " << obj_value << endl;
}
