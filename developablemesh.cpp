#include "developablemesh.h"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Geometry>
#include <set>

#include <fadbad.h>
#include <fadiff.h>
#include "coin/IpIpoptApplication.hpp"

using namespace std;
using namespace Eigen;
using namespace fadbad;
using namespace Ipopt;

const double PI = 3.1415926535898;


void MaterialMesh::setMaterialEdge(int embeddedEdge, int materialEdge)
{
    assert(embeddedEdge != -1);
    assert(materialEdge != -1);
    embedge2matedge_[embeddedEdge] = materialEdge;
}

int MaterialMesh::materialEdge(int embeddedEdge)
{
    map<int, int>::iterator it = embedge2matedge_.find(embeddedEdge);
    if(it != embedge2matedge_.end())
        return it->second;
    return -1;
}


DevelopableMesh::DevelopableMesh() : material_(new MaterialMesh(0,0))
{

}

DevelopableMesh::~DevelopableMesh()
{
    delete material_;
}

bool DevelopableMesh::loadMesh(const string &)
{
    assert(false);
    return false;
}

void DevelopableMesh::buildSchwarzLantern(double r, double h, int n, int m, double angle)
{
    mesh_ = OMMesh();
    double alpha = angle;
    double a = 2*r*sin(PI/n);
    double b = sqrt(4*r*r*sin(alpha/2.0)*sin(alpha/2.0) + h*h/m/m);
    double c = sqrt(4*r*r*sin(PI/n + alpha/2.0)*sin(PI/n + alpha/2.0) + h*h/m/m);
    double costilt = (a*a+b*b-c*c)/(2*a*b);
    double dpx = -b*costilt;


    double s = 0.5*(a+b+c);
    double area = sqrt(s*(s-a)*(s-b)*(s-c));
    double dpy = 2*area/a;

    double W = n*a;
    double H = 2*n*m*area/W;

    assert(fabs(H-dpy*m) < 1e-6);

    delete material_;
    material_ = new MaterialMesh(W,H);
    boundaries_.clear();

    Boundary bottom, top;

    for(int i=0; i<=m; i++)
    {
        double z = h*i/m;
        double basepx = i*dpx;
        double py = i*dpy;

        int firstvid = material_->getMesh().n_vertices();

        for(int j=0; j<n; j++)
        {
            double x = r*cos(2*PI*(j/double(n)) + i*alpha);
            double y = r*sin(2*PI*(j/double(n)) + i*alpha);

            OMMesh::Point newpt(x,y,z);
            if(i == 0)
            {
                bottom.bdryVerts.push_back(mesh_.n_vertices());
                bottom.bdryPos.push_back(Vector3d(x,y,z));
            }
            else if(i==m)
            {
                top.bdryVerts.push_back(mesh_.n_vertices());
                top.bdryPos.push_back(Vector3d(x,y,z));
            }
            mesh_.add_vertex(newpt);

            double px = j*a+basepx;
            OMMesh::Point newmatpt(px,py,0);
            if(i == 0 || i == m)
            {
                material_->getBoundaryVerts().push_back(pair<int, double>(material_->getMesh().n_vertices(), py));
            }
            material_->getMesh().add_vertex(newmatpt);

        }

        double px = n*a+basepx;
        OMMesh::Point newmatpt(px,py,0);
        int lastvid = material_->getMesh().n_vertices();
        material_->getMesh().add_vertex(newmatpt);
        material_->getIdentifiedVerts().push_back(pair<int,int>(firstvid,lastvid));
    }

    boundaries_.push_back(bottom);
    boundaries_.push_back(top);

    for(int i=0; i<m; i++)
    {
        for(int j=0; j<n; j++)
        {
            int fidx1 = i*n+j;
            int fidx2 = i*n + ((j+1) % n);
            int fidx3 = (i+1)*n+ ((j+1)%n);

            int midx1 = i*(n+1)+j;
            int midx2 = i*(n+1) + j + 1;
            int midx3 = (i+1)*(n+1) + j + 1;
            vector<OMMesh::VertexHandle> newface;
            vector<OMMesh::VertexHandle> newmatface;
            newface.push_back(mesh_.vertex_handle(fidx1));
            newface.push_back(mesh_.vertex_handle(fidx2));
            newface.push_back(mesh_.vertex_handle(fidx3));
            newmatface.push_back(material_->getMesh().vertex_handle(midx1));
            newmatface.push_back(material_->getMesh().vertex_handle(midx2));
            newmatface.push_back(material_->getMesh().vertex_handle(midx3));

            mesh_.add_face(newface);
            material_->getMesh().add_face(newmatface);

            material_->setMaterialEdge(findEdge(fidx1,fidx2), material_->findEdge(midx1,midx2));
            material_->setMaterialEdge(findEdge(fidx2,fidx3), material_->findEdge(midx2,midx3));
            material_->setMaterialEdge(findEdge(fidx3,fidx1), material_->findEdge(midx3,midx1));

            fidx2 = fidx3;
            midx2 = midx3;
            newface[1] = mesh_.vertex_handle(fidx2);
            newmatface[1] = material_->getMesh().vertex_handle(midx2);
            fidx3 = (i+1)*n + ((n + j) % n);
            midx3 = (i+1)*(n+1) + j;
            newface[2] = mesh_.vertex_handle(fidx3);
            newmatface[2] = material_->getMesh().vertex_handle(midx3);
            mesh_.add_face(newface);
            material_->getMesh().add_face(newmatface);

            material_->setMaterialEdge(findEdge(fidx1,fidx2), material_->findEdge(midx1,midx2));
            material_->setMaterialEdge(findEdge(fidx2,fidx3), material_->findEdge(midx2,midx3));
            material_->setMaterialEdge(findEdge(fidx3,fidx1), material_->findEdge(midx3,midx1));

        }
    }

    for(int i=0; i<(int)mesh_.n_edges(); i++)
    {
        double len1 = mesh_.calc_edge_length(mesh_.edge_handle(i));

        int medge = material_->materialEdge(i);

        double len2 = material_->getMesh().calc_edge_length(material_->getMesh().edge_handle(medge));

        assert(fabs(len1-len2) < 1e-6);
    }
}

Vector3d DevelopableMesh::point2Vector(OMMesh::Point pt)
{
    Vector3d result;
    for(int i=0; i<3; i++)
        result[i] = pt[i];
    return result;
}

void DevelopableMesh::centerCylinder()
{
    double totx=0;
    double toty=0;
    for(int i=0; i<(int)mesh_.n_vertices(); i++)
    {
        OMMesh::Point pt = mesh_.point(mesh_.vertex_handle(i));
        totx += pt[0];
        toty += pt[1];
    }
    totx /= mesh_.n_vertices();
    toty /= mesh_.n_vertices();

    for(int i=0; i<(int)mesh_.n_vertices(); i++)
    {
        mesh_.point(mesh_.vertex_handle(i))[0] -= totx;
        mesh_.point(mesh_.vertex_handle(i))[1] -= toty;
    }
}

void DevelopableMesh::deformLantern(int maxiters)
{
    int nembverts = mesh_.n_vertices();
    int nmatverts = material_->getMesh().n_vertices();

    VectorXd q(3*nembverts+2*nmatverts);
    for(int i=0; i<nembverts; i++)
    {
        for(int j=0; j<3; j++)
        {
            q[3*i+j] = mesh_.point(mesh_.vertex_handle(i))[j];
        }
    }
    for(int i=0; i<nmatverts; i++)
    {
        for(int j=0; j<2; j++)
        {
            q[3*nembverts+2*i+j] = material_->getMesh().point(material_->getMesh().vertex_handle(i))[j];
        }
    }

    VectorXd g;
    vector<T> Dg;
    vector<vector<T> > Hg;
    buildConstraints(q, g, Dg, Hg);
    cout << "Initial violation: " << g.norm() << endl;
    double f;
    VectorXd Df;
    vector<T> Hf;
    buildObjective(q, f, Df, Hf);
    cout << "Initial energy: " << f << endl;

    DevTLNP *mynlp = new DevTLNP(*this, q);
    IpoptApplication app;
    app.Options()->SetNumericValue("tol", 1e-9);
    //app.Options()->SetStringValue("mu_strategy", "adaptive");
    //app.Options()->SetStringValue("derivative_test", "second-order");
    app.Options()->SetStringValue("check_derivatives_for_naninf", "yes");
    ApplicationReturnStatus status;
    status = app.Initialize();
    status = app.OptimizeTNLP(mynlp);

    q = mynlp->getFinalQ();

    for(int i=0; i<nembverts; i++)
    {
        for(int j=0; j<3; j++)
        {
            mesh_.point(mesh_.vertex_handle(i))[j] = q[3*i+j];
        }
    }

    for(int i=0; i<nmatverts; i++)
    {
        for(int j=0; j<2; j++)
        {
            material_->getMesh().point(material_->getMesh().vertex_handle(i))[j] = q[3*nembverts+2*i+j];
        }
    }

    centerCylinder();
}

void DevelopableMesh::buildObjective(const Eigen::VectorXd &q, double &f, Eigen::VectorXd &Df, std::vector<T> &Hf)
{
    double nembverts = mesh_.n_vertices();
    double nmatverts = material_->getMesh().n_vertices();
    assert(q.size() == 3*nembverts + 2*nmatverts);
    VectorXd embq = q.segment(0,3*nembverts);

    f = 0;
    Df.resize(3*nembverts+2*nmatverts);
    Df.setZero();
    Hf.clear();

    for(int i=0; i<(int)mesh_.n_edges(); i++)
    {
        OMMesh::EdgeHandle eh = mesh_.edge_handle(i);
        if(mesh_.is_boundary(eh))
            continue;
        OMMesh::HalfedgeHandle heh = mesh_.halfedge_handle(eh, 0);
        int vertidx[4];
        OMMesh::VertexHandle v[4];
        F<double> vx[4][3];
        F<F<double> > vxx[4][3];

        int curvert = 0;
        v[curvert] = mesh_.from_vertex_handle(heh);
        vertidx[curvert] = v[curvert].idx();
        for(int i=0; i<3; i++)
        {
            vx[curvert][i] = embq[3*vertidx[curvert]+i];
            vx[curvert][i].diff(3*curvert + i, 12);
            vxx[curvert][i] = vx[curvert][i];
            vxx[curvert][i].diff(3*curvert+ i, 12);
        }

        curvert++;

        v[curvert] = mesh_.to_vertex_handle(heh);
        vertidx[curvert] = v[curvert].idx();
        for(int i=0; i<3; i++)
        {
            vx[curvert][i] = embq[3*vertidx[curvert]+i];
            vx[curvert][i].diff(3*curvert + i, 12);
            vxx[curvert][i] = vx[curvert][i];
            vxx[curvert][i].diff(3*curvert+ i, 12);
        }

        curvert++;

        v[curvert] = mesh_.to_vertex_handle(mesh_.next_halfedge_handle(heh));
        vertidx[curvert] = v[curvert].idx();
        for(int i=0; i<3; i++)
        {
            vx[curvert][i] = embq[3*vertidx[curvert]+i];
            vx[curvert][i].diff(3*curvert + i, 12);
            vxx[curvert][i] = vx[curvert][i];
            vxx[curvert][i].diff(3*curvert+ i, 12);
        }

        curvert++;

        v[curvert] = mesh_.to_vertex_handle(mesh_.next_halfedge_handle(mesh_.opposite_halfedge_handle(heh)));
        vertidx[curvert] = v[curvert].idx();
        for(int i=0; i<3; i++)
        {
            vx[curvert][i] = embq[3*vertidx[curvert]+i];
            vx[curvert][i].diff(3*curvert + i, 12);
            vxx[curvert][i] = vx[curvert][i];
            vxx[curvert][i].diff(3*curvert+ i, 12);
        }

        F<F<double> > e01[3];
        diff(vxx[1],vxx[0],e01);

        F<F<double> > L = norm(e01);

        F<F<double> > e02[3];
        F<F<double> > e03[3];
        diff(vxx[2],vxx[0],e02);
        diff(vxx[3],vxx[0],e03);

        F<F<double> > n0[3];
        F<F<double> > n1[3];
        cross(e01,e02,n0);
        cross(e03,e01,n1);
        F<F<double> > A0 = norm(n0);
        F<F<double> > A1 = norm(n1);
        normalize(n0);
        normalize(n1);
        F<F<double> > costheta = dot(n0,n1);

        //F<F<double> > Lpart = L*L/(A0+A1);
        F<F<double> > Lpart = pow(L, 1.0/3.0);
        //F<F<double> > anglepart = acos(costheta)*acos(costheta);
        F<F<double> > anglepart = pow(acos(costheta), 7.0/3.0);
        F<F<double> > stencilf = Lpart*anglepart;

        f += stencilf.val().val();
        for(int i=0; i<4; i++)
        {
            for(int j=0; j<3; j++)
            {
                Df[3*vertidx[i]+j] += stencilf.d(3*i+j).val();

                for(int k=0; k<4; k++)
                {
                    for(int l=0; l<3; l++)
                    {
                        if(3*vertidx[i]+j <= 3*vertidx[k]+l)
                            Hf.push_back(T(3*vertidx[i]+j, 3*vertidx[k]+l, stencilf.d(3*i+j).d(3*k+l)));
                    }
                }
            }
        }
    }
}

void DevelopableMesh::buildConstraints(const Eigen::VectorXd &q, Eigen::VectorXd &g, std::vector<T> &Dg, vector<vector<T> > &Hg)
{
    double nembverts = mesh_.n_vertices();
    double nmatverts = material_->getMesh().n_vertices();
    assert(q.size() == 3*nembverts + 2*nmatverts);
    VectorXd embq = q.segment(0,3*nembverts);
    VectorXd matq = q.segment(3*nembverts, 2*nmatverts);

    int numconstraints = 0;
    // coupling constraints between edge lengths
    //numconstraints += mesh_.n_edges();
    for(int i=0; i<(int)mesh_.n_edges(); i++)
    {
        if(!mesh_.is_boundary(mesh_.edge_handle(i)))
            numconstraints++;
    }

    // identified verts stay identified
    numconstraints += 2*material_->getIdentifiedVerts().size();

    // top and bottom verts have y values fixed
    numconstraints += material_->getBoundaryVerts().size();

    // top and bottom of embedded cylinder have y values fixed
    for(int i=0; i<(int)boundaries_.size(); i++)
        numconstraints += 3*boundaries_[i].bdryVerts.size();


    g.resize(numconstraints);
    g.setZero();
    Dg.clear();
    Hg.clear();

    int row=0;

    for(int i=0; i<(int)mesh_.n_edges(); i++)
    {
        OMMesh::EdgeHandle eh = mesh_.edge_handle(i);
        if(mesh_.is_boundary(eh))
            continue;
        OMMesh::HalfedgeHandle heh = mesh_.halfedge_handle(eh,0);
        int embids[2];
        embids[0] = mesh_.from_vertex_handle(heh).idx();
        embids[1] = mesh_.to_vertex_handle(heh).idx();

        int mateidx = material_->materialEdge(i);
        assert(mateidx != -1);
        OMMesh::EdgeHandle meh = material_->getMesh().edge_handle(mateidx);
        OMMesh::HalfedgeHandle mheh = material_->getMesh().halfedge_handle(meh,0);
        int matids[2];
        matids[0] = material_->getMesh().from_vertex_handle(mheh).idx();
        matids[1] = material_->getMesh().to_vertex_handle(mheh).idx();

        Vector3d embedge = embq.segment<3>(3*embids[0]) - embq.segment<3>(3*embids[1]);
        double len1 = embedge.squaredNorm();
        Vector2d matedge = matq.segment<2>(2*matids[0]) - matq.segment<2>(2*matids[1]);
        double len2 = matedge.squaredNorm();
        g[row] = (len1-len2);

        vector<T> Hgentry;

        for(int j=0; j<3; j++)
        {
            Dg.push_back(T(row, 3*embids[0]+j,2.0*embedge[j]));
            Dg.push_back(T(row, 3*embids[1]+j,-2.0*embedge[j]));

            Hgentry.push_back((T(3*embids[0]+j, 3*embids[0]+j, 2.0)));
            if(embids[0] < embids[1])
                Hgentry.push_back((T(3*embids[0]+j, 3*embids[1]+j, -2.0)));
            else
                Hgentry.push_back((T(3*embids[1]+j, 3*embids[0]+j, -2.0)));
            Hgentry.push_back((T(3*embids[1]+j, 3*embids[1]+j, 2.0)));
        }
        for(int j=0; j<2; j++)
        {
            Dg.push_back(T(row, 3*nembverts + 2*matids[0]+j, -2.0*matedge[j]));
            Dg.push_back(T(row, 3*nembverts + 2*matids[1]+j, 2.0*matedge[j]));

            Hgentry.push_back((T(3*nembverts + 2*matids[0]+j, 3*nembverts + 2*matids[0]+j, -2.0)));
            if(matids[0] < matids[1])
                Hgentry.push_back((T(3*nembverts + 2*matids[0]+j, 3*nembverts + 2*matids[1]+j, 2.0)));
            else
                Hgentry.push_back((T(3*nembverts + 2*matids[1]+j, 3*nembverts + 2*matids[0]+j, 2.0)));
            Hgentry.push_back((T(3*nembverts + 2*matids[1]+j, 3*nembverts + 2*matids[1]+j, -2.0)));
        }
        Hg.push_back(Hgentry);
        row++;
    }

    for(int i=0; i<(int)material_->getIdentifiedVerts().size(); i++)
    {
        pair<int,int> verts = material_->getIdentifiedVerts()[i];
        Vector2d pt1 = matq.segment<2>(2*verts.first);
        Vector2d pt2 = matq.segment<2>(2*verts.second);
        for(int j=0; j<2; j++)
        {
            g[row+j] = pt1[j]-pt2[j]+(j == 0 ? material_->getW() : 0.0);
            Dg.push_back(T(row+j, 3*nembverts+2*verts.first+j, 1.0));
            Dg.push_back(T(row+j, 3*nembverts+2*verts.second+j, -1.0));
        }
        vector<T> Hgentry;
        Hg.push_back(Hgentry);
        Hg.push_back(Hgentry);
        row += 2;
    }

    for(int i=0; i<(int)material_->getBoundaryVerts().size(); i++)
    {
        pair<int, double> vert = material_->getBoundaryVerts()[i];
        g[row] = matq[2*vert.first+1] - vert.second;
        Dg.push_back(T(row, 3*nembverts+2*vert.first+1, 1.0));
        vector<T> Hgentry;
        Hg.push_back(Hgentry);
        row++;
    }

    for(int i=0; i<(int)boundaries_.size(); i++)
    {
        for(int j=0; j<(int)boundaries_[i].bdryVerts.size(); j++)
        {
            Vector3d targetpt = boundaries_[i].bdryPos[j];
            int vidx = boundaries_[i].bdryVerts[j];
            for(int k=0; k<3; k++)
            {

                g[row] = embq[3*vidx+k] - targetpt[k];
                Dg.push_back(T(row, 3*vidx+k, 1.0));
                vector<T> Hgentry;
                Hg.push_back(Hgentry);
                row++;
            }
        }
    }
    assert(row == numconstraints);
    assert((int)Hg.size() == numconstraints);
}


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
    m = g.size();
    nnz_jac_g = Dg.size();

    for(int i=0; i<(int)Hg.size(); i++)
        nnz_h_lag += Hg[i].size();

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

    for(int i=0; i<m; i++)
    {
        g_l[i] = 0;
        g_u[i] = 0;
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

bool DevTLNP::eval_g(Ipopt::Index n, const Ipopt::Number *x, bool , Ipopt::Index m, Ipopt::Number *g)
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
    for(int i=0; i<m; i++)
        g[i] = Vg[i];

    return true;
}

bool DevTLNP::eval_jac_g(Ipopt::Index n, const Ipopt::Number *x, bool , Ipopt::Index, Ipopt::Index nele_jac, Ipopt::Index *iRow, Ipopt::Index *jCol, Ipopt::Number *values)
{
    if(iRow != NULL && jCol != NULL)
    {
        VectorXd g;
        vector<T> Dg;
        vector<vector<T> > Hg;
        dm_.buildConstraints(startq_, g, Dg, Hg);
        assert(nele_jac == (int)Dg.size());
        for(int i=0; i<nele_jac; i++)
        {
            iRow[i] = Dg[i].row();
            jCol[i] = Dg[i].col();
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
        assert(nele_jac == (int)Dg.size());
        for(int i=0; i<nele_jac; i++)
            values[i] = Dg[i].value();
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
        assert((int)Hg.size() == m);

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
        for(int i=0; i<(int)Hg.size(); i++)
        {
            for(int j=0; j<(int)Hg[i].size(); j++)
            {
                values[row] = lambda[i]*Hg[i][j].value();
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
