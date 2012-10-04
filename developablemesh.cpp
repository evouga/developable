#include "developablemesh.h"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Geometry>
#include <set>

#include <fadbad.h>
#include <fadiff.h>
#include "coin/IpIpoptApplication.hpp"
#include "devnlp.h"
#include <QGLWidget>

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

    for(int i=0; i<(int)material_->getMesh().n_faces(); i++)
    {
        OMMesh::FaceHandle fh = material_->getMesh().face_handle(i);
        for(OMMesh::FaceHalfedgeIter fhi = material_->getMesh().fh_iter(fh); fhi; ++fhi)
        {
            OMMesh::HalfedgeHandle heh = fhi.handle();
            OMMesh::HalfedgeHandle nextheh = material_->getMesh().next_halfedge_handle(heh);
            OMMesh::VertexHandle centv = material_->getMesh().to_vertex_handle(heh);
            OMMesh::VertexHandle v1 = material_->getMesh().from_vertex_handle(heh);
            OMMesh::VertexHandle v2 = material_->getMesh().to_vertex_handle(nextheh);

            OMMesh::Point e0 = material_->getMesh().point(v1) - material_->getMesh().point(centv);
            OMMesh::Point e1 = material_->getMesh().point(v2) - material_->getMesh().point(centv);
            Vector3d ve0 = point2Vector(e0);
            Vector3d ve1 = point2Vector(e1);
            Vector3d cr = ve0.cross(ve1);
            assert(cr[2] < 0);
        }
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
//    app.Options()->SetStringValue("derivative_test", "second-order");
    app.Options()->SetStringValue("check_derivatives_for_naninf", "yes");
    app.Initialize();
    app.OptimizeTNLP(mynlp);

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
    checkInversion(q);

/*
    f = 0;
    Df.resize(q.size());
    Df.setZero();
    Hf.clear();

    double r = 0.53033;
    int n = 4;
    int m = 3;
    //double h = 0.6;
    //double alpha = -2.54968;
    double h = 0.8;
    double alpha = -2.8174;
    double a = 2*r*sin(PI/n);
    double b = sqrt(4*r*r*sin(alpha/2.0)*sin(alpha/2.0) + h*h/m/m);
    double c = sqrt(4*r*r*sin(PI/n + alpha/2.0)*sin(PI/n + alpha/2.0) + h*h/m/m);


    double s = 0.5*(a+b+c);
    double area = sqrt(s*(s-a)*(s-b)*(s-c));
    double dpy = 2*area/a;

    double W = n*a;
    double H = 2*n*m*area/W;

    assert(fabs(H-dpy*m) < 1e-6);
    int vidx = 0;
    Vector3d targetpt;

    for(int i=0; i<=m; i++)
    {
        targetpt[2] = h*i/m;

        for(int j=0; j<n; j++)
        {
            targetpt[0] = r*cos(2*PI*(j/double(n)) + i*alpha);
            targetpt[1] = r*sin(2*PI*(j/double(n)) + i*alpha);

            Vector3d vpos = q.segment<3>(3*vidx);
            for(int k=0; k<3; k++)
            {
                f += (vpos[k]-targetpt[k])*(vpos[k]-targetpt[k]);
                Df[3*vidx+k] = 2.0*(vpos[k]-targetpt[k]);
                Hf.push_back(T(3*vidx+k, 3*vidx+k,2.0));
            }
            vidx++;
        }
    }
*/

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

        /*F<F<double> > Lpart = 0;
        if(L.val().val() < 1e-8)
            Lpart = 0.5 * (norm(e02) + norm(e03));
        else
            Lpart = L*L/(A0+A1);*/
        F<F<double> > Lpart = pow(L, 1.0/3.0);
        //F<F<double> > anglepart = acos(costheta)*acos(costheta);
        F<F<double> > anglepart = pow(acos(costheta), 7.0/3.0);
        F<F<double> > stencilf = Lpart*anglepart;
        f += stencilf.val().val();
        assert(!isnan(f));
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

void DevelopableMesh::buildInversionConstraints(const Eigen::VectorXd &q, Eigen::VectorXd &h, std::vector<T> &Dh, std::vector<std::vector<T> > &Hh)
{
    int f = 0;
    for(OMMesh::FaceIter fi = material_->getMesh().faces_begin(); fi != material_->getMesh().faces_end(); ++fi)
        for(OMMesh::FaceHalfedgeIter fhi = material_->getMesh().fh_iter(fi); fhi; ++fhi)
            f++;
    h.resize(f);
    h.setZero();
    Dh.clear();
    Hh.clear();

    int offset = 3*mesh_.n_vertices();

    int row = 0;
    for(OMMesh::FaceIter fi = material_->getMesh().faces_begin(); fi != material_->getMesh().faces_end(); ++fi)
    {
        for(OMMesh::FaceHalfedgeIter fhi = material_->getMesh().fh_iter(fi); fhi; ++fhi)
        {
            OMMesh::HalfedgeHandle heh1 = fhi.handle();
            OMMesh::HalfedgeHandle heh2 = material_->getMesh().next_halfedge_handle(heh1);
            int centv = material_->getMesh().to_vertex_handle(heh1).idx();
            assert(centv == material_->getMesh().from_vertex_handle(heh2).idx());
            int v1 = material_->getMesh().from_vertex_handle(heh1).idx();
            int v2 = material_->getMesh().to_vertex_handle(heh2).idx();
            Vector2d e1 = q.segment<2>(offset+2*v1) - q.segment<2>(offset+2*centv);
            Vector2d e2 = q.segment<2>(offset+2*v2) - q.segment<2>(offset+2*centv);
            Vector2d emid = q.segment<2>(offset+2*v2) - q.segment<2>(offset+2*v1);
            h[row] = -e1[0]*e2[1] + e1[1]*e2[0];
            Dh.push_back(T(row, offset+2*v1, -e2[1]));
            Dh.push_back(T(row, offset+2*v1+1, e2[0]));
            Dh.push_back(T(row, offset+2*v2, e1[1]));
            Dh.push_back(T(row, offset+2*v2+1, -e1[0]));
            Dh.push_back(T(row, offset+2*centv, emid[1]));
            Dh.push_back(T(row, offset+2*centv+1, -emid[0]));
            vector<T> Hhpart;
            if(2*v1 <= 2*v2+1)
                Hhpart.push_back(T(offset+2*v1, offset+2*v2+1, -1.0));
            else
                Hhpart.push_back(T(offset+2*v2+1, offset+2*v1, -1.0));

            if(2*v1+1 <= 2*v2)
                Hhpart.push_back(T(offset+2*v1+1, offset+2*v2, 1.0));
            else
                Hhpart.push_back(T(offset+2*v2, offset+2*v1+1, 1.0));

            if(2*v1 <= 2*centv+1)
                Hhpart.push_back(T(offset+2*v1, offset+2*centv+1, 1.0));
            else
                Hhpart.push_back(T(offset+2*centv+1, offset+2*v1, 1.0));

            if(2*v1+1 <= 2*centv)
                Hhpart.push_back(T(offset+2*v1+1, offset+2*centv, -1.0));
            else
                Hhpart.push_back(T(offset+2*centv, offset+2*v1+1, -1.0));

            if(2*v2 <= 2*centv+1)
                Hhpart.push_back(T(offset+2*v2, offset+2*centv+1, -1.0));
            else
                Hhpart.push_back(T(offset+2*centv+1, offset+2*v2, -1.0));

            if(2*v2+1 <= 2*centv)
                Hhpart.push_back(T(offset+2*v2+1, offset+2*centv, 1.0));
            else
                Hhpart.push_back(T(offset+2*centv, offset+2*v2+1, 1.0));
            Hh.push_back(Hhpart);
            row++;
        }
    }
    assert(row == (int)h.size());
    assert(row == (int)Hh.size());
}

void DevelopableMesh::checkInversion(const Eigen::VectorXd &q)
{
    int n = mesh_.n_vertices();
    int matoffset = 3*n;
    /*for(int i=0; i<(int)material_->getMesh().n_faces(); i++)
    {
        OMMesh::FaceHandle fh = material_->getMesh().face_handle(i);
        for(OMMesh::FaceHalfedgeIter fhi = material_->getMesh().fh_iter(fh); fhi; ++fhi)
        {
            OMMesh::HalfedgeHandle heh1 = fhi.handle();
            OMMesh::HalfedgeHandle heh2 = material_->getMesh().next_halfedge_handle(heh1);
            int centv = material_->getMesh().to_vertex_handle(heh1).idx();
            int v1 = material_->getMesh().from_vertex_handle(heh1).idx();
            int v2 = material_->getMesh().to_vertex_handle(heh2).idx();
            Vector2d e1 = q.segment<2>(matoffset + 2*v1) - q.segment<2>(matoffset+2*centv);
            Vector2d e2 = q.segment<2>(matoffset + 2*v2) - q.segment<2>(matoffset+2*centv);
            double cr = e1[0]*e2[1]-e1[1]*e2[0];
            if(cr > 0)
                cout << "Triangle Inversion" << endl;
        }
    }*/

    double minlength = std::numeric_limits<double>::infinity();
    for(int i=0; i<(int)material_->getMesh().n_edges(); i++)
    {
        OMMesh::EdgeHandle eh = material_->getMesh().edge_handle(i);
        OMMesh::HalfedgeHandle heh = material_->getMesh().halfedge_handle(eh,0);
        int v1 = material_->getMesh().to_vertex_handle(heh).idx();
        int v2 = material_->getMesh().from_vertex_handle(heh).idx();
        Vector2d vv1 = q.segment<2>(matoffset+2*v1);
        Vector2d vv2 = q.segment<2>(matoffset+2*v2);
        double len = (vv1-vv2).norm();
        minlength = min(minlength, len);
    }
    if(minlength < 1e-6)
        cout << "Edge collapse" << endl;
}

Vector2d DevelopableMesh::materialCenter()
{
    Vector2d result(0,0);
    for(int i=0; i<(int)material_->getMesh().n_vertices(); i++)
    {
        OMMesh::VertexHandle vh = material_->getMesh().vertex_handle(i);
        OMMesh::Point pt = material_->getMesh().point(vh);
        result[0] += pt[0];
        result[1] += pt[1];
    }
    result /= material_->getMesh().n_vertices();
    return result;
}

double DevelopableMesh::materialRadius()
{
    Vector2d center = materialCenter();
    double maxradius = 0;
    for(int i=0; i<(int)material_->getMesh().n_vertices(); i++)
    {
        OMMesh::VertexHandle vh = material_->getMesh().vertex_handle(i);
        OMMesh::Point pt = material_->getMesh().point(vh);
        Vector2d vpt(pt[0],pt[1]);
        maxradius = max(maxradius, (center-vpt).norm());
    }
    return maxradius;
}

void DevelopableMesh::renderMaterial()
{
    glBegin(GL_LINES);
    for(int i=0; i<(int)material_->getMesh().n_edges(); i++)
    {
        OMMesh::EdgeHandle eh = material_->getMesh().edge_handle(i);
        OMMesh::HalfedgeHandle heh = material_->getMesh().halfedge_handle(eh, 0);
        OMMesh::Point pt1 = material_->getMesh().point(material_->getMesh().to_vertex_handle(heh));
        OMMesh::Point pt2 = material_->getMesh().point(material_->getMesh().from_vertex_handle(heh));
        glVertex2d(pt1[0], pt1[1]);
        glVertex2d(pt2[0], pt2[1]);
    }
    glEnd();
}
