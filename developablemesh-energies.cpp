#include "developablemesh.h"
#include <Eigen/Core>
#include <FADBAD++/fadbad.h>
#include <FADBAD++/fadiff.h>

using namespace Eigen;
using namespace fadbad;
using namespace std;


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
        bool onboundary = false;//(mesh_.is_boundary(mesh_.from_vertex_handle(heh)) || mesh_.is_boundary(mesh_.to_vertex_handle(heh)));
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
//        normalize(n0);
//        normalize(n1);
//        F<F<double> > costheta = dot(n0,n1);
//        F<F<double> > n0crossn1[3];
//        cross(n0, n1, n0crossn1);
//        F<F<double> > sintheta = norm(n0crossn1);
//        F<F<double> > theta = atan2(sintheta,costheta);

        F<F<double> > Lpart = 0;

        if(onboundary)
            Lpart = L*L/(A0+A1);
        else
            Lpart = pow(L, 1.0/3.0);

//        F<F<double> > anglepart = 0;
//        if(onboundary)
//            anglepart = theta*theta;
//        else
//            anglepart = pow(theta, 7.0/3.0);
        double anglepart = 1;
        F<F<double> > stencilf = Lpart*anglepart;
        f += stencilf.val().val();
        assert(!isnan(f));
        for(int i=0; i<4; i++)
        {
            for(int j=0; j<3; j++)
            {
                double dfval = stencilf.d(3*i+j).val();
                Df[3*vertidx[i]+j] += dfval;

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
        numconstraints++;
    }

//    // top and bottom verts have y values fixed
//    numconstraints += material_->getBoundaryVerts().size();

//    // bottom verts have x values fixed
//    int numbottom = 0;
//    for(int i=0; i<(int)material_->getBoundaryVerts().size(); i++)
//    {
//        if(material_->getBoundaryVerts()[i].onBottom)
//            numbottom++;
//    }
//    numconstraints += numbottom;

    int numemb = 0;
    // top and bottom of embedded cylinder have values fixed
    for(int i=0; i<(int)boundaries_.size(); i++)
        numemb += 3*boundaries_[i].bdryVerts.size();
    numconstraints += numemb;

    // check on the following

    g.resize(numconstraints);
    g.setZero();
    Dg.clear();
    Hg.clear();

    int row=0;

    for(int i=0; i<(int)mesh_.n_edges(); i++)
    {
        OMMesh::EdgeHandle eh = mesh_.edge_handle(i);
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
        Vector3d offset = material_->getOffset(mheh);

        Vector3d embedge = embq.segment<3>(3*embids[1]) - embq.segment<3>(3*embids[0]);
        double len1 = embedge.squaredNorm();
        Vector2d matedge = matq.segment<2>(2*matids[1]) - matq.segment<2>(2*matids[0]);
        matedge += offset.segment<2>(0);
        double len2 = matedge.squaredNorm();
        g[row] = (len1-len2);

        vector<T> Hgentry;

        for(int j=0; j<3; j++)
        {
            Dg.push_back(T(row, 3*embids[0]+j,-2.0*embedge[j]));
            Dg.push_back(T(row, 3*embids[1]+j,2.0*embedge[j]));

            Hgentry.push_back((T(3*embids[0]+j, 3*embids[0]+j, 2.0)));
            if(embids[0] < embids[1])
                Hgentry.push_back((T(3*embids[0]+j, 3*embids[1]+j, -2.0)));
            else
                Hgentry.push_back((T(3*embids[1]+j, 3*embids[0]+j, -2.0)));
            Hgentry.push_back((T(3*embids[1]+j, 3*embids[1]+j, 2.0)));
        }
        for(int j=0; j<2; j++)
        {
            Dg.push_back(T(row, 3*nembverts + 2*matids[0]+j, 2.0*matedge[j]));
            Dg.push_back(T(row, 3*nembverts + 2*matids[1]+j, -2.0*matedge[j]));

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

void DevelopableMesh::radiusOverlapConstraint(const Eigen::VectorXd &q, OMMesh::HalfedgeHandle heh1, double &g, std::vector<std::pair<int, double> > &Dg, std::vector<T> &Hg)
{
//    const double k = 0.001;

//    int offset = 3*mesh_.n_vertices();

//    OMMesh::HalfedgeHandle heh2 = material_->getMesh().next_halfedge_handle(heh1);
//    int vinds[3];
//    vinds[0] = material_->getMesh().to_vertex_handle(heh1).idx();
//    assert(vinds[0] == material_->getMesh().from_vertex_handle(heh2).idx());

//    vinds[1] = material_->getMesh().from_vertex_handle(heh1).idx();
//    vinds[2] = material_->getMesh().to_vertex_handle(heh2).idx();

//    F<double> vx[3][2];
//    F<F<double> > vxx[3][2];
//    for(int i=0; i<3; i++)
//    {
//        for(int j=0; j<2; j++)
//        {
//            vx[i][j] = q[offset+2*vinds[i]+j];
//            vx[i][j].diff(2*i+j, 6);
//            vxx[i][j] = vx[i][j];
//            vxx[i][j].diff(2*i+j, 6);
//        }
//    }

//    F<F<double> > e1[2];
//    F<F<double> > e2[2];
//    diff2D(vxx[1],vxx[0],e1);
//    diff2D(vxx[2],vxx[0],e2);
//    for(int i=0; i<2; i++)
//    {
//        e1[i] -= material_->getOffset(heh1)[i];
//        e2[i] += material_->getOffset(heh2)[i];
//    }

//    F<F<double> > Ldot = dot2D(e1, e2);

//    F<F<double> > L1 = norm2D(e1);
//    F<F<double> > L2 = norm2D(e2);

//    F<F<double> > num = L1*L2*(pow(L1, 1.0/3.0) * pow(L2, 1.0/3.0) - 4*k*k);
//    F<F<double> > denom = sqrt( (pow(L1, 2.0/3.0)+4*k*k)*(pow(L2,2.0/3.0)+4*k*k));
//    F<F<double> > val = num/denom - Ldot;

//    g = val.val().val();

//    for(int i=0; i<3; i++)
//    {
//        for(int j=0; j<2; j++)
//        {
//            Dg.push_back(std::pair<int, double>(offset+2*vinds[i]+j, val.deriv(2*i+j).val()));
//            for(int k=0; k<3; k++)
//            {
//                for(int l=0; l<2; l++)
//                {
//                    if(2*vinds[i]+j <= 2*vinds[k]+l)
//                    {
//                        Hg.push_back(T(offset+2*vinds[i]+j, offset+2*vinds[k]+l, val.deriv(2*i+j).deriv(2*k+l)));
//                    }
//                }
//            }
//        }
//    }
}

void DevelopableMesh::buildInversionConstraints(const Eigen::VectorXd &q, Eigen::VectorXd &h, std::vector<T> &Dh, std::vector<std::vector<T> > &Hh)
{
//    int f = 0;
//    for(OMMesh::FaceIter fi = material_->getMesh().faces_begin(); fi != material_->getMesh().faces_end(); ++fi)
//        for(OMMesh::FaceHalfedgeIter fhi = material_->getMesh().fh_iter(fi); fhi; ++fhi)
//            f += 2;
//    h.resize(f);
//    h.setZero();
//    Dh.clear();
//    Hh.clear();

//    int offset = 3*mesh_.n_vertices();

//    int row = 0;
//    for(OMMesh::FaceIter fi = material_->getMesh().faces_begin(); fi != material_->getMesh().faces_end(); ++fi)
//    {
//        for(OMMesh::FaceHalfedgeIter fhi = material_->getMesh().fh_iter(fi); fhi; ++fhi)
//        {
//            OMMesh::HalfedgeHandle heh1 = fhi.handle();
//            OMMesh::HalfedgeHandle heh2 = material_->getMesh().next_halfedge_handle(heh1);
//            int centv = material_->getMesh().to_vertex_handle(heh1).idx();
//            assert(centv == material_->getMesh().from_vertex_handle(heh2).idx());
//            int v1 = material_->getMesh().from_vertex_handle(heh1).idx();
//            int v2 = material_->getMesh().to_vertex_handle(heh2).idx();
//            Vector2d e1 = q.segment<2>(offset+2*v1) - q.segment<2>(offset+2*centv);
//            e1 -= material_->getOffset(heh1).segment<2>(0);
//            Vector2d e2 = q.segment<2>(offset+2*v2) - q.segment<2>(offset+2*centv);
//            e2 += material_->getOffset(heh2).segment<2>(0);
//            Vector2d emid = e2 - e1;
//            h[row] = -e1[0]*e2[1] + e1[1]*e2[0];
//            Dh.push_back(T(row, offset+2*v1, -e2[1]));
//            Dh.push_back(T(row, offset+2*v1+1, e2[0]));
//            Dh.push_back(T(row, offset+2*v2, e1[1]));
//            Dh.push_back(T(row, offset+2*v2+1, -e1[0]));
//            Dh.push_back(T(row, offset+2*centv, emid[1]));
//            Dh.push_back(T(row, offset+2*centv+1, -emid[0]));
//            vector<T> Hhpart;
//            if(2*v1 <= 2*v2+1)
//                Hhpart.push_back(T(offset+2*v1, offset+2*v2+1, -1.0));
//            else
//                Hhpart.push_back(T(offset+2*v2+1, offset+2*v1, -1.0));

//            if(2*v1+1 <= 2*v2)
//                Hhpart.push_back(T(offset+2*v1+1, offset+2*v2, 1.0));
//            else
//                Hhpart.push_back(T(offset+2*v2, offset+2*v1+1, 1.0));

//            if(2*v1 <= 2*centv+1)
//                Hhpart.push_back(T(offset+2*v1, offset+2*centv+1, 1.0));
//            else
//                Hhpart.push_back(T(offset+2*centv+1, offset+2*v1, 1.0));

//            if(2*v1+1 <= 2*centv)
//                Hhpart.push_back(T(offset+2*v1+1, offset+2*centv, -1.0));
//            else
//                Hhpart.push_back(T(offset+2*centv, offset+2*v1+1, -1.0));

//            if(2*v2 <= 2*centv+1)
//                Hhpart.push_back(T(offset+2*v2, offset+2*centv+1, -1.0));
//            else
//                Hhpart.push_back(T(offset+2*centv+1, offset+2*v2, -1.0));

//            if(2*v2+1 <= 2*centv)
//                Hhpart.push_back(T(offset+2*v2+1, offset+2*centv, 1.0));
//            else
//                Hhpart.push_back(T(offset+2*centv, offset+2*v2+1, 1.0));
//            Hh.push_back(Hhpart);
//            row++;
//        }
//    }

//    for(OMMesh::FaceIter fi = material_->getMesh().faces_begin(); fi != material_->getMesh().faces_end(); ++fi)
//    {
//        for(OMMesh::FaceHalfedgeIter fhi = material_->getMesh().fh_iter(fi); fhi; ++fhi)
//        {
//            OMMesh::HalfedgeHandle heh1 = fhi.handle();
//            vector<pair<int, double> > Dg;
//            vector<T> Hgpart;
//            radiusOverlapConstraint(q, heh1, h[row], Dg, Hgpart);
//            for(vector<pair<int, double> >::iterator it = Dg.begin(); it != Dg.end(); ++it)
//            {
//                Dh.push_back(T(row, it->first, it->second));
//            }
//            Hh.push_back(Hgpart);
//            row++;
//        }
//    }
//    assert(row == (int)h.size());
//    assert(row == (int)Hh.size());
}

