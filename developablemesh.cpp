#include "developablemesh.h"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Geometry>
#include <set>
#include <fadbad.h>
#include <fadiff.h>

using namespace std;
using namespace Eigen;
using namespace fadbad;


const double PI = 3.1415926535898;

typedef Eigen::Triplet<double> T;

DevelopableMesh::DevelopableMesh()
{

}

bool DevelopableMesh::loadMesh(const string &filename)
{
    if(!Mesh::loadMesh(filename))
        return false;
    identifyBoundaries();
    calculateSurfaceArea();
    return true;
}

void DevelopableMesh::buildSchwarzLantern(double r, double h, int n, int m)
{
    mesh_ = OMMesh();

    for(int i=0; i<m; i++)
    {
        double z = h*i/(m-1);
        for(int j=0; j<n; j++)
        {
            double x = r*cos(2*PI*(j/double(n) + (i%2)/double(2*n)));
            double y = r*sin(2*PI*(j/double(n) + (i%2)/double(2*n)));

            OMMesh::Point newpt(x,y,z);

            mesh_.add_vertex(newpt);
        }
    }
    for(int i=0; i<m-1; i++)
    {
        for(int j=0; j<n; j++)
        {
            int fidx1 = i*n+j;
            int fidx2 = i*n + ((j+1) % n);
            int fidx3 = (i+1)*n+ ((j + i%2) % n);
            vector<OMMesh::VertexHandle> newface;
            newface.push_back(mesh_.vertex_handle(fidx1));
            newface.push_back(mesh_.vertex_handle(fidx2));
            newface.push_back(mesh_.vertex_handle(fidx3));
            mesh_.add_face(newface);

            fidx2 = fidx3;
            newface[1] = mesh_.vertex_handle(fidx2);
            fidx3 = (i+1)*n + ((n + j-1 + i%2) % n);
            newface[2] = mesh_.vertex_handle(fidx3);
            mesh_.add_face(newface);

        }
    }

    identifyBoundaries();
    calculateSurfaceArea();
}

void DevelopableMesh::identifyBoundaries()
{
    boundaries_.clear();
    int numedges = mesh_.n_edges();
    bool *visited = new bool[numedges];

    for(int i=0; i<numedges; i++)
        visited[i] = false;


    for(int i=0; i<numedges; i++)
    {
        if(visited[i])
            continue;

        if(mesh_.is_boundary(mesh_.edge_handle(i)))
        {
            BoundaryCurve bc;
            bc.arclength = 0;
            OMMesh::EdgeHandle eh = mesh_.edge_handle(i);
            OMMesh::HalfedgeHandle heh = mesh_.halfedge_handle(eh,0);
            if(!mesh_.is_boundary(heh))
                heh = mesh_.opposite_halfedge_handle(heh);
            bc.height = mesh_.point(mesh_.from_vertex_handle(heh))[2];
            bc.targetheight = bc.height;
            while(!visited[eh.idx()])
            {
                visited[eh.idx()] = true;
                bc.edges.push_back(eh.idx());
                bc.arclength += mesh_.calc_edge_length(eh);

                heh = mesh_.next_halfedge_handle(heh);
                eh = mesh_.edge_handle(heh);

                assert(mesh_.point(mesh_.from_vertex_handle(heh))[2] == bc.height);
            }

            boundaries_.push_back(bc);

        }

    }

    delete[] visited;
}

void DevelopableMesh::calculateSurfaceArea()
{
    surfacearea_ = 0;
    for(int i=0; i< (int)mesh_.n_vertices(); i++)
    {
        OMMesh::VertexHandle vh = mesh_.vertex_handle(i);

        int valence = mesh_.valence(vh);
        OMMesh::VertexOHalfedgeIter voh = mesh_.voh_iter(vh);
        for(int i=0; i<valence; i++)
        {
            OMMesh::Point e1pt = mesh_.point(mesh_.to_vertex_handle(voh)) - mesh_.point(vh);
            Vector3d e1(e1pt[0], e1pt[1], e1pt[2]);

            ++voh;
            if(!voh)
                voh = mesh_.voh_iter(vh);

            if(mesh_.is_boundary(voh))
                continue;

            OMMesh::Point e2pt = mesh_.point(mesh_.to_vertex_handle(voh)) - mesh_.point(vh);
            Vector3d e2(e2pt[0], e2pt[1], e2pt[2]);
            surfacearea_ += 0.5/3.0 * (e1.cross(e2)).norm();
        }
    }
}

void DevelopableMesh::getBoundaryHeights(std::vector<double> &heights)
{
    heights.clear();
    for(int i=0; i<(int)boundaries_.size(); i++)
        heights.push_back(boundaries_[i].height);
}

void DevelopableMesh::deformLantern(int maxiters)
{
    int numverts = mesh_.n_vertices();

    int numedges = mesh_.n_edges();
    int interiore = 0;
    for(int i=0; i<numedges; i++)
    {
        OMMesh::EdgeHandle eh = mesh_.edge_handle(i);
        if(!mesh_.is_boundary(eh))
            interiore++;
    }

    int numdofs = 3*numverts;    

    VectorXd q(numdofs);
    for(int i=0; i<numverts; i++)
    {
        OMMesh::Point &pt = mesh_.point(mesh_.vertex_handle(i));
        for(int k=0; k<3; k++)
            q[3*i+k] = pt[k];
    }

    VectorXd g;
    SparseMatrix<double> Dg;
    double f;
    VectorXd Df;
    SparseMatrix<double> Hf;

    buildObjective(q, f, Df, Hf);
    buildConstraints(q, g, Dg);

    for(int iter=0; iter<20; iter++)
    {
        VectorXd startq = q;
        VectorXd dq(numdofs);
        dq.setZero();

        vector<T> bigcoeffs;

        int numconstraints = g.size();

        SparseMatrix<double> M(numdofs+numconstraints, numdofs+numconstraints);
        M.setZero();
        for(int i=0; i<numdofs; i++)
            for(int j=0; j<numdofs; j++)
            {
                double coeff = Hf.coeffRef(i,j);
                if(coeff != 0.0)
                    bigcoeffs.push_back(T(i,j,coeff));
            }
        for(int i=0; i<numconstraints; i++)
        {
            for(int j=0; j<numdofs; j++)
            {
                double coeff = Dg.coeffRef(i,j);
                if(coeff != 0.0)
                {
                    bigcoeffs.push_back(T(numdofs+i,j,coeff));
                    bigcoeffs.push_back(T(j, numdofs+i, coeff));
                }
            }
        }

        M.setFromTriplets(bigcoeffs.begin(), bigcoeffs.end());

        BiCGSTAB<SparseMatrix<double> > solver(M);
        MatrixXd dense(M);
        JacobiSVD<MatrixXd> svd(dense, ComputeFullU | ComputeFullV);

        VectorXd rhs(numdofs+numconstraints);
        rhs.setZero();
        rhs.segment(0, numdofs) = -Df;


        VectorXd tmp = svd.matrixU().transpose()*rhs;
        for(int i=0; i<numdofs+numconstraints; i++)
        {
            double div = svd.singularValues()[i];
            if(div < 1e-4)
                tmp[i] = 0;
            else
                tmp[i] /= div;
        }

        VectorXd result = svd.matrixV() * tmp;

        double residual = (M*result-rhs).norm();

        cout << " residual " << residual;

        dq = result.segment(0, numdofs);
        q += dq;
        projectOntoConstraints(q, q);
        cout << "Step " << (q-startq).norm() << endl;
    }

    return;
}

Vector3d DevelopableMesh::point2Vector(OMMesh::Point pt)
{
    Vector3d result;
    for(int i=0; i<3; i++)
        result[i] = pt[i];
    return result;
}

void DevelopableMesh::buildConstraints(const VectorXd &q, VectorXd &g, SparseMatrix<double> &Dg)
{
    int numverts = mesh_.n_vertices();
    int numdofs = 3*numverts;
    int boundaryverts = 0;

    for(int i=0; i<numverts; i++)
    {
        OMMesh::VertexHandle vh = mesh_.vertex_handle(i);
        if(mesh_.is_boundary(vh))
            boundaryverts++;
    }

    // sum of interior vertex angles = 2pi
    // sum of boundary vertex angles = pi
    // arc length of boundaries stays fixed
    // heights of boundaries vertices are fixed
    // total area is fixed
    int numconstraints = numverts + 3*boundaryverts + boundaries_.size() + 1;

    // build g and Dg
    int row = 0;
    g.resize(numconstraints);
    g.setZero();
    std::vector<T> coefficients;


    // interior angles sum to 2 PI, boundary angles sum to PI
    for(int i=0; i<numverts; i++)
    {
        OMMesh::VertexHandle vh = mesh_.vertex_handle(i);
        bool bdry = mesh_.is_boundary(vh);

        F<double> deficit = 0;

        int valence = mesh_.valence(vh);

        vector<int> vidxs;
        for(int j=0; j<=valence; j++)
        {
            vidxs.push_back(0);
        }

        vidxs[valence] = i;
        F<double> centpt[3];
        for(int k=0; k<3; k++)
        {
            centpt[k] = q[3*vh.idx()+k];
            centpt[k].diff(3*valence+k, 3*(valence+1));
        }

        OMMesh::VertexOHalfedgeIter voh = mesh_.voh_iter(vh);
        for(int j=0; j<valence; j++)
        {
            OMMesh::VertexHandle tov = mesh_.to_vertex_handle(voh);
            vidxs[j] = tov.idx();
            F<double> tovpt[3];
            for(int k=0; k<3; k++)
            {
                tovpt[k] = q[3*tov.idx()+k];
                tovpt[k].diff(3*j+k, 3*(valence+1));
            }
            F<double> e1[3];
            for(int k=0; k<3; k++)
            {
                e1[k] = tovpt[k] - centpt[k];
            }

            ++voh;
            if(!voh)
                voh = mesh_.voh_iter(vh);

            // Ignore boundary "faces"
            if(mesh_.is_boundary(voh))
                continue;


            //OMMesh::Point e2pt = mesh_.point(mesh_.to_vertex_handle(voh)) - mesh_.point(vh);
            F<double> tovpt2[3];
            for(int k=0; k<3; k++)
            {
                tovpt2[k] = q[3*mesh_.to_vertex_handle(voh).idx()+k];
                tovpt2[k].diff( 3*((j+1)%valence) + k, 3*(valence+1) );
            }
            F<double> e2[3];
            for(int k=0; k<3; k++)
            {
                e2[k] = tovpt2[k]-centpt[k];
            }
            F<double> e1norm = norm(e1);
            F<double> e2norm = norm(e2);
            F<double> theta = acos(dot(e1,e2) / e1norm / e2norm );
            deficit += theta;
        }
        deficit -= (bdry ? PI : 2*PI);
        g[row] = deficit.val();

        for(int j=0; j<=valence; j++)
        {
            for(int k=0; k<3; k++)
            {
                coefficients.push_back(T(row, 3*vidxs[j]+k, deficit.d(3*j+k)));
            }
        }
        row++;
    }
    assert(row == numverts);

    // boundary arc length is constant
    for(int i=0; i<(int)boundaries_.size(); i++)
    {
        F<double> len = 0;

        int numedges = boundaries_[i].edges.size();
        map<int, int> vertids;
        int curid=0;

        for(int j=0; j<numedges; j++)
        {
            //double curlen = mesh_.calc_edge_length(mesh_.edge_handle(boundaries_[i].edges[j]));
            //len += curlen;
            OMMesh::HalfedgeHandle heh = mesh_.halfedge_handle(mesh_.edge_handle(boundaries_[i].edges[j]), 0);
            int idx1 = mesh_.to_vertex_handle(heh).idx();
            int idx2 = mesh_.from_vertex_handle(heh).idx();
            if(vertids.find(idx1) == vertids.end())
                vertids[idx1] = curid++;
            if(vertids.find(idx2) == vertids.end())
                vertids[idx2] = curid++;

            F<double> pt1[3];
            F<double> pt2[3];
            for(int k=0; k<3; k++)
            {
                pt1[k] = q[3*idx1+k];
                pt1[k].diff(3*vertids[idx1]+k, 3*numedges);
                pt2[k] = q[3*idx2+k];
                pt2[k].diff(3*vertids[idx2]+k, 3*numedges);
            }

            F<double> edge[3];
            for(int k=0; k<3; k++)
            {
                edge[k] = pt2[k] - pt1[k];
            }

            len += norm(edge);

        }

        len -= boundaries_[i].arclength;
        g[row] = len.val();

        for(map<int,int>::iterator it = vertids.begin(); it != vertids.end(); ++it)
        {
            for(int k=0; k<3; k++)
                coefficients.push_back(T(row, 3*it->first+k, len.d(3*it->second+k)));
        }

        row++;
    }

    assert(row == numverts + (int)boundaries_.size());

    // boundary vertices have fixed height
    for(int i=0; i<(int)boundaries_.size(); i++)
    {
        set<int> bdverts;
        for(int j=0; j<(int)boundaries_[i].edges.size(); j++)
        {
            OMMesh::HalfedgeHandle heh = mesh_.halfedge_handle(mesh_.edge_handle(boundaries_[i].edges[j]),0);
            bdverts.insert(mesh_.from_vertex_handle(heh).idx());
            bdverts.insert(mesh_.to_vertex_handle(heh).idx());
        }
        for(set<int>::iterator it = bdverts.begin(); it != bdverts.end(); ++it)
        {
            for(int k=0;k<3;k++)
            {
                OMMesh::Point pt = mesh_.point(mesh_.vertex_handle(*it));
                g[row] = q[3* *it+k] - pt[k];//boundaries_[i].targetheight;
                coefficients.push_back(T(row, 3* *it + k, 1.0));
                row++;
            }
        }
    }

    assert(row == numverts + (int)boundaries_.size() + 3*boundaryverts);

    // total area is constant
    F<double> area = 0;

    for(int i=0; i<numverts; i++)
    {
        OMMesh::VertexHandle vh = mesh_.vertex_handle(i);
        F<double> centpt[3];
        for(int k=0; k<3; k++)
        {
            centpt[k] = q[3*vh.idx()+k];
            centpt[k].diff(3*i+k, 3*numverts);
        }
        int valence = mesh_.valence(vh);
        OMMesh::VertexOHalfedgeIter voh = mesh_.voh_iter(vh);
        for(int j=0; j<valence; j++)
        {
            OMMesh::VertexHandle tov = mesh_.to_vertex_handle(voh);
            int idx1 = tov.idx();
            F<double> tovpt[3];
            for(int k=0; k<3; k++)
            {
                tovpt[k] = q[3*tov.idx() + k];
                tovpt[k].diff(3*idx1+k, 3*numverts);
            }
            F<double> e1[3];
            for(int k=0; k<3; k++)
                e1[k] = tovpt[k]-centpt[k];

            ++voh;
            if(!voh)
                voh = mesh_.voh_iter(vh);

            if(mesh_.is_boundary(voh))
                continue;

            tov = mesh_.to_vertex_handle(voh);
            int idx2 = tov.idx();
            F<double> tovpt2[3];
            for(int k=0; k<3; k++)
            {
                tovpt2[k] = q[3*tov.idx()+k];
                tovpt2[k].diff(3*idx2+k, 3*numverts);
            }
            F<double> e2[3];
            for(int k=0; k<3; k++)
                e2[k] = tovpt2[k]-centpt[k];

            F<double> cr[3];
            cross(e1,e2,cr);
            area += 0.5/3.0 * norm(cr);
        }
    }

    area -= surfacearea_;
    g[row] = area.val();

    for(int i=0; i<numverts; i++)
    {
        for(int k=0; k<3; k++)
            coefficients.push_back(T(row, 3*i+k, area.d(3*i+k)));
    }
    row++;

    assert(row == numconstraints);

    Dg.resize(numconstraints, numdofs);
    Dg.setZero();

    Dg.setFromTriplets(coefficients.begin(), coefficients.end());

}

void DevelopableMesh::buildObjective(const VectorXd &q, double &f, Eigen::VectorXd &Df, Eigen::SparseMatrix<double> &Hf)
{
    int numverts = mesh_.n_vertices();
    int numdofs = 3*numverts;

    f = 0;
    Df.resize(numdofs);
    Df.setZero();
    map<pair<int, int>, double> nonzeros;

    for(int i=0; i<(int)mesh_.n_edges(); i++)
    {
        if(mesh_.is_boundary(mesh_.edge_handle(i)))
            continue;

        F<F<double> > Ff = 0;

        OMMesh::HalfedgeHandle heh = mesh_.halfedge_handle(mesh_.edge_handle(i), 0);
        int inds[4];
        F<F<double> > x[4][3];
        F<double> xx[4][3];

        inds[0] = mesh_.from_vertex_handle(heh).idx();
        inds[1] = mesh_.to_vertex_handle(heh).idx();
        OMMesh::HalfedgeHandle nextheh = mesh_.next_halfedge_handle(heh);
        inds[2] = mesh_.to_vertex_handle(nextheh).idx();
        OMMesh::HalfedgeHandle prevheh = mesh_.next_halfedge_handle(mesh_.opposite_halfedge_handle(heh));
        inds[3] = mesh_.to_vertex_handle(prevheh).idx();

        for(int j=0; j<4; j++)
            for(int k=0; k<3; k++)
            {
                xx[j][k] = q[3*inds[j]+k];
                xx[j][k].diff(3*j+k, 12);
                x[j][k] = xx[j][k];
                x[j][k].diff(3*j+k, 12);
            }

        F<F<double> > e0[3], e1[3], e2[3];
        for(int k=0; k<3; k++)
        {
            e0[k] = x[1][k]-x[0][k];
            e1[k] = x[2][k]-x[0][k];
            e2[k] = x[3][k]-x[0][k];
        }

        F< F<double> > e0n = norm(e0);

        F<F<double> > n0[3], n1[3];

        cross(e0, e1, n0);
        cross(e2, e0, n1);

        normalize(n0);
        normalize(n1);

        F<F<double> > theta = acos(dot(n0,n1)-1e-10);

        Ff = pow(e0n, 1.0/3.0) * pow(theta, 7.0/3.0);

        //cout << Ff.val().val() << endl;

        // NAN alert
        if(theta.val().val() == 0.0)
            continue;

        f += Ff.val().val();

        for(int j=0; j<4; j++)
        {
            for(int k=0; k<3; k++)
            {
                F<double> ideriv = Ff.d(3*j+k);
                Df[3*inds[j]+k] += ideriv.val();

                for(int l=0; l<4; l++)
                {
                    for(int m=0; m<3; m++)
                        nonzeros[pair<int,int>(3*inds[j]+k,3*inds[l]+m)] += ideriv.d(3*l+m);
                }
            }
        }
    }

    vector<T> eCoefficients;
    for(map<pair<int, int>, double>::iterator it = nonzeros.begin(); it != nonzeros.end(); ++it)
    {
        eCoefficients.push_back(T(it->first.first, it->first.second, it->second));
    }


    Hf.resize(numdofs, numdofs);
    Hf.setZero();
    Hf.setFromTriplets(eCoefficients.begin(), eCoefficients.end());

    // sanity check
    /*for(int i=0; i<numdofs; i++)
    {
        for(int j=0; j<numdofs; j++)
        {
            assert(fabs(Hf.coeffRef(i,j) - Hf.coeffRef(j,i)) < 1e-6);
        }
    }*/
}

double DevelopableMesh::lineSearch(const Eigen::VectorXd &q, const Eigen::VectorXd &dq, double mu)
{
    double alpha = 1.0;
    double c1 = 1e-4;
    //double c2 = 0.1;
    double fac = 0.9;

    double f;
    VectorXd g;
    VectorXd Df;
    SparseMatrix<double> Dg;
    SparseMatrix<double> dummy;
    buildObjective(q, f, Df, dummy);
    buildConstraints(q, g, Dg);

    double obj =f + 0.5*mu*g.dot(g);
    VectorXd Dobj = Df + mu*Dg.transpose()*g;

    assert(Dobj.dot(dq) < 0);


    for(int i=0;; i++)
    {
        double newf;
        VectorXd newDf, newg;
        SparseMatrix<double> newDg;

        VectorXd newq = q + alpha*dq;
        buildObjective(newq, newf, newDf, dummy);
        buildConstraints(newq, newg, newDg);
        double newobj = newf + 0.5*mu*newg.dot(newg);
        VectorXd newDobj = newDf + mu * newDg.transpose()*newg;

        bool armijo = newobj < obj + c1*alpha*Dobj.dot(dq);
        //bool curvature = newDobj.dot(dq) >= c2*Dobj.dot(dq);
        if(armijo)
        {
            return alpha;
        }
        alpha *= fac;
    }
    cout << "Line search failed!" << endl;
    return 0;
}

bool DevelopableMesh::canCollapseEdge(int edgeidx)
{
    OMMesh::EdgeHandle eh = mesh_.edge_handle(edgeidx);
    bool bdry = false;
    set<int> ignore;

    for(int i=0; i<2; i++)
    {
        OMMesh::HalfedgeHandle heh = mesh_.halfedge_handle(eh, i);
        OMMesh::HalfedgeHandle flip = mesh_.opposite_halfedge_handle(heh);
        if(mesh_.is_boundary(heh))
        {
            if(mesh_.is_boundary(flip))
                return false;
            if(mesh_.is_boundary(mesh_.next_halfedge_handle(flip)))
                return false;
            if(mesh_.is_boundary(mesh_.prev_halfedge_handle(flip)))
                return false;

            bdry = true;
        }
        else
            ignore.insert(mesh_.to_vertex_handle(mesh_.next_halfedge_handle(heh)).idx());
    }

    if(!bdry)
    {
        for(int i=0; i<2; i++)
        {
            OMMesh::HalfedgeHandle heh = mesh_.halfedge_handle(eh, i);
            OMMesh::HalfedgeHandle flip = mesh_.opposite_halfedge_handle(heh);


            if(mesh_.is_boundary(mesh_.next_halfedge_handle(heh)))
            {
                if(mesh_.is_boundary(mesh_.next_halfedge_handle(flip))
                        || mesh_.is_boundary(mesh_.prev_halfedge_handle(heh)))
                {
                    return false;
                }
            }
        }
    }

    OMMesh::HalfedgeHandle heh = mesh_.halfedge_handle(eh,0);
    OMMesh::VertexHandle v1 = mesh_.from_vertex_handle(heh);
    set<int> adj1;
    for(OMMesh::VertexVertexIter vvi = mesh_.vv_iter(v1); vvi; ++vvi)
    {
        adj1.insert(vvi.handle().idx());
    }
    OMMesh::VertexHandle v2 = mesh_.to_vertex_handle(heh);
    for(OMMesh::VertexVertexIter vvi = mesh_.vv_iter(v2); vvi; ++vvi)
    {
        int idx = vvi.handle().idx();
        if(adj1.count(idx) > 0 && ignore.count(idx) == 0)
        {
            return false;
        }
    }
    if(bdry)
        return false;
    return true;
}

void DevelopableMesh::collapseEdge(int edgeidx)
{
    assert(canCollapseEdge(edgeidx));
    OMMesh::EdgeHandle eh = mesh_.edge_handle(edgeidx);
    OMMesh::HalfedgeHandle heh = mesh_.halfedge_handle(eh, 0);
    OMMesh::VertexHandle tov = mesh_.to_vertex_handle(heh);
    OMMesh::VertexHandle frv = mesh_.from_vertex_handle(heh);
    OMMesh::Point midpt = (mesh_.point(tov) + mesh_.point(frv))*0.5;

    OMMesh newmesh;

    map<int, int> old2new;
    int newidx = 0;
    int mididx = -1;
    for(int i=0; i<(int)mesh_.n_vertices(); i++)
    {
        OMMesh::Point newpt;
        if(i == frv.idx() || i == tov.idx())
        {
            if(mididx == -1)
            {
                newpt = midpt;
                mididx = newidx;
            }
            else
            {
                old2new[i] = mididx;
                continue;
            }
        }
        else
            newpt = mesh_.point(mesh_.vertex_handle(i));
        newmesh.add_vertex(newpt);
        old2new[i] = newidx++;
    }
    assert(newidx == (int)mesh_.n_vertices()-1);

    for(int i=0; i<(int)mesh_.n_faces(); i++)
    {
        OMMesh::FaceHandle fh = mesh_.face_handle(i);
        if(fh == mesh_.face_handle(heh)
                || fh == mesh_.face_handle(mesh_.opposite_halfedge_handle(heh)))
            continue;

        vector<OMMesh::VertexHandle> newface;
        for(OMMesh::FaceVertexIter fvi = mesh_.fv_iter(fh); fvi; ++fvi)
            newface.push_back(newmesh.vertex_handle(old2new[fvi.handle().idx()]));
        newmesh.add_face(newface);
    }

    // fix edges
    for(int i=0; i<(int)boundaries_.size(); i++)
    {
        vector<int> newedges;
        for(int j=0; j<(int)boundaries_[i].edges.size(); j++)
        {
            int oldedgeidx = boundaries_[i].edges[j];
            OMMesh::EdgeHandle oldeh = mesh_.edge_handle(oldedgeidx);
            if(oldeh.idx() == eh.idx())
                continue;
            OMMesh::HalfedgeHandle oldhe = mesh_.halfedge_handle(oldeh,0);
            assert(mesh_.is_boundary(oldeh));
            int v1 = mesh_.to_vertex_handle(oldhe).idx();
            int v2 = mesh_.from_vertex_handle(oldhe).idx();
            int newv1 = old2new[v1];
            int newv2 = old2new[v2];

            bool found = false;

            for(int k=0; k<(int)newmesh.n_halfedges(); k++)
            {
                OMMesh::HalfedgeHandle newhe = newmesh.halfedge_handle(k);
                if(newv1 == newmesh.to_vertex_handle(newhe).idx()
                        && newv2 == newmesh.from_vertex_handle(newhe).idx())
                {
                    assert(!found);
                    found = true;
                    OMMesh::EdgeHandle newbdry = newmesh.edge_handle(newhe);
                    assert(newmesh.is_boundary(newbdry));
                    newedges.push_back(newbdry.idx());
                }

            }
            assert(found);
        }
        boundaries_[i].edges = newedges;
    }
    mesh_ = newmesh;
}

bool DevelopableMesh::collapseShortEdges()
{
    bool atleastone = false;
    bool collapsed = false;
    do
    {
        collapsed = false;

        vector<double> widths;
        computeCreaseWidths(widths);

        vector<double> vertwidths;
        for(int i=0; i<(int)mesh_.n_vertices(); i++)
        {
            OMMesh::VertexHandle vh = mesh_.vertex_handle(i);
            double maxwidth = 0;
            for(OMMesh::VertexEdgeIter vei = mesh_.ve_iter(vh); vei; ++vei)
            {
                maxwidth = std::max(maxwidth, widths[vei.handle().idx()]);
            }
            vertwidths.push_back(maxwidth);
        }

        for(int i=0; i<(int)mesh_.n_edges(); i++)
        {
            OMMesh::EdgeHandle eh = mesh_.edge_handle(i);
            OMMesh::HalfedgeHandle heh = mesh_.halfedge_handle(eh, 0);
            double maxlen = vertwidths[mesh_.from_vertex_handle(heh).idx()];
            maxlen = std::max(maxlen, vertwidths[mesh_.to_vertex_handle(heh).idx()]);

            double len = mesh_.calc_edge_length(eh);
            if(len < 1e-1)//maxlen)
            {
                if(canCollapseEdge(i))
                {
                    collapseEdge(i);
                    collapsed=true;
                    atleastone = true;
                    break;
                }
            }
        }
    } while(collapsed);

    return atleastone;
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

void DevelopableMesh::computeCreaseWidths(std::vector<double> &widths)
{
    widths.clear();

    for(int i=0; i<(int)mesh_.n_edges(); i++)
    {
        if(mesh_.is_boundary(mesh_.edge_handle(i)))
        {
            widths.push_back(0);
        }
        else
        {
            double theta = turningAngle(i);
            double L = mesh_.calc_edge_length(mesh_.edge_handle(i));
            widths.push_back(theta*pow(L, 2.0/3.0));
        }
    }
}

double DevelopableMesh::turningAngle(int edge)
{
    OMMesh::EdgeHandle eh = mesh_.edge_handle(edge);
    OMMesh::HalfedgeHandle heh = mesh_.halfedge_handle(eh,0);
    if(mesh_.is_boundary(eh))
    {
        cout << "boundary" << endl;
        return 0;
    }
    OMMesh::VertexHandle v0h = mesh_.from_vertex_handle(heh);
    OMMesh::VertexHandle v00h = mesh_.to_vertex_handle(heh);
    OMMesh::VertexHandle v1h = mesh_.from_vertex_handle(mesh_.prev_halfedge_handle(heh));
    OMMesh::VertexHandle v2h = mesh_.to_vertex_handle(mesh_.next_halfedge_handle(mesh_.opposite_halfedge_handle(heh)));

    Vector3d e0 = point2Vector(mesh_.point(v00h) - mesh_.point(v0h));
    Vector3d e1 = point2Vector(mesh_.point(v1h) - mesh_.point(v0h));
    Vector3d e2 = point2Vector(mesh_.point(v2h) - mesh_.point(v0h));

    Vector3d n0 = e0.cross(e1);
    n0 /= n0.norm();

    Vector3d n1 = e2.cross(e0);
    n1 /= n1.norm();
    return acos(n0.dot(n1));
}

void DevelopableMesh::projectVectorOntoConstraints(const VectorXd &q, const VectorXd &v, VectorXd &newv)
{
    // min 0.5||v - v0||^2 s.t. \grad G v = 0
    // v - v0 + \grad G^T lambda = 0
    // \grad G v = 0
    //
    // \grad G \grad G^T lambda = \grad G v0
    // v = v0 - \grad G^T lambda

    VectorXd g;
    SparseMatrix<double> Dg;
    newv = v;
    buildConstraints(q, g, Dg);
    SparseMatrix<double> M(Dg*Dg.transpose());
    VectorXd rhs = Dg * v;
    MatrixXd dense(M);

    JacobiSVD<MatrixXd> svd(dense, ComputeFullU | ComputeFullV);
    VectorXd tmp = svd.matrixU().transpose()*rhs;
    for(int i=0; i<g.size(); i++)
    {
        double val = svd.singularValues()[i];
        if(val < 1e-6)
            tmp[i] = 0;
        else
            tmp[i] = tmp[i]/val;
    }
    VectorXd lambda = svd.matrixV()*tmp;
    VectorXd dv = -Dg.transpose()*lambda;
    newv += dv;
}

void DevelopableMesh::projectOntoConstraints(const VectorXd &q, VectorXd &newq)
{
    VectorXd g;
    SparseMatrix<double> Dg;
    newq = q;
    while(true)
    {
        buildConstraints(newq, g, Dg);
        if(g.norm() < 1e-5)
        {
            cout << "Final violation " << g.norm() << endl;
            return;
        }
        SparseMatrix<double> M(Dg*Dg.transpose());
        VectorXd rhs = g;
        MatrixXd dense(M);

        JacobiSVD<MatrixXd> svd(dense, ComputeFullU | ComputeFullV);
        VectorXd tmp = svd.matrixU().transpose()*rhs;
        for(int i=0; i<g.size(); i++)
        {
            double val = svd.singularValues()[i];
            if(val < 1e-6)
                tmp[i] = 0;
            else
                tmp[i] = tmp[i]/val;
        }
        VectorXd lambda = svd.matrixV()*tmp;
        //residual = (M*lambda+rhs).norm();
        VectorXd dq = -Dg.transpose()*lambda;
        double relchange = dq.norm() / newq.norm();
        newq += dq;
        if( relchange < 1e-6)
        {
            cout << "Final violation " << g.norm() << endl;
            return;
        }
    }
}
