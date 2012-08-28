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

void DevelopableMesh::deformLantern(const std::vector<double> &newheights, int maxiters)
{
    assert(newheights.size() == boundaries_.size());
    for(int i=0; i<(int)boundaries_.size(); i++)
        boundaries_[i].height = newheights[i];

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
        OMMesh::Point pt = mesh_.point(mesh_.vertex_handle(i));
        for(int k=0; k<3; k++)
            q[3*i+k] = pt[k];
    }


    for(int iter = 0; iter<maxiters; iter++)
    {
        // build the bending energy and Jacobian

        SparseMatrix<double> Dg;
        VectorXd g;
        buildConstraints(q, g, Dg, newheights);
        int numconstraints = g.size();

        double f;
        VectorXd Df;
        SparseMatrix<double> Hf;
        buildObjective(q, f, Df, Hf);


//        cout << Df.transpose() << endl;
        VectorXd projV;
        projectOntoConstraints(Df, Dg, projV);
        cout << "Iter " << iter << ": proj Df: " << projV.norm() << " g: " << g.norm() << " f: " << f;

        if(projV.norm() < 1e-6 && g.norm() < 1e-6)
        {
            cout << endl << endl << endl;
            break;
        }

//        MatrixXd wtf(Hf);
//        cout << wtf << endl;

        vector<T> finalcoeffs;

        for(int i=0; i<numdofs; i++)
        {
            for(int j=0; j<numdofs; j++)
            {
                if(Hf.coeffRef(i,j) != 0)
                {
                    double coeff = Hf.coeffRef(i,j);
                    if(i == j)
                        coeff += 0.001;
                    finalcoeffs.push_back(T(i, j, coeff));
                }
            }
            for(int j=0; j<numconstraints; j++)
            {
                if(Dg.coeffRef(j,i) != 0)
                {
                    finalcoeffs.push_back(T(numdofs+j,i, Dg.coeffRef(j,i)));
                    finalcoeffs.push_back(T(i, numdofs+j, Dg.coeffRef(j,i)));
                }
            }
        }

        VectorXd rhs1 = -Df;
        VectorXd rhs2 = -g;
        VectorXd rhs(numdofs+numconstraints);
        rhs.segment(0, numdofs) = rhs1;
        rhs.segment(numdofs, numconstraints) = rhs2;

        SparseMatrix<double> M(numdofs+numconstraints, numdofs+numconstraints);
        M.setFromTriplets(finalcoeffs.begin(), finalcoeffs.end());
        MatrixXd dense(M);
        MatrixXd MTM(dense.transpose()*dense);
        VectorXd result = MTM.ldlt().solve(M.transpose()*rhs);
        VectorXd dq = result.segment(0, numdofs);
        VectorXd lambda = result.segment(numdofs, numconstraints);

        double residual = (M*result-rhs).norm();

        cout << " gradient: " << Df.dot(projV) << " lambda: " << lambda.norm() << endl;
        if(Df.dot(projV) > 0)
        {
            double alpha = lineSearch(q, -projV, lambda, newheights);
            q -= alpha*projV;
        }
        projectPositionsOntoConstraints(q, q, newheights);

    }

    for(int i=0; i<numverts; i++)
    {
        for(int k=0; k<3; k++)
            mesh_.point(mesh_.vertex_handle(i))[k] = q[3*i+k];
    }
}

Vector3d DevelopableMesh::point2Vector(OMMesh::Point pt)
{
    Vector3d result;
    for(int i=0; i<3; i++)
        result[i] = pt[i];
    return result;
}

void DevelopableMesh::buildConstraints(const VectorXd &q, VectorXd &g, SparseMatrix<double> &Dg, const vector<double> &targetHeights)
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
    int numconstraints = numverts + boundaryverts + boundaries_.size() + 1;

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

        int numedges = boundaries_[0].edges.size();
        map<int, int> vertids;
        int curid=0;

        for(int j=0; j<numedges; j++)
        {
            //double curlen = mesh_.calc_edge_length(mesh_.edge_handle(boundaries_[i].edges[j]));
            //len += curlen;
            OMMesh::HalfedgeHandle heh = mesh_.halfedge_handle(mesh_.edge_handle(boundaries_[0].edges[j]), 0);
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

        len -= boundaries_[0].arclength;
        g[row] = len.val();

        for(map<int,int>::iterator it = vertids.begin(); it != vertids.end(); ++it)
        {
            for(int k=0; k<3; k++)
                coefficients.push_back(T(row, 3*it->first+k, len.d(3*it->second+k)));
        }

        row++;
    }

    assert(row == numverts + boundaries_.size());

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
            g[row] = q[3* *it+2] - targetHeights[i];
            coefficients.push_back(T(row, 3* *it + 2, 1.0));
            row++;
        }
    }

    assert(row == numverts + boundaries_.size() + boundaryverts);

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

void DevelopableMesh::projectOntoConstraints(const Eigen::VectorXd &v, const Eigen::SparseMatrix<double> &Dg, Eigen::VectorXd &result)
{
    // Min 0.5 || v-v_0 ||^2 s.t. Dg v = 0

    VectorXd rhs = Dg*v;

    //cout << v.norm() << " " << rhs.norm() << endl;

    MatrixXd dense(Dg.transpose());
    JacobiSVD<MatrixXd> svd(dense, ComputeFullU | ComputeFullV);
    //cout << svd.singularValues().transpose() << endl;
    //VectorXd lambda = svd.solve(v);
    VectorXd tmp = svd.matrixU().transpose()*v;
    //cout << "got temp" << endl;
    cout << "Largest SV: " << svd.singularValues()[0] << "   Smallest SV: " << svd.singularValues()[rhs.size()-1] << endl;
    VectorXd scaled(rhs.size());
    for(int i=0; i<rhs.size(); i++)
    {
        double val = svd.singularValues()[i];
        if(val < 1e-6)
            scaled[i] = 0;
        else
            scaled[i] = tmp[i]/val;
    }
    //cout << "rescaled" << endl;
    VectorXd lambda = svd.matrixV()*scaled;
    //cout << "got lambda" << endl;
    /*SparseMatrix<double> temp(Dg*Dg.transpose());
    MatrixXd dense(temp);
    VectorXd lambda = dense.ldlt().solve(Dg*v);*/



    result = v - Dg.transpose()*lambda;
    double residual = (Dg*result).norm();

    //cout << lambda.norm() << " " << result.norm() << " . " << residual << " " << result.dot(v) << endl;

    //assert(residual < 1e-6);
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

        F<F<double> > theta = acos(dot(n0,n1));

        Ff += pow(e0n, 1.0/3.0) * pow(theta, 7.0/3.0);

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

double DevelopableMesh::lineSearch(const Eigen::VectorXd &q, const Eigen::VectorXd &dq, const Eigen::VectorXd &lambda, const vector<double> &targetHeights)
{
    double alpha = 1.0;
    double c1 = 1e-4;
    double c2 = 0.1;
    double fac = 0.9;

    double f;
    VectorXd Df;
    SparseMatrix<double> dummy;
    buildObjective(q, f, Df, dummy);

    cout << "Gradient " << Df.dot(dq) << endl;
    assert(Df.dot(dq) < 0);
    for(int i=0;; i++)
    {
        VectorXd newq = q + alpha*dq;
        projectPositionsOntoConstraints(newq, newq, targetHeights);
        double newf;
        VectorXd newDf;
        buildObjective(newq, newf, newDf, dummy);

        bool armijo = newf <= f + c1*alpha*Df.dot(dq);
        bool curvature = true;//newDf.dot(dq) >= c2*Df.dot(dq);

//        cout << alpha << " " << f + c1*alpha*Df.dot(dq)-newf << " " << newDf.dot(dq) - c2*Df.dot(dq) << endl;
        if(armijo && curvature)
            return alpha;
        alpha *= fac;
    }
    cout << "Line search failed!" << endl;
    return alpha;
}

void DevelopableMesh::projectPositionsOntoConstraints(const Eigen::VectorXd &q, Eigen::VectorXd &result, const vector<double> &targetHeights)
{
    VectorXd g;
    SparseMatrix<double> Dg;
    VectorXd newq = q;
    for(int i=1; i<100; i++)
    {
        buildConstraints(newq, g, Dg, targetHeights);        
        double oldg = g.norm();

        if(oldg < 1e-6)
        {
            result = newq;
            return;
        }

        MatrixXd dense(Dg);
        JacobiSVD<MatrixXd> svd(dense, ComputeFullU | ComputeFullV);

        VectorXd tmp = svd.matrixU().transpose() * -g;
        for(int i=0; i<g.size(); i++)
        {
            double val = svd.singularValues()[i];
            if(val < 1e-6)
                tmp[i] = 0;
            else
                tmp[i] /= val*val;
        }
        VectorXd lambda = svd.matrixU() * tmp;

        //SparseMatrix<double> DgDgT = Dg*Dg.transpose();
        //MatrixXd dense(DgDgT);
        //BiCGSTAB<SparseMatrix<double> > solver(DgDgT);
        //VectorXd lambda = solver.solve(-g);
        //VectorXd lambda = dense.ldlt().solve(-g);
        newq += Dg.transpose()*lambda;
    }
    result = newq;
}
