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

void DevelopableMesh::deformLantern(const std::vector<double> &newheights)
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

    for(int iter = 0; iter<100; iter++)
    {
        // build the bending energy and Jacobian

        SparseMatrix<double> Dg;
        VectorXd g;
        buildConstraints(g, Dg, newheights);
        int numconstraints = g.size();

        double f;
        VectorXd Df;
        SparseMatrix<double> Hf;
        buildObjective(f, Df, Hf);


        VectorXd projV;
        projectOntoConstraints(Df, Dg, projV);
        cout << projV.transpose() << endl;

        vector<T> finalcoeffs;

        for(int i=0; i<numdofs; i++)
        {
            for(int j=0; j<numdofs; j++)
            {
                if(Hf.coeffRef(i,j) != 0)
                {
                    double coeff = Hf.coeffRef(i,j);
                    if(i == j)
                        coeff += 0.0;
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
        rhs.setZero();
        for(int i=0; i<numdofs; i++)
            rhs[i] = rhs1[i];
        for(int i=0; i<numconstraints; i++)
            rhs[i+numdofs] = rhs2[i];

        SparseMatrix<double> M(numdofs+numconstraints, numdofs+numconstraints);
        M.setFromTriplets(finalcoeffs.begin(), finalcoeffs.end());
        MatrixXd dense(M);
        MatrixXd MTM(dense.transpose()*dense);
        VectorXd result = MTM.ldlt().solve(M.transpose()*rhs);
        VectorXd dq(numdofs);
        for(int i=0; i<numdofs; i++)
            dq[i] = result[i];
        VectorXd lambda(numconstraints);
        for(int i=0; i<numconstraints; i++)
            lambda[i] = result[i+numdofs];

        double residual = (M*result-rhs).norm();

        cout << "Iter " << iter << " residual " << residual << " E: " << f << " g: " << g.norm() << " lambda: " << lambda.norm() << endl;
        // Apply delta
        for(int i=0; i<numverts; i++)
        {
            for(int k=0; k<3; k++)
                mesh_.point(mesh_.vertex_handle(i))[k] += dq[3*i+k];
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

void DevelopableMesh::buildConstraints(VectorXd &g, SparseMatrix<double> &Dg, const vector<double> &targetHeights)
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

        double deficit = 0;

        int valence = mesh_.valence(vh);

        vector<int> vidxs;
        vector<Vector3d> grads;
        for(int j=0; j<=valence; j++)
        {
            vidxs.push_back(0);
            grads.push_back(Vector3d(0,0,0));
        }

        vidxs[valence] = i;

        OMMesh::VertexOHalfedgeIter voh = mesh_.voh_iter(vh);
        for(int j=0; j<valence; j++)
        {
            OMMesh::VertexHandle tov = mesh_.to_vertex_handle(voh);
            vidxs[j] = tov.idx();
            OMMesh::Point e1pt = mesh_.point(tov) - mesh_.point(vh);
            Vector3d e1(e1pt[0], e1pt[1], e1pt[2]);

            ++voh;
            if(!voh)
                voh = mesh_.voh_iter(vh);

            // Ignore boundary "faces"
            if(mesh_.is_boundary(voh))
                continue;


            OMMesh::Point e2pt = mesh_.point(mesh_.to_vertex_handle(voh)) - mesh_.point(vh);
            Vector3d e2(e2pt[0], e2pt[1], e2pt[2]);
            double e1norm = e1.norm();
            double e2norm = e2.norm();
            double theta = acos( e1.dot(e2) / e1norm / e2norm );
            deficit += theta;

            Vector3d de1 = (e1.dot(e2)/e1norm * e1/e1norm - e2) / sqrt( e1norm*e1norm * e2norm*e2norm - e1.dot(e2)*e1.dot(e2));
            Vector3d de2 = (e1.dot(e2)/e2norm * e2/e2norm - e1) / sqrt( e1norm*e1norm * e2norm*e2norm - e1.dot(e2)*e1.dot(e2));
            grads[j] += de1;
            grads[(j+1)%valence] += de2;
            grads[valence] += -de1 - de2;
        }
        deficit -= (bdry ? PI : 2*PI);
        g[row] = deficit;

        for(int j=0; j<=valence; j++)
        {
            for(int k=0; k<3; k++)
            {
                coefficients.push_back(T(row, 3*vidxs[j]+k, grads[j][k]));
            }
        }
        row++;
    }
    assert(row == numverts);

    // boundary arc lengths are constant
    for(int i=0; i<(int)boundaries_.size(); i++)
    {
        double len = 0;

        map<int, Vector3d> grads;
        int numedges = boundaries_[i].edges.size();

        for(int j=0; j<numedges; j++)
        {
            double curlen = mesh_.calc_edge_length(mesh_.edge_handle(boundaries_[i].edges[j]));
            len += curlen;
            OMMesh::HalfedgeHandle heh = mesh_.halfedge_handle(mesh_.edge_handle(boundaries_[i].edges[j]), 0);
            int idx1 = mesh_.to_vertex_handle(heh).idx();
            int idx2 = mesh_.from_vertex_handle(heh).idx();
            OMMesh::Point pt1pt = mesh_.point(mesh_.vertex_handle(idx1));
            OMMesh::Point pt2pt = mesh_.point(mesh_.vertex_handle(idx2));
            Vector3d pt1(pt1pt[0], pt1pt[1], pt1pt[2]);
            Vector3d pt2(pt2pt[0], pt2pt[1], pt2pt[2]);

            if(grads.find(idx1) == grads.end())
                grads[idx1] = (pt1-pt2)/curlen;
            else
                grads[idx1] += (pt1-pt2)/curlen;
            if(grads.find(idx2) == grads.end())
                grads[idx2] = -(pt1-pt2)/curlen;
            else
                grads[idx2] += -(pt1-pt2)/curlen;
        }

        len -= boundaries_[i].arclength;
        g[row] = len;

        for(map<int, Vector3d>::iterator it = grads.begin(); it != grads.end(); ++it)
        {
            for(int k=0; k<3; k++)
                coefficients.push_back(T(row, 3*it->first+k, it->second[k]));
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
            g[row] = mesh_.point(mesh_.vertex_handle(*it))[2] - targetHeights[i];
            coefficients.push_back(T(row, 3* *it + 2, 1.0));
            row++;
        }
    }

    assert(row == numverts + (int)boundaries_.size() + boundaryverts);

    // total area is constant
    double area = 0;

    vector<Vector3d> gradA;
    for(int i=0; i<numverts; i++)
        gradA.push_back(Vector3d(0,0,0));

    for(int i=0; i<numverts; i++)
    {
        OMMesh::VertexHandle vh = mesh_.vertex_handle(i);

        int valence = mesh_.valence(vh);
        OMMesh::VertexOHalfedgeIter voh = mesh_.voh_iter(vh);
        for(int j=0; j<valence; j++)
        {
            OMMesh::VertexHandle tov = mesh_.to_vertex_handle(voh);
            int idx1 = tov.idx();
            OMMesh::Point e1pt = mesh_.point(tov) - mesh_.point(vh);
            Vector3d e1(e1pt[0], e1pt[1], e1pt[2]);

            ++voh;
            if(!voh)
                voh = mesh_.voh_iter(vh);

            if(mesh_.is_boundary(voh))
                continue;

            tov = mesh_.to_vertex_handle(voh);
            int idx2 = tov.idx();
            OMMesh::Point e2pt = mesh_.point(tov) - mesh_.point(vh);
            Vector3d e2(e2pt[0], e2pt[1], e2pt[2]);
            area += 0.5/3.0 * (e1.cross(e2)).norm();

            Vector3d de1 = -0.5/3.0 * (e1.cross(e2)).cross(e2)/(e1.cross(e2)).norm();
            Vector3d de2 = 0.5/3.0 * (e1.cross(e2)).cross(e1)/(e1.cross(e2)).norm();
            gradA[idx1] += de1;
            gradA[idx2] += de2;
            gradA[i] += -de1-de2;
        }
    }

    area -= surfacearea_;
    g[row] = area;

    for(int i=0; i<numverts; i++)
    {
        for(int k=0; k<3; k++)
            coefficients.push_back(T(row, 3*i+k, gradA[i][k]));
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
    SparseMatrix<double> DgDgT(Dg*Dg.transpose());
    SimplicialLDLT<SparseMatrix<double> > solver(DgDgT);
    VectorXd lambda = solver.solve(rhs);
    double residual = (DgDgT*lambda-rhs).norm();
    cout << "Projection residual " << residual << endl;

    result = v - Dg.transpose()*lambda;
}

void DevelopableMesh::buildObjective(double &f, Eigen::VectorXd &Df, Eigen::SparseMatrix<double> &Hf)
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

        F<double> Ff = 0;

        OMMesh::HalfedgeHandle heh = mesh_.halfedge_handle(mesh_.edge_handle(i), 0);
        int inds[4];
        F<double> x[4][3];

        inds[0] = mesh_.from_vertex_handle(heh).idx();
        for(int k=0; k<3;k++)
        {
            x[0][k] = mesh_.point(mesh_.from_vertex_handle(heh))[k];
            x[0][k].diff(3*0+k, 12);
        }

        inds[1] = mesh_.to_vertex_handle(heh).idx();
        for(int k=0; k<3; k++)
        {
            x[1][k] = mesh_.point(mesh_.to_vertex_handle(heh))[k];
            x[1][k].diff(3*1+k, 12);
        }

        OMMesh::HalfedgeHandle nextheh = mesh_.next_halfedge_handle(heh);
        inds[2] = mesh_.to_vertex_handle(nextheh).idx();
        for(int k=0; k<3; k++)
        {
            x[2][k] = mesh_.point(mesh_.to_vertex_handle(nextheh))[k];
            x[2][k].diff(3*2+k, 12);
        }

        OMMesh::HalfedgeHandle prevheh = mesh_.next_halfedge_handle(mesh_.opposite_halfedge_handle(heh));
        inds[3] = mesh_.to_vertex_handle(prevheh).idx();
        for(int k=0; k<3; k++)
        {
            x[3][k] = mesh_.point(mesh_.to_vertex_handle(prevheh))[k];
            x[3][k].diff(3*3+k, 12);
        }

        F<double> e0[3], e1[3], e2[3];
        for(int k=0; k<3; k++)
        {
            e0[k] = x[1][k]-x[0][k];
            e1[k] = x[2][k]-x[0][k];
            e2[k] = x[3][k]-x[0][k];
        }

        F<double> e0n = norm(e0);

        F<double> n0[3], n1[3];

        cross(e0, e1, n0);
        cross(e2, e0, n1);

        normalize(n0);
        normalize(n1);

        F<double> theta = acos(dot(n0,n1));

        Ff += pow(e0n, 1.0/3.0) * pow(theta, 7.0/3.0);

        f += Ff.val();

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

}

F<double> DevelopableMesh::norm(fadbad::F<double> *v)
{
    F<double> result = 0;
    for(int i=0; i<3; i++)
        result += v[i]*v[i];
    return sqrt(result);
}

void DevelopableMesh::cross(fadbad::F<double> *v1, fadbad::F<double> *v2, fadbad::F<double> *result)
{
    result[0] = v1[1]*v2[2] - v1[2]*v2[1];
    result[1] = v1[2]*v2[0] - v1[0]*v2[2];
    result[2] = v1[0]*v2[1] - v1[1]*v2[0];
}

void DevelopableMesh::normalize(fadbad::F<double> *v)
{
    F<double> n = norm(v);
    for(int i=0; i<3; i++)
        v[i] /= n;
}

F<double> DevelopableMesh::dot(fadbad::F<double> *v1, fadbad::F<double> *v2)
{
    F<double> result = 0;
    for(int i=0; i<3; i++)
        result += v1[i]*v2[i];
    return result;
}
