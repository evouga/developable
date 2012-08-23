#include "developablemesh.h"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Geometry>
#include <set>

using namespace std;
using namespace Eigen;

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
    int boundaryverts = 0;

    for(int i=0; i<numverts; i++)
    {
        OMMesh::VertexHandle vh = mesh_.vertex_handle(i);
        if(mesh_.is_boundary(vh))
            boundaryverts++;
    }

    int numedges = mesh_.n_edges();
    int interiore = 0;
    for(int i=0; i<numedges; i++)
    {
        OMMesh::EdgeHandle eh = mesh_.edge_handle(i);
        if(!mesh_.is_boundary(eh))
            interiore++;
    }

    int numdofs = 3*numverts;

    // sum of interior vertex angles = 2pi
    // sum of boundary vertex angles = pi
    // arc length of boundaries stays fixed
    // heights of boundaries vertices are fixed
    // total area is fixed
    int numconstraints = numverts + boundaryverts + boundaries_.size() + 1;

    VectorXd g(numconstraints);
    VectorXd E(interiore);

    for(int iter = 0; iter<100; iter++)
    {
        // build the bending energy and Jacobian

        E.setZero();
        vector<T> eCoefficients;

        int row = 0;
        for(int i=0; i<(int)mesh_.n_edges(); i++)
        {
            if(mesh_.is_boundary(mesh_.edge_handle(i)))
                continue;
            OMMesh::HalfedgeHandle heh = mesh_.halfedge_handle(mesh_.edge_handle(i), 0);
            int x0i = mesh_.from_vertex_handle(heh).idx();
            Vector3d x0 = point2Vector(mesh_.point(mesh_.from_vertex_handle(heh)));
            int x1i = mesh_.to_vertex_handle(heh).idx();
            Vector3d x1 = point2Vector(mesh_.point(mesh_.to_vertex_handle(heh)));

            OMMesh::HalfedgeHandle nextheh = mesh_.next_halfedge_handle(heh);
            int x2i = mesh_.to_vertex_handle(nextheh).idx();
            Vector3d x2 = point2Vector(mesh_.point(mesh_.to_vertex_handle(nextheh)));
            assert(x0i == mesh_.to_vertex_handle(mesh_.next_halfedge_handle(nextheh)).idx());

            OMMesh::HalfedgeHandle prevheh = mesh_.next_halfedge_handle(mesh_.opposite_halfedge_handle(heh));
            int x3i = mesh_.to_vertex_handle(prevheh).idx();
            Vector3d x3 = point2Vector(mesh_.point(mesh_.to_vertex_handle(prevheh)));

            Vector3d e0 = x1-x0;
            double e0n = e0.norm();
            Vector3d e1 = x2-x0;
            double e1n = e1.norm();
            Vector3d e2 = x3-x0;
            double e2n = e2.norm();
            Vector3d e3 = x2-x1;
            double e3n = e3.norm();
            Vector3d e4 = x3-x1;
            double e4n = e4.norm();


            Vector3d n0 = e0.cross(e1);
            n0 /= n0.norm();

            Vector3d n1 = e2.cross(e0);
            n1 /= n1.norm();

            double theta = acos(n0.dot(n1));

            double cos01 = e0.dot(e1)/e0n/e1n;
            double cos02 = e0.dot(e2)/e0n/e2n;
            double cos03 = -e0.dot(e3)/e0n/e3n;
            double cos04 = -e0.dot(e4)/e0n/e4n;

            double cot01 = cos01/sqrt(1-cos01*cos01);
            double cot02 = cos02/sqrt(1-cos02*cos02);
            double cot03 = cos03/sqrt(1-cos03*cos03);
            double cot04 = cos04/sqrt(1-cos04*cos04);

            double A0 = 0.5 * (e0.cross(e1)).norm();
            double A1 = 0.5 * (e0.cross(e2)).norm();
            double A = A0 + A1;

            Vector3d t1 = n0.cross(e1);
            t1 *= e1n/t1.norm();

            Vector3d t2 = -n1.cross(e2);
            t2 *= e2n/t2.norm();

            Vector3d t3 = -n0.cross(e3);
            t3 *= e3n/t3.norm();

            Vector3d t4 = n1.cross(e4);
            t4 *= e4n/t4.norm();

            Vector3d t00 = -n0.cross(e0);
            t00 *= e0n/t00.norm();

            Vector3d t01 = n1.cross(e0);
            t01 *= e0n/t01.norm();

//            double rtE = sqrt(3 / (A) ) * e0n * theta;
            double rtE = pow(e0n, 1.0/6.0) * pow(theta, 7.0/6.0);

            Vector3d dtheta0 = -1.0/e0n * (cot03*n0 + cot04*n1);
            Vector3d dtheta1 = -1.0/e0n * (cot01*n0 + cot02*n1);
            Vector3d dtheta2 = e0n/(2*A0)*n0;
            Vector3d dtheta3 = e0n/(2*A1)*n1;

//            Vector3d dP0 = -2.0/A*e0 + e0n*e0n/(2*A*A)*(t3+t4);
//            Vector3d dP1 = 2.0/A*e0 + e0n*e0n/(2*A*A)*(t1+t2);
//            Vector3d dP2 = e0n*e0n/(2*A*A)*t00;
//            Vector3d dP3 = e0n*e0n/(2*A*A)*t01;

//            double thetacoeff = sqrt(3.0/A)*e0n;
//            double Pcoeff = sqrt(3.0)/2.0 / sqrt(e0n*e0n/A) * theta;

//            Vector3d drtE0 = thetacoeff*dtheta0 + Pcoeff * dP0;
//            Vector3d drtE1 = thetacoeff*dtheta1 + Pcoeff * dP1;
//            Vector3d drtE2 = thetacoeff*dtheta2 + Pcoeff * dP2;
//            Vector3d drtE3 = thetacoeff*dtheta3 + Pcoeff * dP3;

            Vector3d lencoeff = 1.0/6.0*pow(e0n, -5.0/6.0) * e0/e0n;
            double dthetacoeff = -7.0/6.0 * pow(e0n, 1.0/6.0) * pow(theta, 1.0/6.0);

            Vector3d drtE0 = -lencoeff + dthetacoeff*dtheta0;
            Vector3d drtE1 = lencoeff + dthetacoeff*dtheta1;
            Vector3d drtE2 = dthetacoeff*dtheta2;
            Vector3d drtE3 = dthetacoeff*dtheta3;

            E[row] = rtE;
            for(int k=0; k<3; k++)
            {
                eCoefficients.push_back(T(row, 3*x0i+k, drtE0[k]));
                eCoefficients.push_back(T(row, 3*x1i+k, drtE1[k]));
                eCoefficients.push_back(T(row, 3*x2i+k, drtE2[k]));
                eCoefficients.push_back(T(row, 3*x3i+k, drtE3[k]));
            }

            row++;
        }

        assert(row == interiore);

        // build g and Dg
        row = 0;

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
                g[row] = mesh_.point(mesh_.vertex_handle(*it))[2] - newheights[i];
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

        SparseMatrix<double> Dg(numconstraints, numdofs);
        Dg.setFromTriplets(coefficients.begin(), coefficients.end());

        SparseMatrix<double> DE(interiore, numdofs);
        DE.setFromTriplets(eCoefficients.begin(), eCoefficients.end());

        SparseMatrix<double> DETDE(DE.transpose()*DE);

        vector<T> finalcoeffs;

        for(int i=0; i<numdofs; i++)
        {
            for(int j=0; j<numdofs; j++)
            {
                if(DETDE.coeffRef(i,j) != 0)
                {
                    double coeff = DETDE.coeffRef(i,j);
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

        VectorXd rhs1 = -DE.transpose()*E;
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

        cout << "Iter " << iter << " residual " << residual << " E: " << E.squaredNorm() << " g: " << g.norm() << " lambda: " << lambda.norm() << endl;
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
