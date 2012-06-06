#include "meshcurvature.h"

#include "mesh.h"
#include <set>
#include <Eigen/Cholesky>
#include <Eigen/Dense>
#include <GL/gl.h>
#include <fstream>

using namespace std;
using namespace Eigen;

CurvatureInfo::CurvatureInfo()
{
    principalCurvature[0] = 0;
    principalCurvature[1] = 0;
    curvatureDir[0].setZero();
    curvatureDir[1].setZero();
}

MeshCurvature::MeshCurvature(Mesh &mesh)
{
    computeCurvatures(mesh);
}

void MeshCurvature::computeCurvatures(Mesh &mesh)
{
    curvature_.clear();

    for(int i=0; i<(int)mesh.getMesh().n_vertices(); i++)
    {
        curvature_.push_back(computeCurvature(mesh, i));
    }
}

void MeshCurvature::computeFrame(const Eigen::Vector3d &normal, Eigen::Vector3d &u, Eigen::Vector3d &v)
{
    Vector3d seed(0,0,0);
    if(fabs(normal[0]) < 0.5)
        seed[0] = 1.0;
    else
        seed[1] = 1.0;
    u = normal.cross(seed);
    u.normalize();
    v = normal.cross(u);
}

CurvatureInfo MeshCurvature::computeCurvature(Mesh &mesh, int vidx)
{
    Vector3d oldnormal = mesh.vertexNormal(vidx);
    Vector3d normal = oldnormal;
    CurvatureInfo result;
    estimateCurvature(mesh, vidx, normal, result);
    int iters = 1;
    while( (oldnormal-normal).norm() > 1e-8 )
    {
        oldnormal = normal;
        estimateCurvature(mesh, vidx, normal, result);
        iters++;
    }
    return result;
}

void MeshCurvature::estimateCurvature(Mesh &mesh, int vidx, Vector3d &normal, CurvatureInfo &curvature)
{
    Vector3d u, v;
    computeFrame(normal, u, v);

    set<int> twoneighbors;
    twoneighbors.insert(vidx);
    for(int i=0; i<2; i++)
    {
        set<int> toadd;
        for(set<int>::iterator it = twoneighbors.begin(); it != twoneighbors.end(); ++it)
        {
            OMMesh::VertexHandle vh = mesh.getMesh().vertex_handle(*it);
            for(OMMesh::VertexVertexIter vvi = mesh.getMesh().vv_iter(vh); vvi; ++vvi)
            {
                toadd.insert(vvi.handle().idx());
            }
        }
        for(set<int>::iterator it = toadd.begin(); it != toadd.end(); ++it)
            twoneighbors.insert(*it);
    }

    int numneighbors = twoneighbors.size();
    curvature = CurvatureInfo();
    if(numneighbors < 6)
        return;

    //a u^2 + b v^2 + c uv

    MatrixXd Mat(numneighbors, 6);
    VectorXd rhs(numneighbors);
    OMMesh::Point centpt = mesh.getMesh().point(mesh.getMesh().vertex_handle(vidx));
    int row=0;
    for(set<int>::iterator it = twoneighbors.begin(); it != twoneighbors.end(); ++it)
    {
        OMMesh::Point adjpt = mesh.getMesh().point(mesh.getMesh().vertex_handle(*it));
        Vector3d diff;
        for(int i=0; i<3; i++)
            diff[i] = adjpt[i]-centpt[i];
        double uval = diff.dot(u);
        double vval = diff.dot(v);
        rhs[row] = diff.dot(normal);
        Mat.coeffRef(row,0) = uval*uval;
        Mat.coeffRef(row,1) = vval*vval;
        Mat.coeffRef(row,2) = uval*vval;
        Mat.coeffRef(row,3) = uval;
        Mat.coeffRef(row,4) = vval;
        Mat.coeffRef(row,5) = 1;
        row++;
    }
    assert(row == numneighbors);
    VectorXd nrhs = Mat.transpose()*rhs;
    MatrixXd MTM = Mat.transpose()*Mat;
    VectorXd coeffs = MTM.ldlt().solve(nrhs);

    Vector3d newnormal = (normal - coeffs[4]*v - coeffs[3]*u)/sqrt(1+coeffs[3]*coeffs[3]+coeffs[4]*coeffs[4]);

    double E = 1 + coeffs[3]*coeffs[3];
    double F = coeffs[3]*coeffs[4];
    double G = 1 + coeffs[4]*coeffs[4];

    double fac = normal.dot(newnormal);
    double L = 2*coeffs[0]*fac;
    double M = coeffs[2]*fac;
    double N = 2*coeffs[1]*fac;

    double det = E*G-F*F;
    if(fabs(det) < 1e-6)
        return;

    Matrix2d shape;
    shape << F*M-G*L, F*L-E*M, F*N-G*M, F*M-E*N;
    shape /= det;

    Matrix2d g;
    g << 1 + coeffs[4]*coeffs[4], -coeffs[3]*coeffs[4], -coeffs[3]*coeffs[4], 1 + coeffs[3]*coeffs[3];
    g /= det;

    shape = g*shape;

    assert(fabs(shape.coeff(0,1)-shape.coeff(1,0)) < 1e-6);

    SelfAdjointEigenSolver<Matrix2d> eigensolver(shape);
    for(int i=0; i<2; i++)
    {
        curvature.principalCurvature[i] = eigensolver.eigenvalues()[i];
        curvature.curvatureDir[i] = eigensolver.eigenvectors().coeffRef(0,i)*u + eigensolver.eigenvectors().coeffRef(1,i)*v;
        curvature.curvatureDir[i].normalize();
    }
    if(fabs(curvature.principalCurvature[1]) < fabs(curvature.principalCurvature[0]))
    {
        swap(curvature.principalCurvature[0],curvature.principalCurvature[1]);
        swap(curvature.curvatureDir[0],curvature.curvatureDir[1]);
    }
    normal = newnormal;
}

void MeshCurvature::renderCurvatureDirs(Mesh &mesh)
{
    glDisable(GL_LIGHTING);
    glEnable(GL_DITHER);

    glPolygonOffset(1.0, 1.0);

    glLineWidth(1.0);
    glBegin(GL_LINES);
    assert(mesh.getMesh().n_vertices() == curvature_.size());
    for(int i=0; i<(int)curvature_.size(); i++)
    {
        {
            glColor3f(1.0, 0.0, 0.0);
            OMMesh::Point pt = mesh.getMesh().point(mesh.getMesh().vertex_handle(i));

            Vector3d ept(pt[0],pt[1],pt[2]);
            double radius = mesh.shortestAdjacentEdge(i);
            Vector3d pt1 = ept - 0.5*radius*curvature_[i].curvatureDir[0];
            Vector3d pt2 = ept + 0.5*radius*curvature_[i].curvatureDir[0];

            glVertex3d(pt1[0], pt1[1], pt1[2]);
            glVertex3d(pt2[0], pt2[1], pt2[2]);
        }
    }
    glEnd();

}

double MeshCurvature::gaussianCurvature(int vidx)
{
    assert(0 <= vidx && vidx < (int)curvature_.size());
    return curvature_[vidx].principalCurvature[0]*curvature_[vidx].principalCurvature[1];
}


double MeshCurvature::meanCurvature(int vidx)
{
    assert(0 <= vidx && vidx < (int)curvature_.size());
    return 0.5*(curvature_[vidx].principalCurvature[0]+curvature_[vidx].principalCurvature[1]);
}
