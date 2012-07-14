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
    vector<ShapeOperator> shapeOperators;
    computeShapeOperators(mesh, shapeOperators);
    vector<ShapeOperator> smoothedOperators;
    smoothShapeOperators(mesh, shapeOperators, smoothedOperators);
    shapeOperators = smoothedOperators;
    smoothShapeOperators(mesh, shapeOperators, smoothedOperators);
    curvature_.clear();
    for(int i=0; i<(int)smoothedOperators.size(); i++)
    {
        CurvatureInfo c;
        computeCurvature(smoothedOperators[i], c);
        curvature_.push_back(c);
    }
}

Frame::Frame(const Eigen::Vector3d &normal) : normal(normal)
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

void MeshCurvature::computeCurvature(const ShapeOperator &shapeOperator, CurvatureInfo &curvature)
{
    SelfAdjointEigenSolver<Matrix2d> eigensolver(shapeOperator.S);
    for(int i=0; i<2; i++)
    {
        curvature.principalCurvature[i] = eigensolver.eigenvalues()[i];
        curvature.curvatureDir[i] = eigensolver.eigenvectors().coeffRef(0,i)*shapeOperator.frame.u + eigensolver.eigenvectors().coeffRef(1,i)*shapeOperator.frame.v;
        curvature.curvatureDir[i].normalize();
    }
    if(fabs(curvature.principalCurvature[1]) < fabs(curvature.principalCurvature[0]))
    {
        swap(curvature.principalCurvature[0],curvature.principalCurvature[1]);
        swap(curvature.curvatureDir[0],curvature.curvatureDir[1]);
    }
}

void MeshCurvature::renderCurvatureDirs(Mesh &mesh)
{
    glDisable(GL_LIGHTING);
    glEnable(GL_DITHER);

    glPolygonOffset(1.0, 1.0);

    glLineWidth(2.0);
    glBegin(GL_LINES);
    assert(mesh.getMesh().n_vertices() == curvature_.size());
    for(int i=0; i<(int)curvature_.size(); i++)
    {
        {
            glColor3f(0.0, 0.0, 0.0);
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

double MeshCurvature::curvatureSpread(int vidx)
{
    assert(0 <= vidx && vidx < (int)curvature_.size());
    return fabs(fabs(curvature_[vidx].principalCurvature[0])-fabs(curvature_[vidx].principalCurvature[1]));
}

void MeshCurvature::computeShapeOperators(Mesh &mesh, vector<ShapeOperator> &operators)
{
    operators.clear();

    for(int i=0; i<(int)mesh.getMesh().n_vertices(); i++)
    {
        operators.push_back(computeShapeOperator(mesh, i));
    }
}

ShapeOperator MeshCurvature::computeShapeOperator(Mesh &mesh, int vidx)
{
    Vector3d oldnormal = mesh.vertexNormal(vidx);
    Vector3d normal = oldnormal;
    ShapeOperator result = estimateShapeOperator(mesh, vidx, normal);
    int iters = 1;
    while( (oldnormal-normal).norm() > 1e-8 )
    {
        oldnormal = normal;
        result = estimateShapeOperator(mesh, vidx, normal);
        iters++;
    }
    return result;
}

ShapeOperator MeshCurvature::estimateShapeOperator(Mesh &mesh, int vidx, Eigen::Vector3d &normal)
{
    Frame frame(normal);

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
    if(numneighbors < 6)
        return ShapeOperator(frame);

    //a u^2 + b v^2 + c uv + d u + e v + f

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
        double uval = diff.dot(frame.u);
        double vval = diff.dot(frame.v);
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

    Vector3d newnormal = (normal - coeffs[4]*frame.v - coeffs[3]*frame.u)/sqrt(1+coeffs[3]*coeffs[3]+coeffs[4]*coeffs[4]);

    double E = 1 + coeffs[3]*coeffs[3];
    double F = coeffs[3]*coeffs[4];
    double G = 1 + coeffs[4]*coeffs[4];

    double fac = normal.dot(newnormal);
    double L = 2*coeffs[0]*fac;
    double M = coeffs[2]*fac;
    double N = 2*coeffs[1]*fac;

    double det = E*G-F*F;
    if(fabs(det) < 1e-6)
        return ShapeOperator(frame);

    Matrix2d shapeOperator;
    shapeOperator << F*M-G*L, F*L-E*M, F*N-G*M, F*M-E*N;
    shapeOperator /= det;

    Matrix2d g;
    g << 1 + coeffs[4]*coeffs[4], -coeffs[3]*coeffs[4], -coeffs[3]*coeffs[4], 1 + coeffs[3]*coeffs[3];
    g /= det;

    shapeOperator = g*shapeOperator;

    assert(fabs(shapeOperator.coeff(0,1)-shapeOperator.coeff(1,0)) < 1e-6);

    normal = newnormal;

    return ShapeOperator(frame, shapeOperator);
}

Matrix2d MeshCurvature::transportShapeOperator(const ShapeOperator &source, const ShapeOperator &dest)
{
    MatrixXd sourceM(3,2);
    sourceM << source.frame.u, source.frame.v;

    MatrixXd destM(3,2);
    destM << dest.frame.u, dest.frame.v;

    Matrix2d M = sourceM.transpose()*sourceM;
    Matrix2d T = M.inverse() * sourceM.transpose()*destM;
    return T.transpose()*source.S*T;
}

void MeshCurvature::smoothShapeOperators(const Mesh &m, const vector<ShapeOperator> &oldOperators, vector<ShapeOperator> &newOperators)
{
    newOperators.clear();
    for(OMMesh::ConstVertexIter cvi = m.getMesh().vertices_begin(); cvi != m.getMesh().vertices_end(); ++cvi)
    {
        int denom = 1;
        Matrix2d S = oldOperators[cvi.handle().idx()].S;
        for(OMMesh::ConstVertexVertexIter cvvi = m.getMesh().cvv_iter(cvi.handle()); cvvi; ++cvvi)
        {
            S += transportShapeOperator(oldOperators[cvvi.handle().idx()], oldOperators[cvi.handle().idx()]);
            denom++;
        }
        S /= denom;
        newOperators.push_back(ShapeOperator(oldOperators[cvi.handle().idx()].frame, S));
    }
}

void MeshCurvature::recomputeSquaredGaussianCurvature(const Mesh &mesh, double cutoff)
{
    assert(curvature_.size() == mesh.getMesh().n_vertices());
    double totabove = 0;
    double totbelow = 0;
    for(int i=0; i<(int)curvature_.size(); i++)
    {
        double cur = gaussianCurvature(i);
        double area = mesh.areaOfInfluence(i);
        if(fabs(cur) > cutoff)
            totabove += area*cur*cur;
        else
            totbelow += area*cur*cur;
    }
    curCutoff_ = cutoff;
    belowSqGaussCurvature_ = totbelow;
    aboveSqGaussCurvature_ = totabove;
}

void MeshCurvature::totalSquaredGaussianCurvature(const Mesh &mesh, pair<double, double> &values, double cutoff)
{
    if(cutoff != curCutoff_)
    {
        recomputeSquaredGaussianCurvature(mesh, cutoff);
    }
    values.first = belowSqGaussCurvature_;
    values.second = aboveSqGaussCurvature_;
}
