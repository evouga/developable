#include "developablemesh.h"
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

void DevelopableMesh::shapeMatch(const std::vector<Vector3d> &sourcePts, const std::vector<Vector3d> &targetPts, std::vector<Vector3d> &resultPts)
{
    int numpts = (int)sourcePts.size();
    MatrixXd A(numpts,3);
    MatrixXd B(numpts,3);
    for(int i=0; i<numpts; i++)
    {
        for(int j=0; j<3; j++)
        {
            A.coeffRef(i,j) = sourcePts[i][j];
            B.coeffRef(i,j) = targetPts[i][j];
        }
    }

    Matrix3d M = A.transpose()*B;
    double det = M.determinant();

    JacobiSVD<Matrix3d> svd(M, ComputeFullU | ComputeFullV );
    Matrix3d newS = Matrix3d::Identity();
    if(det < 0)
        newS.coeffRef(2,2) = -1;
    Matrix3d R = svd.matrixU()*newS*svd.matrixV().transpose();
    resultPts.clear();
    for(int i=0; i<numpts; i++)
        resultPts.push_back(R*sourcePts[i]);
}

OMMesh *DevelopableMesh::warpMaterialToEmbedded(const Eigen::VectorXd &q)
{
    OMMesh *result = new OMMesh();
    int nembverts = mesh_.n_vertices();
    int nmatverts = material_->getMesh().n_vertices();
    assert(nembverts == nmatverts);
    assert(q.size() == 3*nembverts + 2*nmatverts);
    VectorXd embq = q.segment(0,3*nembverts);
    VectorXd matq = q.segment(3*nembverts, 2*nmatverts);

    for(OMMesh::FaceIter fi = mesh_.faces_begin(); fi != mesh_.faces_end(); ++fi)
    {
        OMMesh::FaceHandle fh = fi.handle();
        Vector3d totoffset(0,0,0);
        vector<Vector3d> matverts;
        vector<Vector3d> embverts;

        for(OMMesh::FaceHalfedgeIter fhi = mesh_.fh_iter(fh); fhi; ++fhi)
        {
            int vertid = mesh_.from_vertex_handle(fhi).idx();
            embverts.push_back(embq.segment<3>(3*vertid));
            Vector3d matvert3d;
            matvert3d.segment<2>(0) = matq.segment<2>(2*vertid);
            matvert3d[2] = 0;
            matverts.push_back(matvert3d+totoffset);

            int targetid = mesh_.to_vertex_handle(fhi).idx();
            OMMesh::HalfedgeHandle matheh = material_->getMesh().halfedge_handle(material_->findHalfedge(vertid, targetid));
            totoffset += material_->getOffset(matheh);
        }

        vector<Vector3d> newtri;
        shapeMatch(matverts, embverts, newtri);

        vector<OMMesh::VertexHandle> newvhs;
        for(vector<Vector3d>::iterator it = newtri.begin(); it != newtri.end(); ++it)
        {
            OMMesh::Point newpt;
            for(int j=0; j<3; j++)
                newpt[j] = (*it)[j];
            newvhs.push_back(result->add_vertex(newpt));
        }
        result->add_face(newvhs);
    }
    return result;
}
