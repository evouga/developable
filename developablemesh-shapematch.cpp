#include "developablemesh.h"
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

void DevelopableMesh::shapeMatch(const std::vector<Vector3d> &sourcePts, const std::vector<Vector3d> &targetPts, std::vector<Vector3d> &resultPts)
{
    int numpts = (int)sourcePts.size();
    MatrixXd A(numpts,3);
    MatrixXd B(numpts,3);

    Vector3d sourceCentroid;
    sourceCentroid.setZero();
    for(int i=0; i<numpts; i++)
        sourceCentroid += sourcePts[i];
    sourceCentroid /= numpts;

    Vector3d targetCentroid;
    targetCentroid.setZero();
    for(int i=0; i<numpts; i++)
        targetCentroid += targetPts[i];
    targetCentroid /= numpts;

    for(int i=0; i<numpts; i++)
    {
        for(int j=0; j<3; j++)
        {
            A.coeffRef(i,j) = sourcePts[i][j] - sourceCentroid[j];
            B.coeffRef(i,j) = targetPts[i][j] - targetCentroid[j];
        }
    }

    Matrix3d M = B.transpose()*A;
    double det = M.determinant();

    JacobiSVD<Matrix3d> svd(M, ComputeFullU | ComputeFullV );
    Matrix3d newS = Matrix3d::Identity();
    if(det < 0)
        newS.coeffRef(2,2) = -1;
    Matrix3d R = svd.matrixU()*newS*svd.matrixV().transpose();
    resultPts.clear();
    for(int i=0; i<numpts; i++)
        resultPts.push_back(targetCentroid + R*(sourcePts[i]-sourceCentroid));
}


double DevelopableMesh::averageWarpedMaterialMesh(Eigen::VectorXd &q, WarpedMesh &warped)
{
    int nembverts = mesh_.n_vertices();
    int nmatverts = material_->getMesh().n_vertices();
    map<int, Vector3d> newpos;
    map<int, int> numsamples;

    for(OMMesh::VertexIter vi = warped.mesh.getMesh().vertices_begin(); vi != warped.mesh.getMesh().vertices_end(); ++vi)
    {
        int vid = vi.handle().idx();
        int matid = warped.vertexMap[vid];
        OMMesh::Point pt = warped.mesh.getMesh().point(vi.handle());
        Vector3d ept;
        for(int j=0; j<3; j++)
            ept[j] = pt[j];
        map<int, Vector3d>::iterator mit = newpos.find(matid);
        if(mit == newpos.end())
            newpos[matid] = ept;
        else
            newpos[matid] += ept;
        numsamples[matid]++;
    }

    double error=0;
    for(int i=0; i<nmatverts; i++)
    {
        assert(numsamples[i] > 0);
        Vector3d newpt = newpos[i]/numsamples[i];

        for(map<int,int>::iterator it = warped.vertexMap.begin(); it != warped.vertexMap.end(); ++it)
        {
            if(i == it->second)
            {
                OMMesh::VertexHandle v = warped.mesh.getMesh().vertex_handle(it->first);
                OMMesh::Point oldpt = warped.mesh.getMesh().point(v);
                Vector3d eoldpt;
                for(int j=0; j<3; j++)
                    eoldpt[j] = oldpt[j];
                error += (eoldpt.segment<2>(0)-newpt.segment<2>(0)).squaredNorm();
            }
        }
        q.segment(3*nembverts+2*i, 2) = newpt.segment<2>(0);
    }

    return error;
}

double DevelopableMesh::averageWarpedEmbeddedMesh(Eigen::VectorXd &q, WarpedMesh &warped)
{
    int nembverts = mesh_.n_vertices();
    map<int, Vector3d> newpos;
    map<int, int> numsamples;

    for(OMMesh::VertexIter vi = warped.mesh.getMesh().vertices_begin(); vi != warped.mesh.getMesh().vertices_end(); ++vi)
    {
        int vid = vi.handle().idx();
        int embid = warped.vertexMap[vid];
        OMMesh::Point pt = warped.mesh.getMesh().point(vi.handle());
        Vector3d ept;
        for(int j=0; j<3; j++)
            ept[j] = pt[j];
        map<int, Vector3d>::iterator mit = newpos.find(embid);
        if(mit == newpos.end())
            newpos[embid] = ept;
        else
            newpos[embid] += ept;
        numsamples[embid]++;
    }

    double error=0;
    for(int i=0; i<nembverts; i++)
    {
        assert(numsamples[i] > 0);
        Vector3d newpt = newpos[i]/numsamples[i];
        for(map<int,int>::iterator it = warped.vertexMap.begin(); it != warped.vertexMap.end(); ++it)
        {
            if(i == it->second)
            {
                OMMesh::VertexHandle v = warped.mesh.getMesh().vertex_handle(it->first);
                OMMesh::Point oldpt = warped.mesh.getMesh().point(v);
                Vector3d eoldpt;
                for(int j=0; j<3; j++)
                    eoldpt[j] = oldpt[j];
                error += (eoldpt-newpt).squaredNorm();
            }
        }
        q.segment<3>(3*i) = newpt;
    }

    return error;
}

void DevelopableMesh::warpEmbeddedToMaterial(const VectorXd &q, WarpedMesh &result)
{
    result.mesh.getMesh() = OMMesh();
    result.vertexMap.clear();

    int nembverts = mesh_.n_vertices();
    int nmatverts = material_->getMesh().n_vertices();
    assert(nembverts == nmatverts);
    assert(q.size() == 3*nembverts + 2*nmatverts);
    VectorXd embq = q.segment(0,3*nembverts);
    VectorXd matq = q.segment(3*nembverts, 2*nmatverts);

    int newvertidx=0;

    for(OMMesh::FaceIter fi = material_->getMesh().faces_begin(); fi != material_->getMesh().faces_end(); ++fi)
    {
        OMMesh::FaceHandle fh = fi.handle();
        Vector3d totoffset(0,0,0);
        vector<Vector3d> matverts;
        vector<Vector3d> embverts;
        vector<int> matvertsid;

        for(OMMesh::FaceHalfedgeIter fhi = material_->getMesh().fh_iter(fh); fhi; ++fhi)
        {
            int vertid = material_->getMesh().from_vertex_handle(fhi).idx();
            embverts.push_back(embq.segment<3>(3*vertid));
            matvertsid.push_back(vertid);
            Vector3d matvert3d;
            matvert3d.segment<2>(0) = matq.segment<2>(2*vertid);
            matvert3d[2] = 0;
            matverts.push_back(matvert3d+totoffset);

            totoffset += material_->getOffset(fhi.handle());
        }

        vector<Vector3d> newtri;
        shapeMatch(embverts, matverts, newtri);

        totoffset.setZero();
        vector<OMMesh::VertexHandle> newvhs;
        int i=0;
        for(OMMesh::FaceHalfedgeIter fhi = material_->getMesh().fh_iter(fh); fhi; ++fhi)
        {
            OMMesh::Point newpt;
            for(int j=0; j<3; j++)
                newpt[j] = newtri[i][j] - totoffset[j];
            newvhs.push_back(result.mesh.getMesh().add_vertex(newpt));
            result.vertexMap[newvertidx++] = matvertsid[i];

            totoffset += material_->getOffset(fhi.handle());
            i++;
        }
        result.mesh.getMesh().add_face(newvhs);
    }
}

void DevelopableMesh::warpMaterialToEmbedded(const Eigen::VectorXd &q, WarpedMesh &result)
{
    result.mesh.getMesh() = OMMesh();
    result.vertexMap.clear();

    int nembverts = mesh_.n_vertices();
    int nmatverts = material_->getMesh().n_vertices();
    assert(nembverts == nmatverts);
    assert(q.size() == 3*nembverts + 2*nmatverts);
    VectorXd embq = q.segment(0,3*nembverts);
    VectorXd matq = q.segment(3*nembverts, 2*nmatverts);

    int newvertidx=0;

    for(OMMesh::FaceIter fi = mesh_.faces_begin(); fi != mesh_.faces_end(); ++fi)
    {
        OMMesh::FaceHandle fh = fi.handle();
        Vector3d totoffset(0,0,0);
        vector<Vector3d> matverts;
        vector<Vector3d> embverts;
        vector<int> embvertsid;

        for(OMMesh::FaceHalfedgeIter fhi = mesh_.fh_iter(fh); fhi; ++fhi)
        {
            int vertid = mesh_.from_vertex_handle(fhi).idx();
            embverts.push_back(embq.segment<3>(3*vertid));
            embvertsid.push_back(vertid);
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
        int numtriverts = newtri.size();
        assert(numtriverts == (int)embvertsid.size());
        for(int i=0; i<numtriverts; i++)
        {
            OMMesh::Point newpt;
            for(int j=0; j<3; j++)
                newpt[j] = newtri[i][j];
            newvhs.push_back(result.mesh.getMesh().add_vertex(newpt));
            result.vertexMap[newvertidx++] = embvertsid[i];
        }
        result.mesh.getMesh().add_face(newvhs);
    }
}

void DevelopableMesh::enforceBoundaryConstraints(VectorXd &q)
{
    for(int i=0; i<(int)boundaries_.size(); i++)
    {
        for(int j=0; j<(int)boundaries_[i].bdryVerts.size(); j++)
        {
            Vector3d targetpt = boundaries_[i].bdryPos[j];
            int vidx = boundaries_[i].bdryVerts[j];
            q.segment<3>(3*vidx) = targetpt;
        }
    }
    // TODO material constraints
}
