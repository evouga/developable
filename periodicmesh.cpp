#include "periodicmesh.h"

using namespace std;
using namespace Eigen;

PeriodicMesh::PeriodicMesh()
{

}

void PeriodicMesh::clearOffsets()
{
    offsets_.clear();
}

void PeriodicMesh::addOffset(OMMesh::HalfedgeHandle hei, const Eigen::Vector3d &offset)
{
    offsets_[hei] = offset;
    offsets_[mesh_.opposite_halfedge_handle(hei)] = -offset;
}

Eigen::Vector3d PeriodicMesh::getOffset(OMMesh::HalfedgeHandle hei)
{
    Vector3d offset;
    map<OMMesh::HalfedgeHandle, Eigen::Vector3d>::iterator it = offsets_.find(hei);
    if(it != offsets_.end())
        offset = it->second;
    else
        offset.setZero();

    return offset;
}

double PeriodicMesh::edgeLength(int vid1, int vid2)
{
    OMMesh::VertexHandle vh1 = mesh_.vertex_handle(vid1);
    OMMesh::VertexHandle vh2 = mesh_.vertex_handle(vid2);
    Vector3d diff;
    for(int k=0; k<3; k++)
        diff[k] = mesh_.point(vh2)[k] - mesh_.point(vh1)[k];

    int hehid = findHalfedge(vid1, vid2);
    Vector3d offset = getOffset(mesh_.halfedge_handle(hehid));
    return (diff+offset).norm();
}
