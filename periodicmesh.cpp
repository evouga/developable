#include "periodicmesh.h"

using namespace std;
using namespace Eigen;

PeriodicMesh::PeriodicMesh()
{

}

bool PeriodicMesh::saveToStream(std::ostream &os)
{
    if(!Mesh::saveToStream(os))
        return false;

    int noffsets = offsets_.size();
    writeInt(os, noffsets);
    for(map<OMMesh::HalfedgeHandle, Vector3d>::iterator it = offsets_.begin(); it != offsets_.end(); ++it)
    {
        int sid = mesh_.from_vertex_handle(it->first).idx();
        int did = mesh_.to_vertex_handle(it->first).idx();
        writeInt(os, sid);
        writeInt(os, did);
        Vector3d offset = it->second;
        writeDouble(os, offset[0]);
        writeDouble(os, offset[1]);
        writeDouble(os, offset[2]);
    }
    return os;
}

bool PeriodicMesh::loadFromStream(std::istream &is)
{
    if(!Mesh::loadFromStream(is))
        return false;

    int noffsets = readInt(is);
    if(!is)
        return false;
    offsets_.clear();
    for(int i=0; i<noffsets; i++)
    {
        int sid = readInt(is);
        int did = readInt(is);
        if(!is)
            return false;
        int hid = findHalfedge(sid, did);
        assert(hid >= 0);
        OMMesh::HalfedgeHandle heh = mesh_.halfedge_handle(hid);
        Vector3d offset;
        offset[0] = readDouble(is);
        offset[1] = readDouble(is);
        offset[2] = readDouble(is);
        if(!is)
            return false;
        offsets_[heh] = offset;
    }
    return true;
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
