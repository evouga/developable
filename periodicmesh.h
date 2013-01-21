#ifndef PERIODICMESH_H
#define PERIODICMESH_H

#include "mesh.h"
#include <map>

class PeriodicMesh : public Mesh
{
public:
    PeriodicMesh();

    void clearOffsets();
    void addOffset(OMMesh::HalfedgeHandle hei, const Eigen::Vector3d &offset);
    Eigen::Vector3d getOffset(OMMesh::HalfedgeHandle hei);

    virtual double edgeLength(int vid1, int vid2);
    using Mesh::edgeLength;

private:
    std::map<OMMesh::HalfedgeHandle, Eigen::Vector3d> offsets_;
};

#endif // PERIODICMESH_H
