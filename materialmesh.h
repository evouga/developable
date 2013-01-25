#ifndef MATERIALMESH_H
#define MATERIALMESH_H

#include "periodicmesh.h"
#include <vector>
#include <map>

struct MaterialBoundary
{
    int vertid;
    Eigen::Matrix<double,2,1,Eigen::DontAlign> pos;
};

class MaterialMesh : public PeriodicMesh
{
public:
    MaterialMesh(double H) : H_(H) {}

    virtual bool saveToStream(std::ostream &os);
    virtual bool loadFromStream(std::istream &is);

    double getH() {return H_;}

    std::vector<MaterialBoundary> &getBoundaryVerts() {return bdryverts_;}
    int materialEdge(int embeddedEdge);

    void setMaterialEdge(int embeddedEdge, int materialEdge);

private:
    double H_;
    std::vector<MaterialBoundary> bdryverts_;
    std::map<int, int> embedge2matedge_;
};


#endif // MATERIALMESH_H
