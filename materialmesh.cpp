#include "materialmesh.h"
#include <map>

using namespace std;
using namespace Eigen;

Vector2d MaterialBoundary::getPos()
{
    return Eigen::Vector2d(xpos, onBottom ? 0.0 : m_->getH());
}

bool MaterialMesh::saveToStream(std::ostream &os)
{
    if(!PeriodicMesh::saveToStream(os))
        return false;

    writeDouble(os, H_);

    int nbdryverts = bdryverts_.size();
    writeInt(os, nbdryverts);
    for(int i=0; i<nbdryverts; i++)
    {
        writeInt(os, bdryverts_[i].vertid);
        writeDouble(os, bdryverts_[i].xpos);
        writeBool(os, bdryverts_[i].onBottom);
    }

    int nmapentries = embedge2matedge_.size();
    writeInt(os, nmapentries);
    for(map<int, int>::iterator it = embedge2matedge_.begin(); it != embedge2matedge_.end(); ++it)
    {
        writeInt(os, it->first);
        writeInt(os, it->second);
    }
    return os;
}

bool MaterialMesh::loadFromStream(std::istream &is)
{
    if(!PeriodicMesh::loadFromStream(is))
        return false;

    H_ = readDouble(is);
    if(!is)
        return false;

    int nbdryverts = readInt(is);
    if(!is)
        return false;

    bdryverts_.clear();
    for(int i=0; i<nbdryverts; i++)
    {
        MaterialBoundary bdinfo(this);
        bdinfo.vertid = readInt(is);
        bdinfo.xpos = readDouble(is);
        bdinfo.onBottom = readBool(is);
        if(!is)
            return false;
        bdryverts_.push_back(bdinfo);
    }

    int nmapentries = readInt(is);
    if(!is)
    {
        assert(false);
        return false;
    }
    embedge2matedge_.clear();
    for(int i=0; i<nmapentries; i++)
    {
        int idx1 = readInt(is);
        int idx2 = readInt(is);
        if(!is)
        {
            assert(false);
            return false;
        }
        embedge2matedge_[idx1] = idx2;
    }
    return true;
}

void MaterialMesh::setMaterialEdge(int embeddedEdge, int materialEdge)
{
    assert(embeddedEdge != -1);
    assert(materialEdge != -1);
    embedge2matedge_[embeddedEdge] = materialEdge;
}

int MaterialMesh::materialEdge(int embeddedEdge)
{
    map<int, int>::iterator it = embedge2matedge_.find(embeddedEdge);
    if(it != embedge2matedge_.end())
        return it->second;
    return -1;
}
