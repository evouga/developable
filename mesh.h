#ifndef MESH_H
#define MESH_H

#include <string>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <Eigen/Core>

struct MyTraits : public OpenMesh::DefaultTraits
{
    typedef OpenMesh::Vec3d Point; // use double-values points
    typedef OpenMesh::Vec3d Normal; // use double-values points

    VertexTraits
    {
    private:
        Point vel_;
    public:
        VertexT() : vel_(Point(0.0f, 0.0f, 0.0f) ) {}
        const Point &vel() const {return vel_;}
        void set_vel(const Point &v) {vel_=v;}
    };
};

typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits>  OMMesh;

class GLUquadric;

class Mesh
{
public:
    Mesh();
    virtual ~Mesh();

    bool exportOBJ(const char *filename);
    bool importOBJ(const char *filename);

    virtual bool loadFromStream(std::istream &is);
    virtual bool saveToStream(std::ostream &os);
    OMMesh &getMesh() {return mesh_;}
    const OMMesh &getMesh() const {return mesh_;}
    void render(bool showWireframe, bool smoothShade);

    Eigen::Vector3d centroid();
    double radius();

    Eigen::Vector3d vertexNormal(int vidx) const;
    double shortestAdjacentEdge(int vidx) const;

    double areaOfInfluence(int vidx) const;
    double faceArea(int fidx) const;

    int findEdge(int vid1, int vid2);
    int findHalfedge(int vid1, int vid2);
    virtual double edgeLength(int vid1, int vid2);
    double edgeLength(int eid);


protected:
    OMMesh mesh_;
    void edgeEndpoints(OMMesh::EdgeHandle eh, OMMesh::Point &pt1, OMMesh::Point &pt2);

    void writeInt(std::ostream &os, int i);
    void writeDouble(std::ostream &os, double d);
    void writeBool(std::ostream &os, bool b);
    int readInt(std::istream &is);
    double readDouble(std::istream &is);
    bool readBool(std::istream &is);

private:
    Mesh(const Mesh &other);
    Mesh &operator=(const Mesh &other);
    GLUquadric *quadric_;

    void drawSphere(int vertex);
};

#endif // MESH_H
