#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <Eigen/Core>
#include <vector>

class Mesh;

struct CurvatureInfo
{
    CurvatureInfo();
    double principalCurvature[2];
    Eigen::Vector3d curvatureDir[2];
};

class MeshCurvature
{
public:
    MeshCurvature(Mesh &mesh);
    void renderCurvatureDirs(Mesh &mesh);
    double gaussianCurvature(int vidx);
    double meanCurvature(int vidx);

private:
    std::vector<CurvatureInfo> curvature_;
    void computeFrame(const Eigen::Vector3d &normal, Eigen::Vector3d &u, Eigen::Vector3d &v);

    void computeCurvatures(Mesh &mesh);
    CurvatureInfo computeCurvature(Mesh &mesh, int vidx);
    void estimateCurvature(Mesh &mesh, int vidx, Eigen::Vector3d &normal, CurvatureInfo &curvature);

};

#endif // GEOMETRY_H
