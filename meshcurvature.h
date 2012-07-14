#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <Eigen/Core>
#include <vector>

class Mesh;

struct Frame
{
    Frame(const Eigen::Vector3d &normal);
    Eigen::Vector3d normal;
    Eigen::Vector3d u;
    Eigen::Vector3d v;
};

struct ShapeOperator
{
    ShapeOperator(Frame frame, const Eigen::Matrix2d &S) : frame(frame), S(S) {}
    ShapeOperator(Frame frame) : frame(frame) {S.setZero();}
    Frame frame;
    Eigen::Matrix2d S;
};

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
    double curvatureSpread(int vidx);

    void totalSquaredGaussianCurvature(const Mesh &mesh, std::pair<double, double> &values, double cutoff);

private:
    std::vector<CurvatureInfo> curvature_;

    double curCutoff_;
    double aboveSqGaussCurvature_;
    double belowSqGaussCurvature_;

    void computeCurvatures(Mesh &mesh);
    void computeCurvature(const ShapeOperator &shapeOperator, CurvatureInfo &curvature);
    void computeShapeOperators(Mesh &mesh, std::vector<ShapeOperator> &operators);
    void recomputeSquaredGaussianCurvature(const Mesh &mesh, double cutoff);
    ShapeOperator computeShapeOperator(Mesh &mesh, int vidx);
    ShapeOperator estimateShapeOperator(Mesh &mesh, int vidx, Eigen::Vector3d &normal);

    Eigen::Matrix2d transportShapeOperator(const ShapeOperator &source, const ShapeOperator &dest);
    void smoothShapeOperators(const Mesh &m, const std::vector<ShapeOperator> &oldOperators, std::vector<ShapeOperator> &newOperators);
};

#endif // GEOMETRY_H
