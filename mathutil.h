#ifndef MATHUTIL_H
#define MATHUTIL_H

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>

class MathUtil
{
public:
    const static double PI;

    static bool hasNANs(const Eigen::VectorXd &v);
    static bool hasNANs(const std::vector<Eigen::Triplet<double> > &M);
    static bool hasNANs(const std::vector<std::vector<Eigen::Triplet<double> > > &H);

    static double randomDouble(double min, double max);
    static void symmetricMatrixFromTriplets(const std::vector<Eigen::Triplet<double> > &triplets, Eigen::SparseMatrix<double> &M);

private:
    MathUtil();
};

#endif // MATHUTIL_H
