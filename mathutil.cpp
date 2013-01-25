#include "mathutil.h"

const double MathUtil::PI = 3.1415926535898;

using namespace Eigen;
using namespace std;

typedef Triplet<double> T;

double MathUtil::randomDouble(double min, double max)
{
    double range = max-min;
    return range * (double(rand())/double(RAND_MAX));
}

bool MathUtil::hasNANs(const Eigen::VectorXd &v)
{
    for(int i=0; i<(int)v.size(); i++)
        if(isnan(v[i]))
            return true;
    return false;
}

bool MathUtil::hasNANs(const std::vector<T> &M)
{
    for(int i=0; i<(int)M.size(); i++)
        if(isnan(M[i].value()))
            return true;
    return false;
}

bool MathUtil::hasNANs(const std::vector<std::vector<T> > &H)
{
    for(int i=0; i<(int)H.size(); i++)
        for(int j=0; j<(int)H[i].size(); j++)
        {
            if(isnan(H[i][j].value()))
                return true;
        }
    return false;
}

void MathUtil::symmetricMatrixFromTriplets(const std::vector<T> &triplets, Eigen::SparseMatrix<double> &M)
{
    vector<T> augtriplets = triplets;
    for(vector<T>::const_iterator it = triplets.begin(); it != triplets.end(); ++it)
    {
        if(it->row() < it->col())
            augtriplets.push_back(T(it->col(), it->row(), it->value()));
    }
    M.setFromTriplets(augtriplets.begin(), augtriplets.end());
}
