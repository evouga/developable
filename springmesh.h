#ifndef SPRINGMESH_H
#define SPRINGMESH_H

#include "developablemesh.h"
#include <vector>
#include <string>
#include <map>
#include "mesh.h"
#include <Eigen/Sparse>
#include <FADBAD++/fadiff.h>
#include "autodiffTemplates.h"

class SpringMesh: public DevelopableMesh
{
public:
    double elasticEnergy(const Eigen::VectorXd &q);
    void setupSprings();

 private:
    std::vector<std::vector<int> > springs;
    std::vector<int> baseverts;
    double spring_constant;
};

#endif // SPRINGMESH_H
