#include "developablemesh.h"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Geometry>
#include <set>
#include "coin/IpIpoptApplication.hpp"
#include <QGLWidget>
#include <Eigen/Dense>
#include "mathutil.h"
#include "projectionnlp.h"

using namespace Eigen;
using namespace std;
using namespace Ipopt;

DevelopableMesh::DevelopableMesh() : material_(new MaterialMesh(0))
{

}

DevelopableMesh::~DevelopableMesh()
{
    delete material_;
}

bool DevelopableMesh::saveToStream(std::ostream &os)
{
    if(!Mesh::saveToStream(os))
    {
        return false;
    }
    if(!material_->saveToStream(os))
        return false;

    int nboundaries = boundaries_.size();
    writeInt(os, nboundaries);
    for(int i=0; i<nboundaries; i++)
    {
        int nverts = boundaries_[i].bdryPos.size();
        assert(nverts == (int)boundaries_[i].bdryVerts.size());
        writeInt(os, nverts);
        for(int j=0; j<nverts; j++)
        {
            writeInt(os, boundaries_[i].bdryVerts[j]);
            Vector3d pos = boundaries_[i].bdryPos[j];
            writeDouble(os, pos[0]);
            writeDouble(os, pos[1]);
            writeDouble(os, pos[2]);
        }
    }
    return os;
}

bool DevelopableMesh::loadFromStream(std::istream &is)
{
    if(!Mesh::loadFromStream(is))
    {
        return false;
    }
    if(!material_->loadFromStream(is))
    {
        assert(false);
        return false;
    }

    int nboundaries = readInt(is);
    if(!is)
    {
        assert(false);
        return false;
    }

    boundaries_.clear();

    for(int i=0; i<nboundaries; i++)
    {
        int nverts = readInt(is);
        if(!is)
        {
            assert(false);
            return false;
        }
        Boundary newbd;
        for(int j=0; j<nverts; j++)
        {
            int vertid = readInt(is);
            Vector3d pos;
            pos[0] = readDouble(is);
            pos[1] = readDouble(is);
            pos[2] = readDouble(is);
            cout << pos.transpose() << endl;
            if(!is)
            {
                 return false;
            }
            newbd.bdryVerts.push_back(vertid);
            newbd.bdryPos.push_back(pos);
        }
        boundaries_.push_back(newbd);
    }
    return true;
}

Vector3d DevelopableMesh::point2Vector(OMMesh::Point pt)
{
    Vector3d result;
    for(int i=0; i<3; i++)
        result[i] = pt[i];
    return result;
}

void DevelopableMesh::centerCylinder()
{
    double totx=0;
    double toty=0;
    for(int i=0; i<(int)mesh_.n_vertices(); i++)
    {
        OMMesh::Point pt = mesh_.point(mesh_.vertex_handle(i));
        totx += pt[0];
        toty += pt[1];
    }
    totx /= mesh_.n_vertices();
    toty /= mesh_.n_vertices();

    for(int i=0; i<(int)mesh_.n_vertices(); i++)
    {
        mesh_.point(mesh_.vertex_handle(i))[0] -= totx;
        mesh_.point(mesh_.vertex_handle(i))[1] -= toty;
    }

    double totxm = 0;
    for(int i=0; i<(int)material_->getMesh().n_vertices(); i++)
    {
        OMMesh::Point pt = material_->getMesh().point(material_->getMesh().vertex_handle(i));
        totxm += pt[0];
    }
    totxm /= material_->getMesh().n_vertices();
    for(int i=0; i<(int)material_->getMesh().n_vertices(); i++)
    {
        material_->getMesh().point(material_->getMesh().vertex_handle(i))[0] -= totxm;
    }
}

void DevelopableMesh::projectOntoConstraintManifold(DeformCallback &dc)
{
    VectorXd q;
    VectorXd v;
    gatherDOFs(q,v);

    cout << "Initial equality violation: " << equalityConstraintViolation(q) << endl;

    SmartPtr<ProjectionNLP> pnlp = new ProjectionNLP(*this, dc, q);

    IpoptApplication app;
//    app.Options()->SetNumericValue("tol", 1e-9);
//    app.Options()->SetNumericValue(tol, 1e-6);
//    app.Options()->SetStringValue("derivative_test", "second_order");
//    app.Options()->SetStringValue("check_derivatives_for_naninf", "yes");

    app.Initialize();

    app.OptimizeTNLP(pnlp);

    q = pnlp->getFinalQ();

    flushOutNANs(q);

    cout << "Final equality violation: " << equalityConstraintViolation(q) << endl;

    repopulateDOFs(q,v);
}

void DevelopableMesh::crushLantern(DeformCallback &dc, double dt)
{
    double crushspeed = 0.1;
    int steps = (int)(1.0/crushspeed/dt);
//    int steps = 100;
    int framestep = steps/1000;
//    int framestep = steps;
    VectorXd q, v;
    gatherDOFs(q,v);

    for(int j=0; j<(int)boundaries_[0].bdryPos.size(); j++)
    {
            cout << j << "\t" << boundaries_[0].bdryVerts[j] << " " << boundaries_[0].bdryPos[j][2] << endl;
    }
    for(int j=0; j<(int)boundaries_[1].bdryPos.size(); j++)
    {
        cout << j << "\t" << boundaries_[1].bdryVerts[j] << " " << boundaries_[1].bdryPos[j][2] << endl;
    }
    cout << boundaries_[0].bdryVerts.size() << endl;
    cout << boundaries_[1].bdryVerts.size() << endl;
    cout << "Closing distance error: " << q.segment<3>(3*boundaries_[0].bdryVerts[4]) << " vs " << q.segment<3>(boundaries_[0].bdryVerts[0]) << endl;
    cout << "Constraint violation: " << equalityConstraintViolation(q) << endl;

    cout << "crushing!" << endl;
    for(int i=0; i<steps; i++)
    {
//        v *= 0.9;
//        double k = 100;
//        q += dt*v;
//        repopulateDOFs(q, v);
        for(int j=0; j<(int)boundaries_[1].bdryPos.size(); j++)
        {
                boundaries_[1].bdryPos[j][2] -= crushspeed*dt;
        }
        cout << "Projecting!" << endl;
        projectOntoConstraintManifold(dc);

//        int collapseid = findCollapsibleEdge(q);
//        while(collapseid != -1)
//        {
//            cout << "Collapse!" << endl;
//            assert(canCollapseEdge(collapseid));
//            collapseEdge(collapseid);
//            projectOntoConstraintManifold(dc);
//            collapseid = findCollapsibleEdge(q);
//        }
//        gatherDOFs(q, v);

//        double f;
//        VectorXd Df;
//        vector<T> Hf;
//        buildObjective(q, f, Df, Hf);
//        v += -dt*k*Df;
        if(i%framestep == 0)
            dc.repaintCallback();
    }
//    repopulateDOFs(q,v);
    gatherDOFs(q,v);

    for(int i = 0; i < (int)boundaries_.size(); i++)
    {
        for(int j=0; j<(int)boundaries_[i].bdryVerts.size(); j++)
        {
            cout << j << "th vector should\n" << boundaries_[i].bdryPos[j] << endl;
            cout << j << "th vector is\n" << q.segment<3>(3*boundaries_[i].bdryVerts[j]) << endl;
        }
    }
    cout << "Crushed!" << endl;
    cout << "Closing points error: " << q.segment<3>(3*4) << " vs " << q.segment<3>(0) << endl;
    cout << "Closing distance error: " << (q.segment<3>(3*4)-q.segment<3>(0)).squaredNorm()<< endl;
    cout << "Constraint violation: " << equalityConstraintViolation(q) << endl;
}



bool DevelopableMesh::shouldCollapseEdges(const Eigen::VectorXd &q)
{
    return (findCollapsibleEdge(q) != -1);
}

Vector2d DevelopableMesh::materialCenter()
{
    Vector2d result(0,0);
    for(int i=0; i<(int)material_->getMesh().n_vertices(); i++)
    {
        OMMesh::VertexHandle vh = material_->getMesh().vertex_handle(i);
        OMMesh::Point pt = material_->getMesh().point(vh);
        result[0] += pt[0];
        result[1] += pt[1];
    }
    result /= material_->getMesh().n_vertices();
    return result;
}

double DevelopableMesh::materialRadius()
{
    Vector2d center = materialCenter();
    double maxradius = 0;
    for(int i=0; i<(int)material_->getMesh().n_vertices(); i++)
    {
        OMMesh::VertexHandle vh = material_->getMesh().vertex_handle(i);
        OMMesh::Point pt = material_->getMesh().point(vh);
        Vector2d vpt(pt[0],pt[1]);
        maxradius = max(maxradius, (center-vpt).norm());
        for(OMMesh::VertexOHalfedgeIter voh = material_->getMesh().voh_iter(vh); voh.is_valid(); ++voh)
        {
            OMMesh::HalfedgeHandle heh = voh.handle();
            OMMesh::Point adjompt = material_->getMesh().point(material_->getMesh().to_vertex_handle(heh));
            Vector2d adjpt(adjompt[0],adjompt[1]);
            adjpt += material_->getOffset(heh).segment<2>(0);
            maxradius = max(maxradius, (center-adjpt).norm());
        }
    }
    return maxradius;
}

void DevelopableMesh::renderMaterial()
{
//    glColor3f(0.0,0.0,0.0);
//    glBegin(GL_LINES);
//    for(int i=0; i<(int)material_->getMesh().n_halfedges(); i++)
//    {
//        OMMesh::HalfedgeHandle heh = material_->getMesh().halfedge_handle(i);
//        OMMesh::Point pt1 = material_->getMesh().point(material_->getMesh().from_vertex_handle(heh));
//        OMMesh::Point pt2 = material_->getMesh().point(material_->getMesh().to_vertex_handle(heh));
//        Vector3d offset = material_->getOffset(heh);
//        glVertex2d(pt1[0], pt1[1]);
//        glVertex2d(pt2[0] + offset[0], pt2[1] + offset[1]);
//    }
//    glEnd();

//    glPointSize(6.0);
//    glBegin(GL_POINTS);
//    for(int i=0; i<(int)material_->getMesh().n_vertices(); i++)
//    {
//        glColor3f(0.2, 0.0, 1.0);
//        OMMesh::VertexHandle vh = material_->getMesh().vertex_handle(i);
//        OMMesh::Point pt = material_->getMesh().point(vh);
//        glVertex2d(pt[0], pt[1]);
//    }
//    glColor3f(1.0, 0.0, 0.0);
//    for(int i=0; i<(int)material_->getBoundaryVerts().size(); i++)
//    {
//        OMMesh::VertexHandle vh = material_->getMesh().vertex_handle(material_->getBoundaryVerts()[i].vertid);
//        OMMesh::Point pt = material_->getMesh().point(vh);
//        glVertex2d(pt[0], pt[1]);
//    }
//    glEnd();
//    glPointSize(1.0);
}

void DevelopableMesh::gatherDOFs(Eigen::VectorXd &q, Eigen::VectorXd &v)
{
    int nembverts = mesh_.n_vertices();
//    int nmatverts = material_->getMesh().n_vertices();

//    q.resize(3*nembverts + 2*nmatverts);
//    v.resize(3*nembverts + 2*nmatverts);
    q.resize(3*nembverts);
    v.resize(3*nembverts);

    for(int i=0; i<nembverts; i++)
    {
        OMMesh::Point vel = mesh_.data(mesh_.vertex_handle(i)).vel();
        for(int j=0; j<3; j++)
        {
            q[3*i+j] = mesh_.point(mesh_.vertex_handle(i))[j];
            v[3*i+j] = vel[j];
        }
    }

//    for(int i=0; i<nmatverts; i++)
//    {
//        OMMesh::Point vel = material_->getMesh().data(material_->getMesh().vertex_handle(i)).vel();
//        for(int j=0; j<2; j++)
//        {
//            q[3*nembverts+2*i+j] = material_->getMesh().point(material_->getMesh().vertex_handle(i))[j];
//            v[3*nembverts+2*i+j] = vel[j];
//        }
//    }
}

void DevelopableMesh::repopulateDOFs(const Eigen::VectorXd &q, const Eigen::VectorXd &v)
{
    int nembverts = mesh_.n_vertices();
//    int nmatverts = material_->getMesh().n_vertices();

    for(int i=0; i<nembverts; i++)
    {
        OMMesh::Point newvel;
        for(int j=0; j<3; j++)
        {
            mesh_.point(mesh_.vertex_handle(i))[j] = q[3*i+j];
            newvel[j] = v[3*i+j];
        }
        mesh_.data(mesh_.vertex_handle(i)).set_vel(newvel);
    }

//    for(int i=0; i<nmatverts; i++)
//    {
//        OMMesh::Point newvel;
//        for(int j=0; j<2; j++)
//        {
//            material_->getMesh().point(material_->getMesh().vertex_handle(i))[j] = q[3*nembverts+2*i+j];
//            newvel[j] = v[3*nembverts+2*i+j];
//        }
//        newvel[2] = 0;
//        OMMesh &matmesh = material_->getMesh();
//        matmesh.data(material_->getMesh().vertex_handle(i)).set_vel(newvel);
//    }
}

int DevelopableMesh::activeInequalityConstraints(const Eigen::VectorXd &q)
{
    vector<T> Dh;
    vector<vector<T> > Hh;
    VectorXd h;
    buildInversionConstraints(q, h, Dh, Hh);
    int result=0;
    for(int i=0; i<h.size(); i++)
    {
        if(h[i] < 1e-6)
            result++;
    }
    return result;
}

void DevelopableMesh::buildConstraintBasis(const Eigen::VectorXd &q, std::vector<Eigen::VectorXd> &normalspace, std::vector<Eigen::VectorXd> &tangentspace)
{
    VectorXd g;
    vector<T> Dg;
    vector<vector<T> > Hg;
    buildConstraints(q, g, Dg, Hg);

    SparseMatrix<double> sparse(g.size(), q.size());
    sparse.setFromTriplets(Dg.begin(), Dg.end());
    MatrixXd dense(sparse);
    JacobiSVD<MatrixXd> svd(dense, ComputeFullU | ComputeFullV);
    VectorXd singvals(q.size());
    singvals.setZero();
    for(int i=0; i < svd.singularValues().size() && i < q.size(); i++)
        singvals[i] = svd.singularValues()[i];

    normalspace.clear();
    tangentspace.clear();
    for(int i=0; i<q.size(); i++)
    {
        VectorXd col = svd.matrixV().col(i);
        if(singvals[i] > 1e-6)
            normalspace.push_back(col);
        else
            tangentspace.push_back(col);
    }
}

DevelopableMesh::CriticalPointType DevelopableMesh::checkConstrainedHessian(const Eigen::VectorXd &q, vector<VectorXd> &negdirs)
{
    vector<VectorXd> tspace;
    vector<VectorXd> nspace;
    buildConstraintBasis(q, nspace, tspace);

    MatrixXd Tan(tspace.size(), q.size());
    Tan.setZero();
    for(int i=0; i<(int)tspace.size(); i++)
        Tan.row(i) = tspace[i];

    double f;
    VectorXd Df;
    vector<T> Hf;
    buildObjective(q, f, Df, Hf);
    SparseMatrix<double> sparse(q.size(), q.size());
    MathUtil::symmetricMatrixFromTriplets(Hf, sparse);
    MatrixXd constrainedH = Tan*sparse*Tan.transpose();

    SelfAdjointEigenSolver<MatrixXd> es(constrainedH);
    VectorXd evs = es.eigenvalues();
    double mostneg = -1e-6;
    double mostpos = 1e-6;
    double startneg = mostneg;
    double startpos = mostpos;

    negdirs.clear();

    for(int i=0; i<(int)evs.size(); i++)
    {
        if(evs[i] < startneg)
        {
            negdirs.push_back(Tan.transpose()*es.eigenvectors().col(i));
        }
        if(evs[i] < mostneg)
            mostneg = evs[i];
        if(evs[i] > mostpos)
            mostpos = evs[i];
    }

    if(mostneg < startneg && mostpos > startpos)
    {
        std::cout << "Saddle Point, " << mostneg << " " << mostpos << std::endl;
        return CPT_SADDLE;
    }
    else if(mostneg < startneg)
    {
        std::cout << "Maximum" << std::endl;
        return CPT_MAXIMUM;
    }
    else if(mostpos > startpos)
    {
        std::cout << "Minimum" << std::endl;
        return CPT_MINIMUM;
    }
    else
    {
        std::cout << "Singular" << std::endl;
        return CPT_INDETERMINATE;
    }
}

void DevelopableMesh::flushOutNANs(const Eigen::VectorXd &q)
{
    assert(!MathUtil::hasNANs(q));
    double f;
    VectorXd Df;
    vector<T> Hf;
    buildObjective(q, f, Df, Hf);
    assert(!isnan(f));
    assert(!MathUtil::hasNANs(Df));
    assert(!MathUtil::hasNANs(Hf));

    VectorXd g;
    vector<T> Dg;
    vector<vector<T> > Hg;
    buildConstraints(q, g, Dg, Hg);

    assert(!MathUtil::hasNANs(g));
    assert(!MathUtil::hasNANs(Dg));
    assert(!MathUtil::hasNANs(Hg));
}

double DevelopableMesh::equalityConstraintViolation(const Eigen::VectorXd &q)
{
    VectorXd g;
    vector<T> Dg;
    vector<vector<T> > Hg;
    buildConstraints(q, g, Dg, Hg);
    double max = 0.0;
    for(int i=0; i<(int)g.size(); i++)
        if(fabs(g[i]) > max)
            max = fabs(g[i]);
    return max;
}

void DevelopableMesh::perturbConfiguration(Eigen::VectorXd &q, const std::vector<VectorXd> &dirs, double mag)
{
    cout << "Perturbing along " << dirs.size() << " negative directions" << endl;
    gatherInfo(q);

    double oldf;
    VectorXd Df;
    vector<T> Hf;
    buildObjective(q, oldf, Df, Hf);

    VectorXd perturb(q.size());
    perturb.setZero();

    for(int i=0; i<(int)dirs.size(); i++)
    {
        perturb += dirs[i]*MathUtil::randomDouble(-1.0, 1.0);
    }
    double dt = mag;
    perturb.normalize();
    double f;
    VectorXd newq;
    do
    {
        newq = q + dt * perturb;
        buildObjective(newq, f, Df, Hf);
        dt *= 0.1;
    }
    while(f > oldf);
    q = newq;
    gatherInfo(q);
}

void DevelopableMesh::gatherInfo(const Eigen::VectorXd &q)
{
    double f;
    VectorXd Df;
    vector<T> Hf;
    buildObjective(q, f, Df, Hf);

    cout << "Current objective value: " << f << endl;
    double ecv = equalityConstraintViolation(q);
//    int icv = activeInequalityConstraints(q);

    cout << "Equality constraint violation: " << ecv << endl;
//    cout << icv << " inequality constraints are active" << endl;
    vector<VectorXd> normalspace, tangentspace;
    buildConstraintBasis(q, normalspace, tangentspace);
    cout << "Tangent space is " << tangentspace.size() << " dimensional" << endl;
}
