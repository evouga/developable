#include "developablemesh.h"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Geometry>
#include <set>
#include "coin/IpIpoptApplication.hpp"
#include "devnlp.h"
#include <QGLWidget>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;
using namespace Ipopt;

const double PI = 3.1415926535898;


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


DevelopableMesh::DevelopableMesh() : material_(new MaterialMesh(0))
{

}

DevelopableMesh::~DevelopableMesh()
{
    delete material_;
}

bool DevelopableMesh::loadMesh(const string &)
{
    assert(false);
    return false;
}

void DevelopableMesh::buildSchwarzLantern(double r, double h, int n, int m, double angle)
{
    mesh_ = OMMesh();
    double alpha = angle;
    double a = 2*r*sin(PI/n);
    double b = sqrt(4*r*r*sin(alpha/2.0)*sin(alpha/2.0) + h*h/m/m);
    double c = sqrt(4*r*r*sin(PI/n + alpha/2.0)*sin(PI/n + alpha/2.0) + h*h/m/m);
    double costilt = (a*a+b*b-c*c)/(2*a*b);
    double dpx = -b*costilt;


    double s = 0.5*(a+b+c);
    double area = sqrt(s*(s-a)*(s-b)*(s-c));
    double dpy = 2*area/a;

    double W = n*a;
    double H = 2*n*m*area/W;

    assert(fabs(H-dpy*m) < 1e-6);

    delete material_;
    material_ = new MaterialMesh(H);
    boundaries_.clear();

    Boundary bottom, top;

    for(int i=0; i<=m; i++)
    {
        double z = h*i/m;
        double basepx = i*dpx;
        double py = i*dpy;

        for(int j=0; j<n; j++)
        {
            double x = r*cos(2*PI*(j/double(n)) + i*alpha);
            double y = r*sin(2*PI*(j/double(n)) + i*alpha);

            OMMesh::Point newpt(x,y,z);
            if(i == 0)
            {
                bottom.bdryVerts.push_back(mesh_.n_vertices());
                bottom.bdryPos.push_back(Vector3d(x,y,z));
            }
            else if(i==m)
            {
                top.bdryVerts.push_back(mesh_.n_vertices());
                top.bdryPos.push_back(Vector3d(x,y,z));
            }
            mesh_.add_vertex(newpt);

            double px = j*a+basepx;
            OMMesh::Point newmatpt(px,py,0);
            if(i == 0 || i == m)
            {
                material_->getBoundaryVerts().push_back(pair<int, double>(material_->getMesh().n_vertices(), py));
            }
            material_->getMesh().add_vertex(newmatpt);

        }
    }

    boundaries_.push_back(bottom);
    boundaries_.push_back(top);

    for(int i=0; i<m; i++)
    {
        for(int j=0; j<n; j++)
        {
            int fidx1 = i*n+j;
            int fidx2 = i*n + ((j+1) % n);
            int fidx3 = (i+1)*n+ ((j+1)%n);

            vector<OMMesh::VertexHandle> newface;
            vector<OMMesh::VertexHandle> newmatface;
            newface.push_back(mesh_.vertex_handle(fidx1));
            newface.push_back(mesh_.vertex_handle(fidx2));
            newface.push_back(mesh_.vertex_handle(fidx3));
            newmatface.push_back(material_->getMesh().vertex_handle(fidx1));
            newmatface.push_back(material_->getMesh().vertex_handle(fidx2));
            newmatface.push_back(material_->getMesh().vertex_handle(fidx3));

            mesh_.add_face(newface);
            material_->getMesh().add_face(newmatface);

            if(j == n-1)
            {
                // Wraps around
                Vector3d offset(W,0,0);
                int heid = material_->findHalfedge(fidx1, fidx2);
                assert(heid >= 0);
                material_->addOffset(material_->getMesh().halfedge_handle(heid), offset);
                heid = material_->findHalfedge(fidx1, fidx3);
                assert(heid >= 0);
                material_->addOffset(material_->getMesh().halfedge_handle(heid), offset);
            }

            material_->setMaterialEdge(findEdge(fidx1,fidx2), material_->findEdge(fidx1,fidx2));
            material_->setMaterialEdge(findEdge(fidx2,fidx3), material_->findEdge(fidx2,fidx3));
            material_->setMaterialEdge(findEdge(fidx3,fidx1), material_->findEdge(fidx3,fidx1));

            fidx2 = fidx3;
            newface[1] = mesh_.vertex_handle(fidx2);
            newmatface[1] = material_->getMesh().vertex_handle(fidx2);
            fidx3 = (i+1)*n + j;
            newface[2] = mesh_.vertex_handle(fidx3);
            newmatface[2] = material_->getMesh().vertex_handle(fidx3);
            mesh_.add_face(newface);
            material_->getMesh().add_face(newmatface);

            if(j == n-1)
            {
                // Wraps around
                Vector3d offset(W,0,0);
                int heid = material_->findHalfedge(fidx1, fidx2);
                assert(heid >= 0);
                material_->addOffset(material_->getMesh().halfedge_handle(heid), offset);
                heid = material_->findHalfedge(fidx3, fidx2);
                assert(heid >= 0);
                material_->addOffset(material_->getMesh().halfedge_handle(heid), offset);
            }

            material_->setMaterialEdge(findEdge(fidx1,fidx2), material_->findEdge(fidx1,fidx2));
            material_->setMaterialEdge(findEdge(fidx2,fidx3), material_->findEdge(fidx2,fidx3));
            material_->setMaterialEdge(findEdge(fidx3,fidx1), material_->findEdge(fidx3,fidx1));

        }
    }

    for(int i=0; i<(int)mesh_.n_edges(); i++)
    {
        double len1 = edgeLength(i);

        int medge = material_->materialEdge(i);

        double len2 = material_->edgeLength(medge);

        assert(fabs(len1-len2) < 1e-6);
    }

    for(int i=0; i<(int)material_->getMesh().n_faces(); i++)
    {
        OMMesh::FaceHandle fh = material_->getMesh().face_handle(i);
        for(OMMesh::FaceHalfedgeIter fhi = material_->getMesh().fh_iter(fh); fhi; ++fhi)
        {
            OMMesh::HalfedgeHandle heh = fhi.handle();
            OMMesh::HalfedgeHandle nextheh = material_->getMesh().next_halfedge_handle(heh);
            OMMesh::VertexHandle centv = material_->getMesh().to_vertex_handle(heh);
            OMMesh::VertexHandle v1 = material_->getMesh().from_vertex_handle(heh);
            OMMesh::VertexHandle v2 = material_->getMesh().to_vertex_handle(nextheh);

            OMMesh::Point e0 = material_->getMesh().point(v1) - material_->getMesh().point(centv);
            OMMesh::Point e1 = material_->getMesh().point(v2) - material_->getMesh().point(centv);
            Vector3d ve0 = point2Vector(e0) - material_->getOffset(heh);
            Vector3d ve1 = point2Vector(e1) + material_->getOffset(nextheh);
            Vector3d cr = ve0.cross(ve1);
            assert(cr[2] < 0);
        }
    }

    centerCylinder();
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

void DevelopableMesh::deformLantern(DeformCallback &dc)
{
    bool keepgoing = true;
    while(keepgoing)
    {
        int nembverts = mesh_.n_vertices();
        int nmatverts = material_->getMesh().n_vertices();

        VectorXd q(3*nembverts+2*nmatverts);
        for(int i=0; i<nembverts; i++)
        {
            for(int j=0; j<3; j++)
            {
                q[3*i+j] = mesh_.point(mesh_.vertex_handle(i))[j];
            }
        }
        for(int i=0; i<nmatverts; i++)
        {
            for(int j=0; j<2; j++)
            {
                q[3*nembverts+2*i+j] = material_->getMesh().point(material_->getMesh().vertex_handle(i))[j];
            }
        }

        cout << "Initial equality violation: " << equalityConstraintViolation(q) << endl;

        int inequalities = activeInequalityConstraints(q);
        cout << "Active inequality constraints: " << inequalities << endl;
        double f;
        VectorXd Df;
        vector<T> Hf;
        buildObjective(q, f, Df, Hf);
        cout << "Initial energy: " << f << endl;
        DevTLNP *mynlp = new DevTLNP(*this, dc, q);
        IpoptApplication app;
        app.Options()->SetNumericValue("tol", 1e-6);
        //app.Options()->SetStringValue("mu_strategy", "adaptive");
        //app.Options()->SetStringValue("derivative_test", "second-order");
        app.Options()->SetStringValue("check_derivatives_for_naninf", "yes");
        //app.Options()->SetIntegerValue("print_level", 12);
        app.Initialize();
        ApplicationReturnStatus status = app.OptimizeTNLP(mynlp);

        keepgoing = (status == User_Requested_Stop);
        q = mynlp->getFinalQ();

        flushOutNANs(q);

        inequalities = activeInequalityConstraints(q);
        cout << "Initial equality violation: " << equalityConstraintViolation(q) << endl;
        cout << "Active inequality constraints: " << inequalities << endl;

        checkConstrainedHessian(q);

        int collapseid = findCollapsibleEdge(q);

        repopulateDOFs(q);

        if(collapseid != -1)
        {
            assert(canCollapseEdge(collapseid));
            collapseEdge(collapseid);
        }

        centerCylinder();
    }
}

bool DevelopableMesh::shouldCollapseEdges(const Eigen::VectorXd &q)
{
    return (findCollapsibleEdge(q) != -1);
}

int DevelopableMesh::findCollapsibleEdge(const Eigen::VectorXd &q)
{
    const double constrtol = 0.01;

    for(OMMesh::HalfedgeIter hei = mesh_.halfedges_begin(); hei != mesh_.halfedges_end(); ++hei)
    {
        if(mesh_.is_boundary(hei.handle()))
            continue;
        double g;
        vector<pair<int, double> > dummy1;
        vector<T> dummy2;

        int v1 = mesh_.from_vertex_handle(hei.handle()).idx();
        int v2 = mesh_.to_vertex_handle(hei.handle()).idx();

        int hidx = material_->findHalfedge(v1, v2);
        assert(hidx != -1);

        OMMesh::HalfedgeHandle heh1 = material_->getMesh().halfedge_handle(hidx);

        radiusOverlapConstraint(q, heh1, g, dummy1, dummy2);

        if(fabs(g) < constrtol)
        {
            int candidate = mesh_.edge_handle(mesh_.prev_halfedge_handle(hei.handle())).idx();
            if(canCollapseEdge(candidate))
                return candidate;
        }
    }

    return -1;
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
        for(OMMesh::VertexOHalfedgeIter voh = material_->getMesh().voh_iter(vh); voh; ++voh)
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
    glColor3f(0.0,0.0,0.0);
    glBegin(GL_LINES);
    for(int i=0; i<(int)material_->getMesh().n_halfedges(); i++)
    {
        OMMesh::HalfedgeHandle heh = material_->getMesh().halfedge_handle(i);
        OMMesh::Point pt1 = material_->getMesh().point(material_->getMesh().from_vertex_handle(heh));
        OMMesh::Point pt2 = material_->getMesh().point(material_->getMesh().to_vertex_handle(heh));
        Vector3d offset = material_->getOffset(heh);
        glVertex2d(pt1[0], pt1[1]);
        glVertex2d(pt2[0] + offset[0], pt2[1] + offset[1]);
    }
    glEnd();

    glPointSize(6.0);
    glBegin(GL_POINTS);
    for(int i=0; i<(int)material_->getMesh().n_vertices(); i++)
    {
        glColor3f(0.2, 0.0, 1.0);
        OMMesh::VertexHandle vh = material_->getMesh().vertex_handle(i);
        OMMesh::Point pt = material_->getMesh().point(vh);
        glVertex2d(pt[0], pt[1]);
    }
    glColor3f(1.0, 0.0, 0.0);
    for(int i=0; i<(int)material_->getBoundaryVerts().size(); i++)
    {
        OMMesh::VertexHandle vh = material_->getMesh().vertex_handle(material_->getBoundaryVerts()[i].first);
        OMMesh::Point pt = material_->getMesh().point(vh);
        glVertex2d(pt[0], pt[1]);
    }
    glEnd();
    glPointSize(1.0);
}

void DevelopableMesh::collapseEdge(int eid)
{
    OMMesh::EdgeHandle eh = mesh_.edge_handle(eid);
    OMMesh::HalfedgeHandle heh = mesh_.halfedge_handle(eh,0);
    if(mesh_.is_boundary(mesh_.from_vertex_handle(heh)))
        heh = mesh_.opposite_halfedge_handle(heh);
    assert(!mesh_.is_boundary(mesh_.from_vertex_handle(heh)));

    int v2 = mesh_.to_vertex_handle(heh).idx();
    int v1 = mesh_.from_vertex_handle(heh).idx();

    int mheid = material_->findHalfedge(v1, v2);
    OMMesh::HalfedgeHandle mheh = material_->getMesh().halfedge_handle(mheid);
    Vector3d deletedoffset = material_->getOffset(mheh);

    vector<int> old2new;

    OMMesh newmesh;
    MaterialMesh *newmmesh = new MaterialMesh(material_->getH());

    int mergedid = -1;
    int newid = 0;
    for(int i=0; i<(int)mesh_.n_vertices(); i++)
    {
        if(i == v1)
        {
            old2new.push_back(-1);
            continue;
        }
        if(i == v2)
        {
            mergedid = newid;
        }

        old2new.push_back(newid++);
        OMMesh::Point pt = mesh_.point(mesh_.vertex_handle(i));
        newmesh.add_vertex(pt);
        OMMesh::Point mpt = material_->getMesh().point(material_->getMesh().vertex_handle(i));
        newmmesh->getMesh().add_vertex(mpt);
    }

    old2new[v1] = mergedid;

    // Fix boundary data
    vector<Boundary> newbdries;
    for(int i=0; i<(int)boundaries_.size(); i++)
    {
        bool seen = false;
        Boundary newbd;
        for(int j=0; j<(int)boundaries_[i].bdryVerts.size(); j++)
        {
            if(boundaries_[i].bdryVerts[j] == v1 || boundaries_[i].bdryVerts[j] == v2)
            {
                if(seen)
                    continue;
                else
                    seen = true;
            }
            newbd.bdryPos.push_back(boundaries_[i].bdryPos[j]);
            newbd.bdryVerts.push_back(old2new[boundaries_[i].bdryVerts[j]]);
        }
        newbdries.push_back(newbd);
    }

    vector<pair<int, double> > newBdryVerts;
    bool seen = false;
    for(int i=0; i<(int)material_->getBoundaryVerts().size(); i++)
    {
        pair<int, double> entry = material_->getBoundaryVerts()[i];
        if(entry.first == v1 || entry.first == v2)
        {
            if(seen)
                continue;
            else
                seen = true;
        }
        newBdryVerts.push_back(pair<int, double>(old2new[entry.first], entry.second));
        newmmesh->getBoundaryVerts() = newBdryVerts;
    }

    // ignored faces
    int f1 = mesh_.face_handle(heh).idx();
    int f2 = mesh_.opposite_face_handle(heh).idx();

    for(int i=0; i<(int)mesh_.n_faces(); i++)
    {
        if(i == f1 || i == f2)
            continue;
        OMMesh::FaceHandle fh = mesh_.face_handle(i);
        vector<OMMesh::VertexHandle> newface;
        for(OMMesh::FaceVertexIter fvi = mesh_.fv_iter(fh); fvi; ++fvi)
        {
            int vid = fvi.handle().idx();
            OMMesh::VertexHandle vh = newmesh.vertex_handle(old2new[vid]);
            newface.push_back(vh);
        }
        newmesh.add_face(newface);
        OMMesh::FaceHandle mfh = material_->getMesh().face_handle(i);
        vector<OMMesh::VertexHandle> newmface;
        for(OMMesh::FaceVertexIter fvi = material_->getMesh().fv_iter(mfh); fvi; ++fvi)
        {
            int vid = fvi.handle().idx();
            OMMesh::VertexHandle vh = newmmesh->getMesh().vertex_handle(old2new[vid]);
            newmface.push_back(vh);
        }
        newmmesh->getMesh().add_face(newmface);
    }

    // Reset edge map
    for(OMMesh::EdgeIter ei = newmesh.edges_begin(); ei != newmesh.edges_end(); ++ei)
    {
        OMMesh::HalfedgeHandle heh = newmesh.halfedge_handle(ei.handle(), 0);
        int v1 = newmesh.from_vertex_handle(heh).idx();
        int v2 = newmesh.to_vertex_handle(heh).idx();
        int meid = newmmesh->findEdge(v1,v2);
        newmmesh->setMaterialEdge(ei.handle().idx(), meid);
    }

    // set edge offset (except v2's edges)
    for(OMMesh::HalfedgeIter hei = material_->getMesh().halfedges_begin(); hei != material_->getMesh().halfedges_end(); ++hei)
    {
        int oldv1 = material_->getMesh().from_vertex_handle(hei.handle()).idx();
        int oldv2 = material_->getMesh().to_vertex_handle(hei.handle()).idx();
        if(oldv1 == v2 || oldv2 == v2)
            continue;
        int newv1 = old2new[oldv1];
        int newv2 = old2new[oldv2];
        int newheid = newmmesh->findHalfedge(newv1,newv2);
        assert(newheid != -1);
        OMMesh::HalfedgeHandle heh = newmmesh->getMesh().halfedge_handle(newheid);
        newmmesh->addOffset(heh, material_->getOffset(hei.handle()));
    }

    // last order of business: fix the offsets of the edges coming off of v2
    OMMesh::VertexHandle source = material_->getMesh().vertex_handle(v2);
    for(OMMesh::VertexOHalfedgeIter voh = material_->getMesh().voh_iter(source); voh; ++voh)
    {
        int target = material_->getMesh().to_vertex_handle(voh.handle()).idx();
        if(target == v1)
            continue;
        Vector3d newoffset = deletedoffset + material_->getOffset(voh.handle());
        int newheid = newmmesh->findHalfedge(old2new[v1], old2new[target]);
        assert(newheid != -1);
        OMMesh::HalfedgeHandle newheh = newmmesh->getMesh().halfedge_handle(newheid);
        newmmesh->addOffset(newheh, newoffset);
    }

    mesh_ = newmesh;
    delete material_;
    material_ = newmmesh;
    boundaries_ = newbdries;
}

bool DevelopableMesh::canCollapseEdge(int eid)
{
    // Some sanity checks. First, the edge should not span the boundaries:
    OMMesh::EdgeHandle eh = mesh_.edge_handle(eid);
    OMMesh::HalfedgeHandle heh = mesh_.halfedge_handle(eh,0);
    OMMesh::VertexHandle source = mesh_.from_vertex_handle(heh);
    OMMesh::VertexHandle target = mesh_.to_vertex_handle(heh);

    if(mesh_.is_boundary(source) && mesh_.is_boundary(target))
        return false;

    // Next, the source and target verts should have only two other verts in common (otherwise collapsed mesh is nonmanifold)
    set<int> neighbors1;
    for(OMMesh::VertexVertexIter vvi = mesh_.vv_iter(source); vvi; ++vvi)
        neighbors1.insert(vvi.handle().idx());

    int duplicates=0;
    for(OMMesh::VertexVertexIter vvi = mesh_.vv_iter(target); vvi; ++vvi)
    {
        if(neighbors1.count(vvi.handle().idx()) > 0)
            duplicates++;
    }
    if(duplicates != 2)
        return false;

    // we're ok
    return true;
}

void DevelopableMesh::repopulateDOFs(const Eigen::VectorXd &q)
{
    int nembverts = mesh_.n_vertices();
    int nmatverts = material_->getMesh().n_vertices();

    for(int i=0; i<nembverts; i++)
    {
        for(int j=0; j<3; j++)
        {
            mesh_.point(mesh_.vertex_handle(i))[j] = q[3*i+j];
        }
    }

    for(int i=0; i<nmatverts; i++)
    {
        for(int j=0; j<2; j++)
        {
            material_->getMesh().point(material_->getMesh().vertex_handle(i))[j] = q[3*nembverts+2*i+j];
        }
    }

    centerCylinder();
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

void DevelopableMesh::checkConstrainedHessian(const Eigen::VectorXd &q)
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
    sparse.setFromTriplets(Hf.begin(), Hf.end());
    MatrixXd constrainedH = Tan*sparse*Tan.transpose();

    SelfAdjointEigenSolver<MatrixXd> es;
    es.compute(constrainedH);
    VectorXd evs = es.eigenvalues();
    double mostneg = -1e-6;
    double mostpos = 1e-6;
    double startneg = mostneg;
    double startpos = mostpos;
    for(int i=0; i<(int)evs.size(); i++)
    {
        if(evs[i] < mostneg)
            mostneg = evs[i];
        if(evs[i] > mostpos)
            mostpos = evs[i];
    }
    if(mostneg < startneg && mostpos > startpos)
        std::cout << "Saddle Point, " << mostneg << " " << mostpos << std::endl;
    else if(mostneg < startneg)
        std::cout << "Maximum" << std::endl;
    else if(mostpos > startpos)
        std::cout << "Minimum" << std::endl;
    else
        std::cout << "Singular" << std::endl;
}

bool DevelopableMesh::hasNANs(const Eigen::VectorXd &v)
{
    for(int i=0; i<(int)v.size(); i++)
        if(isnan(v[i]))
            return true;
    return false;
}

bool DevelopableMesh::hasNANs(const std::vector<T> &M)
{
    for(int i=0; i<(int)M.size(); i++)
        if(isnan(M[i].value()))
            return true;
    return false;
}

bool DevelopableMesh::hasNANs(const std::vector<std::vector<T> > &H)
{
    for(int i=0; i<(int)H.size(); i++)
        for(int j=0; j<(int)H[i].size(); j++)
        {
            if(isnan(H[i][j].value()))
                return true;
        }
    return false;
}

void DevelopableMesh::flushOutNANs(const Eigen::VectorXd &q)
{
    assert(!hasNANs(q));
    double f;
    VectorXd Df;
    vector<T> Hf;
    buildObjective(q, f, Df, Hf);
    assert(!isnan(f));
    assert(!hasNANs(Df));
    assert(!hasNANs(Hf));

    VectorXd g;
    vector<T> Dg;
    vector<vector<T> > Hg;
    buildConstraints(q, g, Dg, Hg);

    assert(!hasNANs(g));
    assert(!hasNANs(Dg));
    assert(!hasNANs(Hg));
}

double DevelopableMesh::equalityConstraintViolation(const Eigen::VectorXd &q)
{
    VectorXd g;
    vector<T> Dg;
    vector<vector<T> > Hg;
    buildConstraints(q, g, Dg, Hg);
    return g.norm();
}
