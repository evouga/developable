#include "developablemesh.h"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Geometry>
#include <set>
#include <fadbad.h>
#include <fadiff.h>

using namespace std;

using namespace Eigen;
using namespace fadbad;
using namespace Ipopt;

const double PI = 3.1415926535898;

DevelopableMesh::DevelopableMesh()
{

}

bool DevelopableMesh::loadMesh(const string &filename)
{
    if(!Mesh::loadMesh(filename))
        return false;
    identifyBoundaries();
    calculateSurfaceArea();
    return true;
}

void DevelopableMesh::buildSchwarzLantern(double r, double h, int n, int m)
{
    mesh_ = OMMesh();

    for(int i=0; i<m; i++)
    {
        double z = h*i/(m-1);
        for(int j=0; j<n; j++)
        {
            double x = r*cos(2*PI*(j/double(n) + (i%2)/double(2*n)));
            double y = r*sin(2*PI*(j/double(n) + (i%2)/double(2*n)));

            OMMesh::Point newpt(x,y,z);

            mesh_.add_vertex(newpt);
        }
    }
    for(int i=0; i<m-1; i++)
    {
        for(int j=0; j<n; j++)
        {
            int fidx1 = i*n+j;
            int fidx2 = i*n + ((j+1) % n);
            int fidx3 = (i+1)*n+ ((j + i%2) % n);
            vector<OMMesh::VertexHandle> newface;
            newface.push_back(mesh_.vertex_handle(fidx1));
            newface.push_back(mesh_.vertex_handle(fidx2));
            newface.push_back(mesh_.vertex_handle(fidx3));
            mesh_.add_face(newface);

            fidx2 = fidx3;
            newface[1] = mesh_.vertex_handle(fidx2);
            fidx3 = (i+1)*n + ((n + j-1 + i%2) % n);
            newface[2] = mesh_.vertex_handle(fidx3);
            mesh_.add_face(newface);

        }
    }

    identifyBoundaries();
    calculateSurfaceArea();
}

void DevelopableMesh::identifyBoundaries()
{
    boundaries_.clear();
    int numedges = mesh_.n_edges();
    bool *visited = new bool[numedges];

    for(int i=0; i<numedges; i++)
        visited[i] = false;


    for(int i=0; i<numedges; i++)
    {
        if(visited[i])
            continue;

        if(mesh_.is_boundary(mesh_.edge_handle(i)))
        {
            BoundaryCurve bc;
            bc.arclength = 0;
            OMMesh::EdgeHandle eh = mesh_.edge_handle(i);
            OMMesh::HalfedgeHandle heh = mesh_.halfedge_handle(eh,0);
            if(!mesh_.is_boundary(heh))
                heh = mesh_.opposite_halfedge_handle(heh);
            bc.height = mesh_.point(mesh_.from_vertex_handle(heh))[2];
            bc.targetheight = bc.height;
            while(!visited[eh.idx()])
            {
                visited[eh.idx()] = true;
                bc.edges.push_back(eh.idx());
                bc.arclength += mesh_.calc_edge_length(eh);

                heh = mesh_.next_halfedge_handle(heh);
                eh = mesh_.edge_handle(heh);

                assert(mesh_.point(mesh_.from_vertex_handle(heh))[2] == bc.height);
            }

            boundaries_.push_back(bc);

        }

    }

    delete[] visited;
}

void DevelopableMesh::calculateSurfaceArea()
{
    surfacearea_ = 0;
    for(int i=0; i< (int)mesh_.n_vertices(); i++)
    {
        OMMesh::VertexHandle vh = mesh_.vertex_handle(i);

        int valence = mesh_.valence(vh);
        OMMesh::VertexOHalfedgeIter voh = mesh_.voh_iter(vh);
        for(int i=0; i<valence; i++)
        {
            OMMesh::Point e1pt = mesh_.point(mesh_.to_vertex_handle(voh)) - mesh_.point(vh);
            Vector3d e1(e1pt[0], e1pt[1], e1pt[2]);

            ++voh;
            if(!voh)
                voh = mesh_.voh_iter(vh);

            if(mesh_.is_boundary(voh))
                continue;

            OMMesh::Point e2pt = mesh_.point(mesh_.to_vertex_handle(voh)) - mesh_.point(vh);
            Vector3d e2(e2pt[0], e2pt[1], e2pt[2]);
            surfacearea_ += 0.5/3.0 * (e1.cross(e2)).norm();
        }
    }
}

void DevelopableMesh::getBoundaryHeights(std::vector<double> &heights)
{
    heights.clear();
    for(int i=0; i<(int)boundaries_.size(); i++)
        heights.push_back(boundaries_[i].height);
}

void DevelopableMesh::deformLantern(int maxiters)
{
    int numverts = mesh_.n_vertices();

    int numedges = mesh_.n_edges();
    int interiore = 0;
    for(int i=0; i<numedges; i++)
    {
        OMMesh::EdgeHandle eh = mesh_.edge_handle(i);
        if(!mesh_.is_boundary(eh))
            interiore++;
    }

    int numdofs = 3*numverts;    

    VectorXd q(numdofs);
    for(int i=0; i<numverts; i++)
    {
        OMMesh::Point &pt = mesh_.point(mesh_.vertex_handle(i));
        for(int k=0; k<3; k++)
            q[3*i+k] = pt[k];
    }

    double f;
    VectorXd Df;
    vector<T> Hf;
    VectorXd g;
    vector<T> Dg;
    vector<vector<T> > Hg;

    cout << "build objective" << endl;
    buildObjective(q, f, Df, Hf);
    cout << "build constraints" << endl;
    buildConstraints(q, g, Dg, Hg);
    cout << "done" << endl;

    double nnz_j = Dg.size();
    double nnz_h = Hf.size();
    for(int i=0; i<Hg.size(); i++)
        nnz_h += Hg[i].size();

    double numconstraints = g.size();

    SmartPtr<TNLP> mynlp = new IpoptSolver(numdofs, numconstraints, nnz_j, nnz_h, q, *this);

    cout << "Creating app" << endl;
    IpoptApplication *app = new IpoptApplication();

     // Change some options
     // Note: The following choices are only examples, they might not be
     //       suitable for your optimization problem.
     app->Options()->SetNumericValue("tol", 1e-9);
     app->Options()->SetStringValue("mu_strategy", "adaptive");
     //app->Options()->SetStringValue("derivative_test", "second-order");
     //app->Options()->SetStringValue("output_file", "ipopt.out");

     // Intialize the IpoptApplication and process the options
     ApplicationReturnStatus status;
     status = app->Initialize();
     if (status != Solve_Succeeded) {
       cout << "\n\n*** Error during initialization!\n" << endl;
       return;
     }
     cout << "Starting solve" << endl;

     // Ask Ipopt to solve the problem
     status = app->OptimizeTNLP(mynlp);
}

Vector3d DevelopableMesh::point2Vector(OMMesh::Point pt)
{
    Vector3d result;
    for(int i=0; i<3; i++)
        result[i] = pt[i];
    return result;
}

void DevelopableMesh::buildConstraints(const VectorXd &q, VectorXd &g, vector<T> &Dg, vector<vector<T> > &Hg)
{
    int numverts = mesh_.n_vertices();
    int numdofs = 3*numverts;
    int boundaryverts = 0;
    Dg.clear();
    Hg.clear();

    for(int i=0; i<numverts; i++)
    {
        OMMesh::VertexHandle vh = mesh_.vertex_handle(i);
        if(mesh_.is_boundary(vh))
            boundaryverts++;
    }

    // sum of interior vertex angles = 2pi
    // sum of boundary vertex angles = pi
    // arc length of boundaries stays fixed
    // heights of boundaries vertices are fixed
    // total area is fixed
    int numconstraints = numverts + 3*boundaryverts + boundaries_.size() + 1;

    // build g and Dg
    int row = 0;
    g.resize(numconstraints);
    g.setZero();

    // interior angles sum to 2 PI, boundary angles sum to PI
    for(int i=0; i<numverts; i++)
    {
        OMMesh::VertexHandle vh = mesh_.vertex_handle(i);
        bool bdry = mesh_.is_boundary(vh);

        F<F<double> > deficit = 0;

        int valence = mesh_.valence(vh);

        vector<int> vidxs;
        for(int j=0; j<=valence; j++)
        {
            vidxs.push_back(0);
        }

        vidxs[valence] = i;
        F<F<double> > centpt[3];
        F<double> centptd[3];
        for(int k=0; k<3; k++)
        {
            centptd[k] = q[3*vh.idx()+k];
            centptd[k].diff(3*valence+k, 3*(valence+1));
            centpt[k] = centptd[k];
            centpt[k].diff(3*valence+k, 3*(valence+1));
        }

        OMMesh::VertexOHalfedgeIter voh = mesh_.voh_iter(vh);
        for(int j=0; j<valence; j++)
        {
            OMMesh::VertexHandle tov = mesh_.to_vertex_handle(voh);
            vidxs[j] = tov.idx();
            F<F<double> > tovpt[3];
            F<double> tovptd[3];
            for(int k=0; k<3; k++)
            {
                tovptd[k] = q[3*tov.idx()+k];
                tovptd[k].diff(3*j+k, 3*(valence+1));
                tovpt[k] = tovptd[k];
                tovpt[k].diff(3*j+k, 3*(valence+1));
            }
            F<F<double> > e1[3];
            for(int k=0; k<3; k++)
            {
                e1[k] = tovpt[k] - centpt[k];
            }

            ++voh;
            if(!voh)
                voh = mesh_.voh_iter(vh);

            // Ignore boundary "faces"
            if(mesh_.is_boundary(voh))
                continue;


            //OMMesh::Point e2pt = mesh_.point(mesh_.to_vertex_handle(voh)) - mesh_.point(vh);
            F<F<double> > tovpt2[3];
            F<double> tovpt2d[3];
            for(int k=0; k<3; k++)
            {
                tovpt2d[k] = q[3*mesh_.to_vertex_handle(voh).idx()+k];
                tovpt2d[k].diff( 3*((j+1)%valence) + k, 3*(valence+1) );
                tovpt2[k] = tovpt2d[k];
                tovpt2[k].diff( 3*((j+1)%valence) + k, 3*(valence+1) );
            }
            F<F<double> > e2[3];
            for(int k=0; k<3; k++)
            {
                e2[k] = tovpt2[k]-centpt[k];
            }
            F<F<double> > e1norm = norm(e1);
            F<F<double> > e2norm = norm(e2);
            F<F<double> > theta = acos(dot(e1,e2) / e1norm / e2norm );
            deficit += theta;
        }
        deficit -= (bdry ? PI : 2*PI);
        g[row] = deficit.val().val();

        vector<T> Hgentry;

        for(int j=0; j<=valence; j++)
        {
            for(int k=0; k<3; k++)
            {
                Dg.push_back(T(row, 3*vidxs[j]+k, deficit.d(3*j+k).val()));
                for(int l=0; l<=valence; l++)
                {
                    for(int m=0; m<3; m++)
                    {
                        if(3*j+k <= 3*l+m)
                            Hgentry.push_back(T(3*vidxs[j]+k, 3*vidxs[l]+m, deficit.d(3*j+k).d(3*l+m)));
                    }
                }
            }
        }
        Hg.push_back(Hgentry);
        row++;
    }
    assert(row == numverts);

    // boundary arc length is constant
    for(int i=0; i<(int)boundaries_.size(); i++)
    {
        F<F<double> > len = 0;

        int numedges = boundaries_[i].edges.size();
        map<int, int> vertids;
        int curid=0;

        for(int j=0; j<numedges; j++)
        {
            //double curlen = mesh_.calc_edge_length(mesh_.edge_handle(boundaries_[i].edges[j]));
            //len += curlen;
            OMMesh::HalfedgeHandle heh = mesh_.halfedge_handle(mesh_.edge_handle(boundaries_[i].edges[j]), 0);
            int idx1 = mesh_.to_vertex_handle(heh).idx();
            int idx2 = mesh_.from_vertex_handle(heh).idx();
            if(vertids.find(idx1) == vertids.end())
                vertids[idx1] = curid++;
            if(vertids.find(idx2) == vertids.end())
                vertids[idx2] = curid++;

            F<F<double> > pt1[3];
            F<F<double> > pt2[3];
            F<double> pt1d[3];
            F<double> pt2d[3];
            for(int k=0; k<3; k++)
            {
                pt1d[k] = q[3*idx1+k];
                pt1d[k].diff(3*vertids[idx1]+k, 3*numedges);
                pt1[k] = pt1d[k];
                pt1[k].diff(3*vertids[idx1]+k, 3*numedges);
                pt2d[k] = q[3*idx2+k];
                pt2d[k].diff(3*vertids[idx2]+k, 3*numedges);
                pt2[k] = pt2d[k];
                pt2[k].diff(3*vertids[idx2]+k, 3*numedges);
            }

            F<F<double> > edge[3];
            for(int k=0; k<3; k++)
            {
                edge[k] = pt2[k] - pt1[k];
            }

            len += norm(edge);

        }

        len -= boundaries_[i].arclength;
        g[row] = len.val().val();

        vector<T> Hgentry;
        for(map<int,int>::iterator it = vertids.begin(); it != vertids.end(); ++it)
        {
            for(int k=0; k<3; k++)
            {
                Dg.push_back(T(row, 3*it->first+k, len.d(3*it->second+k).val()));
                for(map<int,int>::iterator it2 = vertids.begin(); it2 != vertids.end(); ++it2)
                {
                    for(int l=0; l<3; l++)
                    {
                        if(3*it->first+k <= 3*it2->first+l)
                            Hgentry.push_back(T(3*it->first+k, 3*it2->first+l, len.d(3*it->second+k).d(3*it2->second+l)));
                    }
                }
            }
        }
        Hg.push_back(Hgentry);

        row++;
    }

    assert(row == numverts + (int)boundaries_.size());

    // boundary vertices have fixed height
    for(int i=0; i<(int)boundaries_.size(); i++)
    {
        set<int> bdverts;
        for(int j=0; j<(int)boundaries_[i].edges.size(); j++)
        {
            OMMesh::HalfedgeHandle heh = mesh_.halfedge_handle(mesh_.edge_handle(boundaries_[i].edges[j]),0);
            bdverts.insert(mesh_.from_vertex_handle(heh).idx());
            bdverts.insert(mesh_.to_vertex_handle(heh).idx());
        }
        for(set<int>::iterator it = bdverts.begin(); it != bdverts.end(); ++it)
        {
            for(int k=0;k<3;k++)
            {
                OMMesh::Point pt = mesh_.point(mesh_.vertex_handle(*it));
                g[row] = q[3* *it+k] - pt[k];//boundaries_[i].targetheight;
                Dg.push_back(T(row, 3* *it + k, 1.0));
                vector<T> empty;
                Hg.push_back(empty);
                row++;
            }
        }
    }

    assert(row == numverts + (int)boundaries_.size() + 3*boundaryverts);

    // total area is constant
    F<F<double> > area = 0;

    for(int i=0; i<numverts; i++)
    {
        OMMesh::VertexHandle vh = mesh_.vertex_handle(i);
        F<F<double> > centpt[3];
        F<double> centptd[3];
        for(int k=0; k<3; k++)
        {
            centptd[k] = q[3*vh.idx()+k];
            centptd[k].diff(3*i+k, 3*numverts);
            centpt[k] = centptd[k];
            centpt[k].diff(3*i+k, 3*numverts);
        }
        int valence = mesh_.valence(vh);
        OMMesh::VertexOHalfedgeIter voh = mesh_.voh_iter(vh);
        for(int j=0; j<valence; j++)
        {
            OMMesh::VertexHandle tov = mesh_.to_vertex_handle(voh);
            int idx1 = tov.idx();
            F<F<double> > tovpt[3];
            F<double> tovptd[3];
            for(int k=0; k<3; k++)
            {
                tovptd[k] = q[3*tov.idx() + k];
                tovptd[k].diff(3*idx1+k, 3*numverts);
                tovpt[k] = tovptd[k];
                tovpt[k].diff(3*idx1+k, 3*numverts);
            }
            F<F<double> > e1[3];
            for(int k=0; k<3; k++)
                e1[k] = tovpt[k]-centpt[k];

            ++voh;
            if(!voh)
                voh = mesh_.voh_iter(vh);

            if(mesh_.is_boundary(voh))
                continue;

            tov = mesh_.to_vertex_handle(voh);
            int idx2 = tov.idx();
            F<F<double> > tovpt2[3];
            F<double> tovpt2d[3];
            for(int k=0; k<3; k++)
            {
                tovpt2d[k] = q[3*tov.idx()+k];
                tovpt2d[k].diff(3*idx2+k, 3*numverts);
                tovpt2[k] = tovpt2d[k];
                tovpt2[k].diff(3*idx2+k, 3*numverts);
            }
            F<F<double> > e2[3];
            for(int k=0; k<3; k++)
                e2[k] = tovpt2[k]-centpt[k];

            F<F<double> > cr[3];
            cross(e1,e2,cr);
            area += 0.5/3.0 * norm(cr);
        }
    }

    area -= surfacearea_;
    g[row] = area.val().val();
    vector<T> Hgpart;

    for(int i=0; i<numverts; i++)
    {
        for(int k=0; k<3; k++)
        {
            Dg.push_back(T(row, 3*i+k, area.d(3*i+k).val()));
            for(int j=0; j<numverts; j++)
            {
                for(int l=0; l<3; l++)
                {
                    if(3*i+k <= 3*j+l)
                        Hgpart.push_back(T(3*i+k, 3*j+l, area.d(3*i+k).d(3*j+l)));
                }
            }
        }
    }
    Hg.push_back(Hgpart);
    row++;
    assert(Hg.size() == numconstraints);
    assert(row == numconstraints);
}

void DevelopableMesh::buildObjective(const VectorXd &q, double &f, Eigen::VectorXd &Df, vector<T> &Hf)
{
    int numverts = mesh_.n_vertices();
    int numdofs = 3*numverts;

    f = 0;
    Df.resize(numdofs);
    Df.setZero();
    map<pair<int, int>, double> nonzeros;

    for(int i=0; i<(int)mesh_.n_edges(); i++)
    {
        if(mesh_.is_boundary(mesh_.edge_handle(i)))
            continue;

        F<F<double> > Ff = 0;

        OMMesh::HalfedgeHandle heh = mesh_.halfedge_handle(mesh_.edge_handle(i), 0);
        int inds[4];
        F<F<double> > x[4][3];
        F<double> xx[4][3];

        inds[0] = mesh_.from_vertex_handle(heh).idx();
        inds[1] = mesh_.to_vertex_handle(heh).idx();
        OMMesh::HalfedgeHandle nextheh = mesh_.next_halfedge_handle(heh);
        inds[2] = mesh_.to_vertex_handle(nextheh).idx();
        OMMesh::HalfedgeHandle prevheh = mesh_.next_halfedge_handle(mesh_.opposite_halfedge_handle(heh));
        inds[3] = mesh_.to_vertex_handle(prevheh).idx();

        for(int j=0; j<4; j++)
            for(int k=0; k<3; k++)
            {
                xx[j][k] = q[3*inds[j]+k];
                xx[j][k].diff(3*j+k, 12);
                x[j][k] = xx[j][k];
                x[j][k].diff(3*j+k, 12);
            }

        F<F<double> > e0[3], e1[3], e2[3];
        for(int k=0; k<3; k++)
        {
            e0[k] = x[1][k]-x[0][k];
            e1[k] = x[2][k]-x[0][k];
            e2[k] = x[3][k]-x[0][k];
        }

        F< F<double> > e0n = norm(e0);

        F<F<double> > n0[3], n1[3];

        cross(e0, e1, n0);
        cross(e2, e0, n1);

        normalize(n0);
        normalize(n1);

        F<F<double> > theta = acos(dot(n0,n1)-1e-10);

        Ff = pow(e0n, 1.0/3.0) * pow(theta, 7.0/3.0);

        //cout << Ff.val().val() << endl;

        // NAN alert
        if(theta.val().val() == 0.0)
            continue;

        f += Ff.val().val();

        for(int j=0; j<4; j++)
        {
            for(int k=0; k<3; k++)
            {
                F<double> ideriv = Ff.d(3*j+k);
                Df[3*inds[j]+k] += ideriv.val();

                for(int l=0; l<4; l++)
                {
                    for(int m=0; m<3; m++)
                        if(3*inds[j]+k <= 3*inds[l]+m)
                            nonzeros[pair<int,int>(3*inds[j]+k,3*inds[l]+m)] += ideriv.d(3*l+m);
                }
            }
        }
    }

    for(map<pair<int, int>, double>::iterator it = nonzeros.begin(); it != nonzeros.end(); ++it)
    {
        Hf.push_back(T(it->first.first, it->first.second, it->second));
    }
}

bool DevelopableMesh::canCollapseEdge(int edgeidx)
{
    OMMesh::EdgeHandle eh = mesh_.edge_handle(edgeidx);
    bool bdry = false;
    set<int> ignore;

    for(int i=0; i<2; i++)
    {
        OMMesh::HalfedgeHandle heh = mesh_.halfedge_handle(eh, i);
        OMMesh::HalfedgeHandle flip = mesh_.opposite_halfedge_handle(heh);
        if(mesh_.is_boundary(heh))
        {
            if(mesh_.is_boundary(flip))
                return false;
            if(mesh_.is_boundary(mesh_.next_halfedge_handle(flip)))
                return false;
            if(mesh_.is_boundary(mesh_.prev_halfedge_handle(flip)))
                return false;

            bdry = true;
        }
        else
            ignore.insert(mesh_.to_vertex_handle(mesh_.next_halfedge_handle(heh)).idx());
    }

    if(!bdry)
    {
        for(int i=0; i<2; i++)
        {
            OMMesh::HalfedgeHandle heh = mesh_.halfedge_handle(eh, i);
            OMMesh::HalfedgeHandle flip = mesh_.opposite_halfedge_handle(heh);


            if(mesh_.is_boundary(mesh_.next_halfedge_handle(heh)))
            {
                if(mesh_.is_boundary(mesh_.next_halfedge_handle(flip))
                        || mesh_.is_boundary(mesh_.prev_halfedge_handle(heh)))
                {
                    return false;
                }
            }
        }
    }

    OMMesh::HalfedgeHandle heh = mesh_.halfedge_handle(eh,0);
    OMMesh::VertexHandle v1 = mesh_.from_vertex_handle(heh);
    set<int> adj1;
    for(OMMesh::VertexVertexIter vvi = mesh_.vv_iter(v1); vvi; ++vvi)
    {
        adj1.insert(vvi.handle().idx());
    }
    OMMesh::VertexHandle v2 = mesh_.to_vertex_handle(heh);
    for(OMMesh::VertexVertexIter vvi = mesh_.vv_iter(v2); vvi; ++vvi)
    {
        int idx = vvi.handle().idx();
        if(adj1.count(idx) > 0 && ignore.count(idx) == 0)
        {
            return false;
        }
    }
    if(bdry)
        return false;
    return true;
}

void DevelopableMesh::collapseEdge(int edgeidx)
{
    assert(canCollapseEdge(edgeidx));
    OMMesh::EdgeHandle eh = mesh_.edge_handle(edgeidx);
    OMMesh::HalfedgeHandle heh = mesh_.halfedge_handle(eh, 0);
    OMMesh::VertexHandle tov = mesh_.to_vertex_handle(heh);
    OMMesh::VertexHandle frv = mesh_.from_vertex_handle(heh);
    OMMesh::Point midpt = (mesh_.point(tov) + mesh_.point(frv))*0.5;

    OMMesh newmesh;

    map<int, int> old2new;
    int newidx = 0;
    int mididx = -1;
    for(int i=0; i<(int)mesh_.n_vertices(); i++)
    {
        OMMesh::Point newpt;
        if(i == frv.idx() || i == tov.idx())
        {
            if(mididx == -1)
            {
                newpt = midpt;
                mididx = newidx;
            }
            else
            {
                old2new[i] = mididx;
                continue;
            }
        }
        else
            newpt = mesh_.point(mesh_.vertex_handle(i));
        newmesh.add_vertex(newpt);
        old2new[i] = newidx++;
    }
    assert(newidx == (int)mesh_.n_vertices()-1);

    for(int i=0; i<(int)mesh_.n_faces(); i++)
    {
        OMMesh::FaceHandle fh = mesh_.face_handle(i);
        if(fh == mesh_.face_handle(heh)
                || fh == mesh_.face_handle(mesh_.opposite_halfedge_handle(heh)))
            continue;

        vector<OMMesh::VertexHandle> newface;
        for(OMMesh::FaceVertexIter fvi = mesh_.fv_iter(fh); fvi; ++fvi)
            newface.push_back(newmesh.vertex_handle(old2new[fvi.handle().idx()]));
        newmesh.add_face(newface);
    }

    // fix edges
    for(int i=0; i<(int)boundaries_.size(); i++)
    {
        vector<int> newedges;
        for(int j=0; j<(int)boundaries_[i].edges.size(); j++)
        {
            int oldedgeidx = boundaries_[i].edges[j];
            OMMesh::EdgeHandle oldeh = mesh_.edge_handle(oldedgeidx);
            if(oldeh.idx() == eh.idx())
                continue;
            OMMesh::HalfedgeHandle oldhe = mesh_.halfedge_handle(oldeh,0);
            assert(mesh_.is_boundary(oldeh));
            int v1 = mesh_.to_vertex_handle(oldhe).idx();
            int v2 = mesh_.from_vertex_handle(oldhe).idx();
            int newv1 = old2new[v1];
            int newv2 = old2new[v2];

            bool found = false;

            for(int k=0; k<(int)newmesh.n_halfedges(); k++)
            {
                OMMesh::HalfedgeHandle newhe = newmesh.halfedge_handle(k);
                if(newv1 == newmesh.to_vertex_handle(newhe).idx()
                        && newv2 == newmesh.from_vertex_handle(newhe).idx())
                {
                    assert(!found);
                    found = true;
                    OMMesh::EdgeHandle newbdry = newmesh.edge_handle(newhe);
                    assert(newmesh.is_boundary(newbdry));
                    newedges.push_back(newbdry.idx());
                }

            }
            assert(found);
        }
        boundaries_[i].edges = newedges;
    }
    mesh_ = newmesh;
}

bool DevelopableMesh::collapseShortEdges()
{
    bool atleastone = false;
    bool collapsed = false;
    do
    {
        collapsed = false;

        vector<double> widths;
        computeCreaseWidths(widths);

        vector<double> vertwidths;
        for(int i=0; i<(int)mesh_.n_vertices(); i++)
        {
            OMMesh::VertexHandle vh = mesh_.vertex_handle(i);
            double maxwidth = 0;
            for(OMMesh::VertexEdgeIter vei = mesh_.ve_iter(vh); vei; ++vei)
            {
                maxwidth = std::max(maxwidth, widths[vei.handle().idx()]);
            }
            vertwidths.push_back(maxwidth);
        }

        for(int i=0; i<(int)mesh_.n_edges(); i++)
        {
            OMMesh::EdgeHandle eh = mesh_.edge_handle(i);
            OMMesh::HalfedgeHandle heh = mesh_.halfedge_handle(eh, 0);
            double maxlen = vertwidths[mesh_.from_vertex_handle(heh).idx()];
            maxlen = std::max(maxlen, vertwidths[mesh_.to_vertex_handle(heh).idx()]);

            double len = mesh_.calc_edge_length(eh);
            if(len < 1e-1)//maxlen)
            {
                if(canCollapseEdge(i))
                {
                    collapseEdge(i);
                    collapsed=true;
                    atleastone = true;
                    break;
                }
            }
        }
    } while(collapsed);

    return atleastone;
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
}

void DevelopableMesh::computeCreaseWidths(std::vector<double> &widths)
{
    widths.clear();

    for(int i=0; i<(int)mesh_.n_edges(); i++)
    {
        if(mesh_.is_boundary(mesh_.edge_handle(i)))
        {
            widths.push_back(0);
        }
        else
        {
            double theta = turningAngle(i);
            double L = mesh_.calc_edge_length(mesh_.edge_handle(i));
            widths.push_back(theta*pow(L, 2.0/3.0));
        }
    }
}

double DevelopableMesh::turningAngle(int edge)
{
    OMMesh::EdgeHandle eh = mesh_.edge_handle(edge);
    OMMesh::HalfedgeHandle heh = mesh_.halfedge_handle(eh,0);
    if(mesh_.is_boundary(eh))
    {
        cout << "boundary" << endl;
        return 0;
    }
    OMMesh::VertexHandle v0h = mesh_.from_vertex_handle(heh);
    OMMesh::VertexHandle v00h = mesh_.to_vertex_handle(heh);
    OMMesh::VertexHandle v1h = mesh_.from_vertex_handle(mesh_.prev_halfedge_handle(heh));
    OMMesh::VertexHandle v2h = mesh_.to_vertex_handle(mesh_.next_halfedge_handle(mesh_.opposite_halfedge_handle(heh)));

    Vector3d e0 = point2Vector(mesh_.point(v00h) - mesh_.point(v0h));
    Vector3d e1 = point2Vector(mesh_.point(v1h) - mesh_.point(v0h));
    Vector3d e2 = point2Vector(mesh_.point(v2h) - mesh_.point(v0h));

    Vector3d n0 = e0.cross(e1);
    n0 /= n0.norm();

    Vector3d n1 = e2.cross(e0);
    n1 /= n1.norm();
    return acos(n0.dot(n1));
}

IpoptSolver::IpoptSolver(int n, int m, int nnz_j, int nnz_h, VectorXd &initq, DevelopableMesh &mesh) : n_(n), m_(m), nnz_j_(nnz_j), nnz_h_(nnz_h), initq_(initq), mesh_(mesh)
{

}

bool IpoptSolver::get_nlp_info(Ipopt::Index &n, Ipopt::Index &m, Ipopt::Index &nnz_jac_g, Ipopt::Index &nnz_h_lag, IndexStyleEnum &index_style)
{
    n = n_;
    m = m_;
    nnz_jac_g = nnz_j_;
    nnz_h_lag = nnz_h_;
    index_style = C_STYLE;
    return true;
}

bool IpoptSolver::get_bounds_info(Ipopt::Index n, Ipopt::Number *x_l, Ipopt::Number *x_u, Ipopt::Index m, Ipopt::Number *g_l, Ipopt::Number *g_u)
{
    for(int i=0; i<n; i++)
    {
        x_l[i] = -std::numeric_limits<double>::infinity();
        x_u[i] = std::numeric_limits<double>::infinity();
    }
    for(int i=0; i<m; i++)
    {
        g_l[i] = 0;
        g_u[i] = 0;
    }
    return true;
}

bool IpoptSolver::get_starting_point(Ipopt::Index n, bool init_x, Ipopt::Number *x, bool init_z, Ipopt::Number *z_L, Ipopt::Number *z_U, Ipopt::Index m, bool init_lambda, Ipopt::Number *lambda)
{
    if(init_x)
    {
        for(int i=0; i<n; i++)
            x[i] = initq_[i];
    }
    assert(!init_z);
    assert(!init_lambda);
    return true;
}

bool IpoptSolver::eval_f(Ipopt::Index n, const Ipopt::Number *x, bool new_x, Ipopt::Number &obj_value)
{
    VectorXd q(n);
    for(int i=0; i<n; i++)
        q[i] = x[i];
    double f;
    VectorXd dummy;
    vector<T> dummyT;
    mesh_.buildObjective(q, f, dummy, dummyT);
    obj_value = f;
    return true;
}

bool IpoptSolver::eval_grad_f(Ipopt::Index n, const Ipopt::Number *x, bool new_x, Ipopt::Number *grad_f)
{
    VectorXd q(n);
    for(int i=0; i<n; i++)
        q[i] = x[i];
    double f;
    VectorXd Df;
    vector<T> dummyT;
    mesh_.buildObjective(q, f, Df, dummyT);
    for(int i=0; i<n; i++)
        grad_f[i] = Df[i];
    return true;
}

bool IpoptSolver::eval_g(Ipopt::Index n, const Ipopt::Number *x, bool new_x, Ipopt::Index m, Ipopt::Number *g)
{
    VectorXd q(n);
    for(int i=0; i<n; i++)
        q[i] = x[i];
    VectorXd myg;
    vector<T> dummyT;
    vector<vector<T> > dummyH;
    mesh_.buildConstraints(q, myg, dummyT, dummyH);
    for(int i=0; i<m; i++)
        g[i] = myg[i];
    return true;
}

bool IpoptSolver::eval_jac_g(Ipopt::Index n, const Ipopt::Number *x, bool new_x, Ipopt::Index m, Ipopt::Index nele_jac, Ipopt::Index *iRow, Ipopt::Index *jCol, Ipopt::Number *values)
{
    if(iRow != NULL && jCol != NULL)
    {
        VectorXd g;
        vector<T> Dg;
        vector<vector<T> > dummyH;
        mesh_.buildConstraints(initq_, g, Dg, dummyH);
        assert(Dg.size() == nele_jac);
        for(int i=0; i<Dg.size(); i++)
        {
            iRow[i] = Dg[i].row();
            jCol[i] = Dg[i].col();
        }
    }
    else
    {
        VectorXd q(n);
        for(int i=0; i<n; i++)
            q[i] = x[i];
        VectorXd g;
        vector<T> Dg;
        vector<vector<T> > dummyH;
        mesh_.buildConstraints(q, g, Dg, dummyH);
        for(int i=0; i<Dg.size(); i++)
            values[i] = Dg[i].value();
    }
    return true;
}

bool IpoptSolver::eval_h(Ipopt::Index n, const Ipopt::Number *x, bool new_x, Ipopt::Number obj_factor, Ipopt::Index m, const Ipopt::Number *lambda, bool new_lambda, Ipopt::Index nele_hess, Ipopt::Index *iRow, Ipopt::Index *jCol, Ipopt::Number *values)
{
    if(values == NULL)
    {
        double f;
        VectorXd Df;
        vector<T> Hf;
        mesh_.buildObjective(initq_, f, Df, Hf);

        VectorXd g;
        vector<T> dummyT;
        vector<vector<T> > Hg;

        mesh_.buildConstraints(initq_, g, dummyT, Hg);
        int idx=0;
        for(int i=0; i<Hf.size(); i++)
        {
            iRow[idx] = Hf[i].row();
            jCol[idx] = Hf[i].col();
            idx++;
        }
        for(int i=0; i<Hg.size();i++)
        {
            for(int j=0; j<Hg[i].size(); j++)
            {
                iRow[idx] = Hg[i][j].row();
                jCol[idx] = Hg[i][j].col();
                idx++;
            }
        }
        assert(idx == nele_hess);
    }
    else
    {
        VectorXd q(n);
        for(int i=0; i<n; i++)
        {
            q[i] = x[i];
        }

        double f;
        VectorXd Df;
        vector<T> Hf;
        mesh_.buildObjective(q, f, Df, Hf);

        VectorXd g;
        vector<T> dummyT;
        vector<vector<T> > Hg;

        mesh_.buildConstraints(q, g, dummyT, Hg);
        int idx=0;
        for(int i=0; i<Hf.size(); i++)
        {
            values[idx] = obj_factor*Hf[i].value();
            idx++;
        }
        for(int i=0; i<Hg.size();i++)
        {
            for(int j=0; j<Hg[i].size(); j++)
            {
                values[idx] = lambda[i]*Hg[i][j].value();
                idx++;
            }
        }
        assert(idx == nele_hess);
    }
    return true;
}

void IpoptSolver::finalize_solution(Ipopt::SolverReturn status, Ipopt::Index n, const Ipopt::Number *x, const Ipopt::Number *z_L, const Ipopt::Number *z_U, Ipopt::Index m, const Ipopt::Number *g, const Ipopt::Number *lambda, Ipopt::Number obj_value, const Ipopt::IpoptData *ip_data, Ipopt::IpoptCalculatedQuantities *ip_cq)
{
    cout << "Status " << status << " final objective " << obj_value << endl;
    for(int i=0; i<n; i++)
        initq_[i] = x[i];
}
