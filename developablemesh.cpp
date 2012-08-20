#include "developablemesh.h"

using namespace std;

const double PI = 3.1415926535898;

DevelopableMesh::DevelopableMesh()
{

}

bool DevelopableMesh::loadMesh(const string &filename)
{
    if(!Mesh::loadMesh(filename))
        return false;
    identifyBoundaries();
    return true;
}

void DevelopableMesh::buildSchwarzLantern(double r, double h, int n, int m)
{
    mesh_ = OMMesh();

    for(int i=0; i<m; i++)
    {
        double z = h*i;
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
            while(!visited[eh.idx()])
            {
                visited[eh.idx()] = true;
                bc.edges.push_back(eh.idx());
                bc.arclength += mesh_.calc_edge_length(eh);

                heh = mesh_.next_halfedge_handle(heh);
                eh = mesh_.edge_handle(heh);
            }

            boundaries_.push_back(bc);

        }

    }

    delete[] visited;
}
