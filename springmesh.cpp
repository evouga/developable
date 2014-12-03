#include "springmesh.h"

void SpringMesh::setupSprings(){
//    std::vector<std::vector<int> >::iterator i = this->springs.begin();

//    for(; i != this->springs.end(); i++)
//        i->clear();
//    this->springs.clear();

//    unsigned int vertex_no = 0;
//    OMMesh::FaceIter fi = mesh_.faces_begin();
//    for(; fi != mesh_.faces_end(); fi++){
//        std::vector<int> current_indices;
//        int current = 0;
//        OMMesh::FaceVertexIter fvi = mesh_.fv_iter(fi);
//        for(;fvi.is_valid(); fvi++){
//            current_indices.push_back(3*(*fi)+current);
//            vertex_no++;
//        }
//    }
}


double SpringMesh::elasticEnergy(const Eigen::VectorXd& q){

    double energy = 0.0;
    OMMesh::EdgeIter ei = mesh_.edges_begin();
    for(; ei != mesh_.edges_end(); ei++){

        if(mesh_.is_boundary(*ei)){
            continue;
        }
        OMMesh::HalfedgeHandle h1= mesh_.halfedge_handle(*ei,0);
        OMMesh::HalfedgeHandle h2= mesh_.halfedge_handle(*ei,1);

        OMMesh::FaceHandle f1 = mesh_.face_handle(h1);
        OMMesh::FaceHandle f2 = mesh_.face_handle(h2);

        OMMesh::VertexHandle v1 = mesh_.from_vertex_handle(h1);
        OMMesh::VertexHandle v2 = mesh_.from_vertex_handle(h2);

        unsigned int local_vertex_index1 = 0;
        unsigned int local_vertex_index2 = 0;

        for(OMMesh::FaceVertexIter fvi = mesh_.fv_iter(f1); *fvi==v1; fvi++){
            local_vertex_index1++;
        }
        for(OMMesh::FaceVertexIter fvi = mesh_.fv_iter(f2); *fvi==v2; fvi++){
            local_vertex_index2++;
        }
        assert((local_vertex_index1 <= 2) && (local_vertex_index2 <= 2));

        Eigen::Vector3d q11 = q.segment<3>(3*f1.idx()+local_vertex_index1);
        Eigen::Vector3d q12 = q.segment<3>(3*f2.idx()+local_vertex_index2);

        local_vertex_index1 = (local_vertex_index1+1)%3;
        local_vertex_index2 = (local_vertex_index2+2)%3;
        Eigen::Vector3d q21 = q.segment<3>(3*f1.idx()+local_vertex_index1);
        Eigen::Vector3d q22 = q.segment<3>(3*f2.idx()+local_vertex_index2);


        energy += (q11-q12).squaredNorm();
        energy += (q21-q22).squaredNorm();
    }
    energy *= spring_constant;
    return energy;
}
