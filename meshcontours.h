#ifndef __meshcontours_h__
#define __meshcontours_h__

#include <vector>

#include <mesh.h>

typedef std::vector< std::pair< OMMesh::Point, OMMesh::Point > > contour_t;

class MeshContours
{
public:
    const void renderContours();
    
    void clearPrecomputedContours();
    void precomputeContoursAtValuesForMeshUsingVertexValues( const std::vector< double >& isovalues, const OMMesh& mesh, const std::vector< double >& vertex_values );
    
private:
    std::vector< contour_t > m_contours;
};

#endif /* __meshcontours_h__ */
