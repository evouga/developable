#include "meshcontours.h"

#include <cassert>
#include <algorithm>

#include <GL/gl.h>

namespace
{
inline void test_and_push( std::vector< OMMesh::Point >& pts, const double isoval, std::pair< double, int > val_and_index0, std::pair< double, int > val_and_index1, const OMMesh& mesh )
{
    if( val_and_index0.first > val_and_index1.first )
    {
        std::swap( val_and_index0, val_and_index1 );
    }
    
    if( isoval < val_and_index1.first && !( isoval < val_and_index0.first ) )
    {
        const double t = ( isoval - val_and_index0.first ) / ( val_and_index1.first - val_and_index0.first );
        const OMMesh::Point p0 = mesh.point( OMMesh::VertexHandle( val_and_index0.second ) );
        const OMMesh::Point p1 = mesh.point( OMMesh::VertexHandle( val_and_index1.second ) );
        pts.push_back( p0 + ( p1 - p0 )*t );
    }
}

std::vector< contour_t >
ContoursAtValuesForMeshUsingVertexValues(
    const std::vector< double >& isovalues,
    const OMMesh& mesh,
    const std::vector< double >& vertex_values
    )
{
    /// 1 For each face, gather vertex indices.
    /// 2 For each iso-value in 'isovalues',
    ///   if the iso-value lies between the 'vertex_values' of two
    ///   of the face's vertices, remember that point.
    /// 3 There should always be zero or two of such points;
    ///   output a line segment when there are two.
    
    assert( mesh.n_vertices() == vertex_values.size() );
    
    std::vector< contour_t > result;
    result.resize( isovalues.size() );
    
    /// 1
    std::vector< std::pair< double, int > > vals_and_indices;
    std::vector< OMMesh::Point > pts;
    
    for( OMMesh::ConstFaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it )
    {
        vals_and_indices.clear();
        vals_and_indices.reserve(3);
        
        for( OMMesh::ConstFaceVertexIter fv_it = mesh.cfv_iter( f_it.handle() ); fv_it; ++fv_it )
        {
             vals_and_indices.push_back( std::make_pair(
                vertex_values.at( fv_it.handle().idx() ),
                fv_it.handle().idx()
                ) );
        }
        assert( vals_and_indices.size() == 3 );
        
        /// 2
        for( unsigned int i = 0; i < isovalues.size(); ++i )
        {
            pts.clear();

            const double val = isovalues.at( i );
            
            for( int c = 0; c < 3; ++c )
            {
                test_and_push( pts, val, vals_and_indices[c], vals_and_indices[(c+1)%3], mesh );
            }
            
            /// 3
            // TODO: There is a bug here where contours are not closed.
            //       When I check if pts.size() > 1, I see fewer closed contours
            //       but some spurious lines.
            if( pts.size() == 2 )
            {
                assert( pts.size() == 0 || pts.size() == 2 );
                result[i].push_back( std::make_pair( pts.front(), pts.back() ) );
            }
        }
    }
    
    return result;
}
}

void MeshContours::clearPrecomputedContours()
{
    m_contours.clear();
}

void MeshContours::precomputeContoursAtValuesForMeshUsingVertexValues(
    const std::vector< double >& isovalues,
    const OMMesh& mesh,
    const std::vector< double >& vertex_values
    )
{
    m_contours = ContoursAtValuesForMeshUsingVertexValues( isovalues, mesh, vertex_values );
}

void MeshContours::renderContours() const
{
    glDisable( GL_LIGHTING );
    
    glLineWidth( 1.0 );
    glBegin( GL_LINES );
    for( unsigned int i = 0; i < m_contours.size(); i++ )
    for( unsigned int j = 0; j < m_contours[i].size(); j++ )
    {
        glColor3f( 0., 0., 0. );
        
        glVertex3fv( m_contours[i][j].first.data() );
        glVertex3fv( m_contours[i][j].second.data() );
    }
    glEnd();
}
