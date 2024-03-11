#pragma once

#include "OpenMeshCraft/NumberTypes/NumberUtils.h"
#include "OpenMeshCraft/Utils/Exception.h"

namespace OMC {

enum class PointInType
{
	STRICTLY_INSIDE  = -1, // strictly inside  the input simplex
	ON_BOUNDARY      = 0,  // point is on the boundary (vert, edge, face...)
	STRICTLY_OUTSIDE = 1   // strictly outside the input simplex
};

// location of intersection points for point_in_{segment | triangle | tet}
// predicates. Elements' orders are compliant with the tables in:
enum class PointInSimplexType
{
	STRICTLY_OUTSIDE = 0,  // strictly outside the input simplex
	STRICTLY_INSIDE  = 1,  // strictly inside  the input simplex
	ON_VERT0         = 2,  // used for segs, tris and tets
	ON_VERT1         = 3,  // used for segs, tris and tets
	ON_VERT2         = 4,  // used for tris and tets
	ON_VERT3         = 5,  // used for tets
	ON_EDGE0         = 6,  // used for tris and tets
	ON_EDGE1         = 7,  // used for tris and tets
	ON_EDGE2         = 8,  // used for tris and tets
	ON_EDGE3         = 9,  // used for tets
	ON_EDGE4         = 10, // used for tets
	ON_EDGE5         = 11, // used for tets
	ON_FACE0         = 12, // used for tets
	ON_FACE1         = 13, // used for tets
	ON_FACE2         = 14, // used for tets
	ON_FACE3         = 15  // used for tets
};

// intersection types
enum class SimplexIntersectionType
{
	DO_NOT_INTERSECT   = 0, // simplices do not intersect
	SIMPLICIAL_COMPLEX = 1, // simplices form a valid simplicial complex (i.e.
	                        // they are coincident or share a sub-simplex)
	INTERSECT          = 2, // simplices intersect in a non conforming way
	OVERLAP = 3 // for corner cases: simplices intersect and partially overlap
};            // (e.g. collinear segments or coplanar triangles)

template <typename NT>
inline bool vec_equals_2d(const NT *v0, const NT *v1)
{
	return ((v0[0] == v1[0]) && (v0[1] == v1[1]));
}

template <typename NT>
inline bool vec_equals_3d(const NT *v0, const NT *v1)
{
	return ((v0[0] == v1[0]) && (v0[1] == v1[1]) && (v0[2] == v1[2]));
}

} // namespace OMC