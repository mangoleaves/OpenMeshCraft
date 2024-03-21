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

#ifdef OMC_INTER_PROFILE

enum class IntersectionNames : size_t
{
	T3T3 = 0,
	T3S3,
	CNT
};

struct IntersectionProfile
{
	static const uint32_t BRANCH_CNT = 32;

	static uint32_t total_count[(size_t)IntersectionNames::CNT];
	static uint32_t reach_count[(size_t)IntersectionNames::CNT][BRANCH_CNT];
	static uint32_t reach_line[(size_t)IntersectionNames::CNT][BRANCH_CNT];

	static void initialize();
	static void inc_total(IntersectionNames name);
	static void inc_reach(IntersectionNames name, uint32_t branch_flag,
	                      uint32_t branch_line);

	static void print();
};

	#define OMC_INTER_PROFILE_INIT OMC::IntersectionProfile::initialize()
	#define OMC_INTER_PROFILE_PRINT OMC::IntersectionProfile::print()

	#define OMC_INTER_PROFILE_INC_TOTAL(inter) \
		OMC::IntersectionProfile::inc_total(inter)

	#define OMC_INTER_PROFILE_INC_REACH(inter, branch_flag) \
		OMC::IntersectionProfile::inc_reach(inter, branch_flag, __LINE__)

#else

	#define OMC_INTER_PROFILE_INIT
	#define OMC_INTER_PROFILE_PRINT

	#define OMC_INTER_PROFILE_INC_TOTAL(inter)
	#define OMC_INTER_PROFILE_INC_REACH(inter, branch_flag)

#endif

} // namespace OMC