#pragma once

#include "OpenMeshCraft/Geometry/Primitives/GenericPoint2T.h"
#include "OpenMeshCraft/Geometry/Primitives/GenericPoint3T.h"

#include "OpenMeshCraft/NumberTypes/NumberUtils.h"

#include "OpenMeshCraft/Utils/Macros.h"

#include <tuple>

namespace OMC {

#define TEMPLATE_DECL template <typename FT, typename IT, typename ET>

/// @brief Dot product 2D
TEMPLATE_DECL
class DotProductSign2D_Indirect
{
public:
	using VecT   = Vec2T<FT>;
	using PointT = GenericPoint2T<IT, ET>;

public:
	/// @brief Dot product between 2D vectors.
	/// @return sign of (p-q).dot(r-q)
	Sign operator()(const PointT &p, const PointT &r, const PointT &q);

	/// @brief Dot product between 2D vectors.
	/// @return sign of (p-q).dot(r-s)
	Sign operator()(const PointT &p, const PointT &r, const PointT &q,
	                const PointT &s);
};

/// @brief Dot product 3D
TEMPLATE_DECL
class DotProductSign3D_Indirect
{
public:
	using VecT   = Vec3T<FT>;
	using PointT = GenericPoint3T<IT, ET>;

public:
	/// @brief Dot product between 3D vectors.
	/// @return sign of (p-q).dot(r-q)
	Sign operator()(const PointT &p, const PointT &r, const PointT &q);

	/// @brief Dot product between 3D vectors.
	/// @return sign of (p-q).dot(r-s)
	Sign operator()(const PointT &p, const PointT &r, const PointT &q,
	                const PointT &s);
};

/// @brief Dot product of 3D points on 2D planes
TEMPLATE_DECL
class DotProductSignOn2D_Indirect
{
public:
	using VecT   = Vec3T<FT>;
	using PointT = GenericPoint3T<IT, ET>;

public:
	/// @brief Dot product between 3D vectors on xy.
	/// @return sign of (p-q).dot(r-q)
	Sign on_xy(const PointT &p, const PointT &r, const PointT &q);

	/// @brief Dot product between 3D vectors on xy.
	/// @return sign of (p-q).dot(r-s)
	Sign on_xy(const PointT &p, const PointT &r, const PointT &q,
	           const PointT &s);

	/// @brief Dot product between 3D vectors on yz.
	/// @return sign of (p-q).dot(r-q)
	Sign on_yz(const PointT &p, const PointT &r, const PointT &q);

	/// @brief Dot product between 3D vectors on yz.
	/// @return sign of (p-q).dot(r-s)
	Sign on_yz(const PointT &p, const PointT &r, const PointT &q,
	           const PointT &s);

	/// @brief Dot product between 3D vectors on zx.
	/// @return sign of (p-q).dot(r-q)
	Sign on_zx(const PointT &p, const PointT &r, const PointT &q);

	/// @brief Dot product between 3D vectors on zx.
	/// @return sign of (p-q).dot(r-s)
	Sign on_zx(const PointT &p, const PointT &r, const PointT &q,
	           const PointT &s);
};

/// @brief  Orinet2d
TEMPLATE_DECL
class Orient2D_Indirect
{
public:
	using VecT   = Vec2T<FT>;
	using PointT = GenericPoint2T<IT, ET>;

public:
	/// @brief test orientation of query with respect to 2D line pq.
	/// left: POSITIVE, on : ZERO, right: NEGATIVE.       (q-p) cross (query-p)
	Sign operator()(const PointT &p, const PointT &q, const PointT &query);

	/// @brief test orientation of query with respect to 2D line pq.
	/// left: POSITIVE, on : ZERO, right: NEGATIVE.       (q-p) cross (query-p)
	Sign operator()(const FT *p, const FT *q, const FT *query);
};

/// @brief  Square distance in 2D.
TEMPLATE_DECL
class SquareDistance2D_Indirect
{
public:
	using VecT   = Vec2T<FT>;
	using PointT = GenericPoint2T<IT, ET>;

public:
	/// @brief Calculate the square distance between \p p and \p q ,
	/// then compare with \p sqr_dis.
	/// @return Sign. NEGATIVE: calculated square distance is less than \p
	/// sqr_dis. ZERO: calculated square distance is equal to \p sqr_dis.
	/// POSITIVE: calculated square distance is larger than \p sqr_dis.
	Sign operator()(const PointT &p, const PointT &q, FT sqr_dis);
};

/// @brief  Square distance in 3D.
TEMPLATE_DECL
class SquareDistance3D_Indirect
{
public:
	using VecT   = Vec3T<FT>;
	using PointT = GenericPoint3T<IT, ET>;

public:
	/// @brief Calculate the square distance between \p p and \p q ,
	/// then compare with \p sqr_dis.
	/// @return Sign. NEGATIVE: calculated square distance is less than \p
	/// sqr_dis. ZERO: calculated square distance is equal to \p sqr_dis.
	/// POSITIVE: calculated square distance is larger than \p sqr_dis.
	Sign operator()(const PointT &p, const PointT &q, FT sqr_dis);
};

/// @brief  InCircle
TEMPLATE_DECL
class InCircle_Indirect
{
public:
	using VecT   = Vec2T<FT>;
	using PointT = GenericPoint2T<IT, ET>;

public:
	/**
	 * @brief In 2D, test a point \p pd is inside the circumcircle of three
	 * points \p pa, \p pb, and \p pc. We assume that pa, pb and pc is given in
	 * couter-clock-wise order.
	 * @return POSITIVE: inside, ZERO: exactly on the circumcircle, NEGATIVE:
	 * outside.
	 */
	Sign operator()(const PointT &pa, const PointT &pb, const PointT &pc,
	                const PointT &pd);

	/**
	 * @brief In 2D, test a point \p pd is inside the circumcircle of three
	 * points \p pa, \p pb, and \p pc. We assume that pa, pb and pc is given in
	 * couter-clock-wise order.
	 * @return POSITIVE: inside, ZERO: exactly on the circumcircle, NEGATIVE:
	 * outside.
	 */
	Sign operator()(const FT *pa, const FT *pb, const FT *pc, const FT *pd);
};

TEMPLATE_DECL
class CollinearPoints2D_Indirect
{
public:
	using VecT     = Vec2T<FT>;
	using PointT   = GenericPoint2T<IT, ET>;
	using Orient2D = Orient2D_Indirect<FT, IT, ET>;

public:
	bool operator()(const PointT &p, const PointT &q, const PointT &r)
	{
		return Orient2D()(p, q, r) == Sign::ZERO;
	}
};

/// @brief ad = b-a, bd = c-a, cd = d-a
/// det = dot(cross(ad,bd), cd), ZERO -> coplanar, POSITIVE -> positive volume,
/// NEGATIVE -> negative volume.
/// @details In geometric sense, a, b, c come from a triangle, and d is the
/// query point.
TEMPLATE_DECL
class Orient3D_Indirect
{
public:
	using VecT   = Vec3T<FT>;
	using PointT = GenericPoint3T<IT, ET>;

public:
	/// @brief ad = b-a, bd = c-a, cd = d-a
	/// det = dot(cross(ad,bd), cd), ZERO -> coplanar, POSITIVE -> positive
	/// volume, NEGATIVE -> negtive volume.
	Sign operator()(const PointT &a, const PointT &b, const PointT &c,
	                const PointT &d);

	/// @brief ad = b-a, bd = c-a, cd = d-a
	/// det = dot(cross(ad,bd), cd), ZERO -> coplanar, POSITIVE -> positive
	/// volume, NEGATIVE -> negtive volume.
	Sign operator()(const FT *a, const FT *b, const FT *c, const FT *d);

public: // for cache
	/// @brief Calculate cached data for three points pa, pb and pc (which
	/// generally come from three points of a triangle)
	static void get_minors(const FT *pa, const FT *pb, const FT *pc, FT *minor,
	                       FT *perm);

	/// @brief Calculate orient3d with cached data. Be careful to put points
	/// in right order!
	/// Cached data are calculated for pa, pb and pc. The query point is pd.
	static Sign with_cached_minors(const FT *pa, const FT *pb, const FT *pc,
	                               const FT *pd, const FT *minor, const FT *perm);

public: // separated filter stage
	Sign filter(const FT *a, const FT *b, const FT *c, const FT *d);
};

/**
 * @brief Check 3D points' orientation on 2D planes (xy, yz or zx plane).
 */
TEMPLATE_DECL
class OrientOn2D_Indirect
{
public:
	using VecT   = Vec3T<FT>;
	using PointT = GenericPoint3T<IT, ET>;

	using OrientOn2D = OrientOn2D_Indirect<FT, IT, ET>;

public:
	// Faster. It assumes that points are coplanar and the dominant normal
	// component is n_max (see maxComponentInTriangleNormal()).
	Sign operator()(const PointT &a, const PointT &b, const PointT &c, int n_max);

	Sign on_xy(const PointT &a, const PointT &b, const PointT &c);
	Sign on_yz(const PointT &a, const PointT &b, const PointT &c);
	Sign on_zx(const PointT &a, const PointT &b, const PointT &c);

	Sign operator()(const FT *a, const FT *b, const FT *c, int n_max);

	Sign on_xy(const FT *a, const FT *b, const FT *c);
	Sign on_yz(const FT *a, const FT *b, const FT *c);
	Sign on_zx(const FT *a, const FT *b, const FT *c);
};

TEMPLATE_DECL
class LessThan2D_Indirect
{
public:
	using VecT   = Vec2T<FT>;
	using PointT = GenericPoint2T<IT, ET>;

public:
	// TODO implement this

	// lessThanOnX (resp. Y)
	// lessThanOnX(a,b) =
	// NEGATIVE - if a.X < b.X
	// ZERO  - if a.X == b.X
	// POSITIVE  - if a.X > b.X
	Sign on_x(const PointT &a, const PointT &b);
	Sign on_y(const PointT &a, const PointT &b);

	// lessThan
	// Input points can be any combination of 2D points
	// lessThan(a,b) =
	// NEGATIVE - if a < b
	// ZERO  - if a == b
	// POSITIVE  - if a > b
	// in lexicographical order
	Sign operator()(const PointT &a, const PointT &b);

	// TRUE if the two points are coincident
	bool coincident(const PointT &a, const PointT &b)
	{
		return operator()(a, b) == Sign::ZERO;
	}
};

TEMPLATE_DECL
class LessThan3D_Indirect
{
public:
	using VecT   = Vec3T<FT>;
	using PointT = GenericPoint3T<IT, ET>;

public:
	// lessThanOnX (resp. Y, Z)
	// Input points can be any combination of 3D points
	// lessThanOnX(a,b) =
	// NEGATIVE - if a.X < b.X
	// ZERO  - if a.X == b.X
	// POSITIVE  - if a.X > b.X
	Sign on_x(const PointT &a, const PointT &b);
	Sign on_y(const PointT &a, const PointT &b);
	Sign on_z(const PointT &a, const PointT &b);

	Sign on(const PointT &a, const PointT &b, size_t dim);

	std::array<Sign, 3> on_all(const PointT &a, const PointT &b);

	Sign on_x(const PointT &a, const FT *b);
	Sign on_y(const PointT &a, const FT *b);
	Sign on_z(const PointT &a, const FT *b);

	Sign on(const PointT &a, const FT *b, size_t dim);

	std::array<Sign, 3> on_all(const PointT &a, const FT *b);

	Sign on_x(const PointT &a, const FT b);
	Sign on_y(const PointT &a, const FT b);
	Sign on_z(const PointT &a, const FT b);

	Sign on(const PointT &a, const FT b, size_t dim);

	// lessThan
	// Input points can be any combination of 3D points
	// lessThan(a,b) =
	// NEGATIVE - if a < b
	// ZERO  - if a == b
	// POSITIVE  - if a > b
	// in lexicographical order
	Sign operator()(const PointT &a, const PointT &b);

	// TRUE if the two points are coincident
	bool coincident(const PointT &a, const PointT &b)
	{
		return operator()(a, b) == Sign::ZERO;
	}
};

TEMPLATE_DECL
class CollinearPoints3D_Indirect
{
public:
	using VecT   = Vec3T<FT>;
	using PointT = GenericPoint3T<IT, ET>;

	using OrientOn2D = OrientOn2D_Indirect<FT, IT, ET>;

public:
	/// @brief check if p, q, r are collinear.
	/// @param p
	/// @param q
	/// @param r
	/// @return true if they are collinear, otherwise false.
	bool operator()(const PointT &p, const PointT &q, const PointT &r)
	{
		return !misaligned(p, q, r);
	}

	bool operator()(const FT *p, const FT *q, const FT *r)
	{
		return !misaligned(p, q, r);
	}

	// TRUE if A-B-C are not collinear
	bool misaligned(const PointT &A, const PointT &B, const PointT &C);

	// Faster. It assumes that points are coplanar and the dominant normal
	// component is n_max (see MaxComponentInTriangleNormal).
	bool misaligned(const PointT &A, const PointT &B, const PointT &C, int n_max);

	bool misaligned(const FT *A, const FT *B, const FT *C);
};

TEMPLATE_DECL
class CollinearSort3D_Indirect
{
public:
	using VecT   = Vec3T<FT>;
	using PointT = GenericPoint3T<IT, ET>;

	using LessThan3D = LessThan3D_Indirect<FT, IT, ET>;

public:
	/**
	 * @brief Sort three collinear points along the line.
	 * @return A tuple contains three references that are ordered.
	 */
	std::tuple<const PointT &, const PointT &, const PointT &>
	operator()(const PointT &p, const PointT &q, const PointT &r);
};

TEMPLATE_DECL
class MaxComponentInTriangleNormal
{
public:
	// Let n = (x,y,z) be the normal of the triangle <v1,v2,v3>
	// and let m be the absolute value of its largest component.
	// That is, m = max(|x|, |y|, |z|).
	// maxComponentInTriangleNormal(v1,v2,v3) returns:
	// 0 - if m == |x|
	// 1 - if m == |y|
	// 2 - if m == |z|
	//
	// Warning: this function assumes that the triangle is not exactly
	// degenerate. It may crash otherwise.
	int operator()(FT v1x, FT v1y, FT v1z, FT v2x, FT v2y, FT v2z, FT v3x, FT v3y,
	               FT v3z);

	int operator()(const FT *v1, const FT *v2, const FT *v3)
	{
		return operator()(v1[0], v1[1], v1[2], v2[0], v2[1], v2[2], v3[0], v3[1],
		                  v3[2]);
	}
};

TEMPLATE_DECL
class LongestAxis_Indirect
{
public:
	using VecT   = Vec3T<FT>;
	using PointT = GenericPoint3T<IT, ET>;

public:
	int operator()(const PointT &a, const PointT &b);
};

TEMPLATE_DECL
class InSphere_Indirect
{
public:
	using VecT   = Vec3T<FT>;
	using PointT = GenericPoint3T<IT, ET>;

public:
	/**
	 * @brief In 3D, test a point \p e is inside the circumsphere of four
	 * points \p a, \p b, \p c, and \p d.
	 * @return POSITIVE: inside, ZERO: exactly on the circumcircle, NEGATIVE:
	 * outside.
	 */
	Sign operator()(const PointT &a, const PointT &b, const PointT &c,
	                const PointT &d, const PointT &e);

	Sign operator()(const FT *a, const FT *b, const FT *c, const FT *d,
	                const FT *e);
};

} // namespace OMC

#ifdef OMC_HAS_IMPL
	#include "IndirectPredicates.inl"
#endif