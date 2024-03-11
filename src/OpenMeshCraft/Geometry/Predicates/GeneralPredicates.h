#pragma once

#include "OpenMeshCraft/NumberTypes/NumberUtils.h"

#include <tuple>

namespace OMC {

/// @brief Dot product 2D
template <typename NT>
class DotProductSign2D_GNR
{
public:
	using VecT   = Vec2T<NT>;
	using PointT = Point2T<NT>;

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
template <typename NT>
class DotProductSign3D_GNR
{
public:
	using VecT   = Vec3T<NT>;
	using PointT = Point3T<NT>;

public:
	/// @brief Dot product between 3D vectors.
	/// @return sign of (p-q).dot(r-q)
	Sign operator()(const PointT &p, const PointT &r, const PointT &q);

	/// @brief Dot product between 3D vectors.
	/// @return sign of (p-q).dot(r-s)
	Sign operator()(const PointT &p, const PointT &r, const PointT &q,
	                const PointT &s);
};

/// @brief Dot product of 3D vectors on 2D
template <typename NT>
class DotProductSignOn2D_GNR
{
public:
	using VecT   = Vec3T<NT>;
	using PointT = Point3T<NT>;

public:
	/// @brief Dot product between 3D vectors on 2D.
	/// @return sign of (p-q).dot(r-q) on xy plane
	Sign on_xy(const PointT &p, const PointT &r, const PointT &q);

	/// @brief Dot product between 3D vectors on 2D.
	/// @return sign of (p-q).dot(r-s) on xy plane
	Sign on_xy(const PointT &p, const PointT &r, const PointT &q,
	           const PointT &s);

	/// @brief Dot product between 3D vectors on 2D.
	/// @return sign of (p-q).dot(r-q) on xy plane
	Sign on_yz(const PointT &p, const PointT &r, const PointT &q);

	/// @brief Dot product between 3D vectors on 2D.
	/// @return sign of (p-q).dot(r-s) on xy plane
	Sign on_yz(const PointT &p, const PointT &r, const PointT &q,
	           const PointT &s);

	/// @brief Dot product between 3D vectors on 2D.
	/// @return sign of (p-q).dot(r-q) on xy plane
	Sign on_zx(const PointT &p, const PointT &r, const PointT &q);

	/// @brief Dot product between 3D vectors on 2D.
	/// @return sign of (p-q).dot(r-s) on xy plane
	Sign on_zx(const PointT &p, const PointT &r, const PointT &q,
	           const PointT &s);
};

/**
 * @brief General orientation check in 2D.
 * Given an oriented line in 2D (the line has a direction, so it has
 * left halfspace and right halfspace), check if a query point is on the
 * left of the line, the right of the line, or exactly the line itself.
 * @tparam NT the number type.
 */
template <typename NT>
class Orient2D_GNR
{
public:
	using VecT   = Vec2T<NT>;
	using PointT = Point2T<NT>;

public:
	/**
	 * @brief \p p and \p q define an oriented line whose direction is from \p p
	 * to \p q . Check the orientation of the query with respect to Line(p,q).
	 * @return 1: in left halfspace. 0: on the plane. -1: in right halfspace.
	 */
	Sign operator()(const PointT &p, const PointT &q, const PointT &query);

	/**
	 * @brief \p p and \p q define an oriented line whose direction is from \p p
	 * to \p q . Check the orientation of the query with respect to Line(p,q).
	 * @return 1: in left halfspace. 0: on the plane. -1: in right halfspace.
	 */
	Sign operator()(const NT &px, const NT &py, const NT &qx, const NT &qy,
	                const NT &query_x, const NT &query_y);

	// OPT: will it benefit from std::move when number type is exact type?
};

/**
 * @brief General orientation check in 3D.
 * Given an oriented plane in 3D (the plane has normal and the normal points to
 * its positive halfspace), check if a query point is in the positive halfspace,
 * in the negative halfspace, or on the plane.
 * @tparam NT the number type.
 */
template <typename NT>
class Orient3D_GNR
{
public:
	using VecT   = Vec3T<NT>;
	using PointT = Point3T<NT>;

public:
	/**
	 * @brief \p p , \p q and \p r define an oriented plane, whose normal is
	 * (q-p).cross(r-p). Check the orientation of the query with respect to
	 * Plane(p,q,r).
	 * @return 1: in positive halfspace. 0: on the plane. -1: in negative
	 * halfspace.
	 */
	Sign operator()(const PointT &p, const PointT &q, const PointT &r,
	                const PointT &query);

	// OPT: will it benefit from std::move when number type is exact type?
	// OPT: provide an interface receving px,py,pz...
};

/**
 * @brief Check 3D points' orientation on 2D planes (xy, yz or zx plane).
 */
template <typename NT>
class OrientOn2D_GNR
{
public:
	using VecT   = Vec3T<NT>;
	using PointT = Point3T<NT>;

	using Orient2D = Orient2D_GNR<NT>;

public:
	Sign on_xy(const PointT &a, const PointT &b, const PointT &c);
	Sign on_yz(const PointT &a, const PointT &b, const PointT &c);
	Sign on_zx(const PointT &a, const PointT &b, const PointT &c);
};

/**
 * @brief Check if three points in 2D are collinear.
 * @tparam NT the number type
 */
template <typename NT>
class CollinearPoints2D_GNR
{
public:
	using Orient2D = Orient2D_GNR<NT>;
	using PointT   = Point2T<NT>;

public:
	/**
	 * @brief Check if \p p , \p q and \p r are collinear.
	 * @return true if they are collinear, otherwise false.
	 */
	bool operator()(const PointT &p, const PointT &q, const PointT &r);
};

/**
 * @brief Check if a 2D point is less than another 2D point.
 * @tparam NT Number type.
 */
template <typename NT>
class LessThan2D_GNR
{
public:
	using PointT = Point2T<NT>;

public:
	/**
	 * @brief Check if all coordinates of a point \p p are less than coordinates
	 * of another point \p q in lexicographic order.
	 * @return Sign. -1: p is less than q. 0: p is equal to q. -1: p is larger
	 * than q.
	 */
	Sign operator()(const PointT &p, const PointT &q);

	/**
	 * @brief Check if \p p.x is less than \p q.x.
	 * @return Sign. -1: p.x < q.x ; 0: p.x==q.x ; -1: p.x > q.x .
	 */
	Sign on_x(const PointT &p, const PointT &q);

	/**
	 * @brief Check if \p p.y is less than \p q.y.
	 * @return Sign. -1: p.y < q.y ; 0: p.y==q.y ; -1: p.y > q.y .
	 */
	Sign on_y(const PointT &p, const PointT &q);

	/**
	 * @brief Check if \p p and \p q is coincident (same on all coordinates).
	 * @return true if they are coincident
	 */
	bool coincident(const PointT &p, const PointT &q)
	{
		return operator()(p, q) == Sign::ZERO;
	}
};

/**
 * @brief Check if a 3D point is less than another 3D point.
 * @tparam NT Number type.
 */
template <typename NT>
class LessThan3D_GNR
{
public:
	using PointT = Point3T<NT>;

public:
	/**
	 * @brief Check if all coordinates of a point \p p are less than coordinates
	 * of another point \p q in lexicographic order.
	 * @return Sign. -1: p is less than q. 0: p is equal to q. -1: p is larger
	 * than q.
	 */
	Sign operator()(const PointT &p, const PointT &q);

	/**
	 * @brief Check if \p p.x is less than \p q.x.
	 * @return Sign. -1: p.x < q.x ; 0: p.x==q.x ; -1: p.x > q.x .
	 */
	Sign on_x(const PointT &p, const PointT &q);

	/**
	 * @brief Check if \p p.y is less than \p q.y.
	 * @return Sign. -1: p.y < q.y ; 0: p.y==q.y ; -1: p.y > q.y .
	 */
	Sign on_y(const PointT &p, const PointT &q);

	/**
	 * @brief Check if \p p.z is less than \p q.z.
	 * @return Sign. -1: p.z < q.z ; 0: p.z==q.z ; -1: p.z > q.z .
	 */
	Sign on_z(const PointT &p, const PointT &q);

	/**
	 * @brief Check if \p p and \p q is coincident (same on all coordinates).
	 * @return true if they are coincident
	 */
	bool coincident(const PointT &p, const PointT &q)
	{
		return operator()(p, q) == Sign::ZERO;
	}
};

/**
 * @brief Check if three points in 3D are collinear.
 * @tparam NT the number type
 */
template <typename NT>
class CollinearPoints3D_GNR
{
public:
	using Orient2D = Orient2D_GNR<NT>;
	using PointT   = Point3T<NT>;

public:
	/**
	 * @brief Check if \p p , \p q and \p r are collinear.
	 * @return true if they are collinear, otherwise false.
	 */
	bool operator()(const PointT &p, const PointT &q, const PointT &r);
};

/**
 * @brief Given three collinear points in 3D, check if they are ordered along
 * the line.
 * @tparam NT the number type
 */
template <typename NT>
class CollinearOrdered3D_GNR
{
public:
	using PointT = Point3T<NT>;

public:
	/**
	 * @brief Given collinear points \p p , \p q and \p r in 3D, check if they are
	 * ordered along the line.
	 * @return true if they are ordered, otherwise false.
	 * @note If you want sort points, call function `sort`.
	 */
	bool operator()(const PointT &p, const PointT &q, const PointT &r);
};

/**
 * @brief Given three collinear points in 3D, sort them along the line.
 * @tparam NT the number type
 */
template <typename NT>
class CollinearSort3D_GNR
{
public:
	using PointT = Point3T<NT>;

public:
	/**
	 * @brief Sort three collinear points along the line.
	 * @return A tuple contains three references that are ordered.
	 */
	std::tuple<const PointT &, const PointT &, const PointT &>
	operator()(const PointT &p, const PointT &q, const PointT &r);
};

/**
 * @brief Check if four points in 3D are coplanar.
 * @tparam NT the number type
 */
template <typename NT>
class CoplanarPoints3D_GNR
{
public:
	using Orient3D = Orient3D_GNR<NT>;
	using PointT   = Point3T<NT>;

public:
	/**
	 * @brief Check if \p p , \p q , \p r and \p r are coplanar.
	 * @return true if they are coplanar, otherwise false.
	 */
	bool operator()(const PointT &p, const PointT &q, const PointT &r,
	                const PointT &s);
};

template <typename NT>
class InCircle_GNR
{
public:
	using PointT = Point2T<NT>;

public:
	/**
	 * @brief In 2D, test a \p query point is inside the circumcircle of three
	 * points \p p, \p q, and \p r. We assume that p, q and r are given in
	 * couter-clock-wise order when being viewed from +z.
	 * @return 1: inside, 0: exactly on the circumcircle, -1: outside.
	 * @note If p, q and r are given in clock-wise order, the result should be
	 * inversed.
	 */
	Sign operator()(const PointT &p, const PointT &q, const PointT &r,
	                const PointT &query);

	/**
	 * @brief In 2D, test a \p query point is inside the circumcircle of three
	 * points \p p, \p q, and \p r. We assume that p, q and r are given in
	 * couter-clock-wise order when being viewed from +z.
	 * @return 1: inside, 0: exactly on the circumcircle, -1: outside.
	 * @note If p, q and r are given in clock-wise order, the result should be
	 * inversed.
	 */
	Sign operator()(const NT &px, const NT &py, const NT &qx, const NT &qy,
	                const NT &rx, const NT &ry, const NT &queryx,
	                const NT &queryy);
};

template <typename NT>
class InSphere_GNR
{
public:
	using PointT = Point3T<NT>;

public:
	/**
	 * @brief In 3D, test a \p query point is inside the circumcircle of four
	 * points \p a, \p b, \p c and \p d. We assume that the tetrahedron formed by
	 * a, b, c and d has a positive volume. In another word, a, b and c form a
	 * triangle (in CCW order), the vector from a to d must point to the same
	 * direction of the normal.
	 * @return 1: inside, 0: exactly on , -1: outside.
	 */
	Sign operator()(const PointT &a, const PointT &b, const PointT &c,
	                const PointT &d, const PointT &query);

	/**
	 * @brief In 3D, test a \p query point is inside the circumcircle of four
	 * points \p a, \p b, \p c and \p d. We assume that the tetrahedron formed by
	 * a, b, c and d has a positive volume. In another word, a, b and c form a
	 * triangle (in CCW order), the vector from a to d must point to the same
	 * direction of the normal.
	 * @return 1: inside, 0: exactly on , -1: outside.
	 */
	Sign operator()(const NT &ax, const NT &ay, const NT &az, const NT &bx,
	                const NT &by, const NT &bz, const NT &cx, const NT &cy,
	                const NT &cz, const NT &dx, const NT &dy, const NT &dz,
	                const NT &queryx, const NT &queryy, const NT &queryz);
};

} // namespace OMC

#ifdef OMC_HAS_IMPL
	#include "GeneralPredicates.inl"
#endif