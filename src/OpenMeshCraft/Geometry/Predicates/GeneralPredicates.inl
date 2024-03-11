#pragma once

#include "GeneralPredicates.h"

#include <tuple>

namespace OMC {

template <typename NT>
Sign DotProductSign2D_GNR<NT>::operator()(const PointT &p, const PointT &r,
                                          const PointT &q)
{
	return OMC::sign((p - q).dot(r - q));
}

template <typename NT>
Sign DotProductSign2D_GNR<NT>::operator()(const PointT &p, const PointT &r,
                                          const PointT &q, const PointT &s)
{
	return OMC::sign((p - q).dot(r - s));
}

template <typename NT>
Sign DotProductSign3D_GNR<NT>::operator()(const PointT &p, const PointT &r,
                                          const PointT &q)
{
	return OMC::sign((p - q).dot(r - q));
}

template <typename NT>
Sign DotProductSign3D_GNR<NT>::operator()(const PointT &p, const PointT &r,
                                          const PointT &q, const PointT &s)
{
	return OMC::sign((p - q).dot(r - s));
}

template <typename NT>
Sign DotProductSignOn2D_GNR<NT>::on_xy(const PointT &p, const PointT &r,
                                       const PointT &q)
{
	return OMC::sign((p.x() - q.x()) * (r.x() - q.x()) +
	                 (p.y() - q.y()) * (r.y() - q.y()));
}

template <typename NT>
Sign DotProductSignOn2D_GNR<NT>::on_xy(const PointT &p, const PointT &r,
                                       const PointT &q, const PointT &s)
{
	return OMC::sign((p.x() - q.x()) * (r.x() - s.x()) +
	                 (p.y() - q.y()) * (r.y() - s.y()));
}

template <typename NT>
Sign DotProductSignOn2D_GNR<NT>::on_yz(const PointT &p, const PointT &r,
                                       const PointT &q)
{
	return OMC::sign((p.y() - q.y()) * (r.y() - q.y()) +
	                 (p.z() - q.z()) * (r.z() - q.z()));
}

template <typename NT>
Sign DotProductSignOn2D_GNR<NT>::on_yz(const PointT &p, const PointT &r,
                                       const PointT &q, const PointT &s)
{
	return OMC::sign((p.y() - q.y()) * (r.y() - s.y()) +
	                 (p.z() - q.z()) * (r.z() - s.z()));
}

template <typename NT>
Sign DotProductSignOn2D_GNR<NT>::on_zx(const PointT &p, const PointT &r,
                                       const PointT &q)
{
	return OMC::sign((p.x() - q.x()) * (r.x() - q.x()) +
	                 (p.z() - q.z()) * (r.z() - q.z()));
}

template <typename NT>
Sign DotProductSignOn2D_GNR<NT>::on_zx(const PointT &p, const PointT &r,
                                       const PointT &q, const PointT &s)
{
	return OMC::sign((p.x() - q.x()) * (r.x() - s.x()) +
	                 (p.z() - q.z()) * (r.z() - s.z()));
}

template <typename NT>
Sign Orient2D_GNR<NT>::operator()(const PointT &p, const PointT &q,
                                  const PointT &query)
{
	return OMC::sign((q - p).determinant(query - p));
}

template <typename NT>
Sign Orient2D_GNR<NT>::operator()(const NT &px, const NT &py, const NT &qx,
                                  const NT &qy, const NT &query_x,
                                  const NT &query_y)
{
	return OMC::sign((qx - px) * (query_y - py) - (qy - py) * (query_x - px));
}

template <typename NT>
Sign Orient3D_GNR<NT>::operator()(const PointT &p, const PointT &q,
                                  const PointT &r, const PointT &query)
{
	return OMC::sign((query - p).dot((q - p).cross(r - p)));
}

template <typename NT>
Sign OrientOn2D_GNR<NT>::on_xy(const PointT &a, const PointT &b,
                               const PointT &c)
{
	return Orient2D()(a.x(), a.y(), b.x(), b.y(), c.x(), c.y());
}

template <typename NT>
Sign OrientOn2D_GNR<NT>::on_yz(const PointT &a, const PointT &b,
                               const PointT &c)
{
	return Orient2D()(a.y(), a.z(), b.y(), b.z(), c.y(), c.z());
}

template <typename NT>
Sign OrientOn2D_GNR<NT>::on_zx(const PointT &a, const PointT &b,
                               const PointT &c)
{
	return Orient2D()(a.z(), a.x(), b.z(), b.x(), c.z(), c.x());
}

template <typename NT>
bool CollinearPoints2D_GNR<NT>::operator()(const PointT &p, const PointT &q,
                                           const PointT &r)
{
	return Orient2D()(p, q, r) == Sign::ZERO;
}

template <typename NT>
Sign LessThan2D_GNR<NT>::operator()(const PointT &p, const PointT &q)
{
	Sign ret = on_x(p, q);
	if (is_sign_posneg(ret))
		return ret;
	return on_y(p, q);
}

template <typename NT>
Sign LessThan2D_GNR<NT>::on_x(const PointT &p, const PointT &q)
{
	return static_cast<Sign>(((p.x() > q.x()) - (p.x() < q.x())));
}

template <typename NT>
Sign LessThan2D_GNR<NT>::on_y(const PointT &p, const PointT &q)
{
	return static_cast<Sign>(((p.y() > q.y()) - (p.y() < q.y())));
}

template <typename NT>
Sign LessThan3D_GNR<NT>::operator()(const PointT &p, const PointT &q)
{
	Sign ret = on_x(p, q);
	if (is_sign_posneg(ret))
		return ret;
	ret = on_y(p, q);
	if (is_sign_posneg(ret))
		return ret;
	return on_z(p, q);
}

template <typename NT>
Sign LessThan3D_GNR<NT>::on_x(const PointT &p, const PointT &q)
{
	return static_cast<Sign>(((p.x() > q.x()) - (p.x() < q.x())));
}

template <typename NT>
Sign LessThan3D_GNR<NT>::on_y(const PointT &p, const PointT &q)
{
	return static_cast<Sign>(((p.y() > q.y()) - (p.y() < q.y())));
}

template <typename NT>
Sign LessThan3D_GNR<NT>::on_z(const PointT &p, const PointT &q)
{
	return static_cast<Sign>(((p.z() > q.z()) - (p.z() < q.z())));
}

template <typename NT>
bool CollinearPoints3D_GNR<NT>::operator()(const PointT &p, const PointT &q,
                                           const PointT &r)
{
	Orient2D orient2d;
	return orient2d(p.x(), p.y(), q.x(), q.y(), r.x(), r.y()) == Sign::ZERO &&
	       orient2d(p.x(), p.z(), q.x(), q.z(), r.x(), r.z()) == Sign::ZERO &&
	       orient2d(p.y(), p.z(), q.y(), q.z(), r.y(), r.z()) == Sign::ZERO;
}

template <typename NT>
bool CollinearOrdered3D_GNR<NT>::operator()(const PointT &p, const PointT &q,
                                            const PointT &r)
{
	if (p.x() < q.x())
		return !(r.x() < q.x());
	if (q.x() < p.x())
		return !(q.x() < r.x());
	if (p.y() < q.y())
		return !(r.y() < q.y());
	if (q.y() < p.y())
		return !(q.y() < r.y());
	if (p.z() < q.z())
		return !(r.z() < q.z());
	if (q.z() < p.z())
		return !(q.z() < r.z());
	return true; // p==q
}

template <typename NT>
auto CollinearSort3D_GNR<NT>::operator()(const PointT &p, const PointT &q,
                                         const PointT &r)
  -> std::tuple<const PointT &, const PointT &, const PointT &>
{
	using CR = const PointT &;
	using CP = const PointT *;

#define SORT_ON_AXIS(axis)   \
	if (a->axis() > b->axis()) \
		std::swap(a, b);         \
	if (b->axis() > c->axis()) \
		std::swap(b, c);         \
	if (a->axis() > c->axis()) \
		std::swap(a, c);

	CP a = &p, b = &q, c = &r;
	if (p.x() != q.x() || p.x() != r.x() || q.x() != r.x())
	{
		SORT_ON_AXIS(x);
	}
	else if (p.y() != q.y() || p.y() != r.y() || q.y() != r.y())
	{
		SORT_ON_AXIS(y);
	}
	else if (p.z() != q.z() || p.z() != r.z() || q.z() != r.z())
	{
		SORT_ON_AXIS(z);
	}
#undef SORT_ON_AXIS

	return std::make_tuple<CR, CR, CR>(*a, *b, *c);
}

template <typename NT>
bool CoplanarPoints3D_GNR<NT>::operator()(const PointT &p, const PointT &q,
                                          const PointT &r, const PointT &s)
{
	return Orient3D()(p, q, r, s) == Sign::ZERO;
}

template <typename NT>
Sign InCircle_GNR<NT>::operator()(const PointT &p, const PointT &q,
                                  const PointT &r, const PointT &query)
{
	return operator()(p.x(), p.y(), q.x(), q.y(), r.x(), r.y(), query.x(),
	                  query.y());
}

template <typename NT>
Sign InCircle_GNR<NT>::operator()(const NT &px, const NT &py, const NT &qx,
                                  const NT &qy, const NT &rx, const NT &ry,
                                  const NT &queryx, const NT &queryy)
{
	// Lift p, q, r and query to parabolla p', q', r' and query' (x, y,
	// x^2+y^2), then test the orientation of query' with respect to the
	// hyperplane formed by p',q' and r'.

	// More efficiently, lift (p-query), (q-query) and (r-query), then
	// test the orientation of origin with respect to the hyperplane.
	NT px_ = px - queryx;
	NT py_ = py - queryy;
	NT pz_ = px_ * px_ + py_ * py_;
	NT qx_ = qx - queryx;
	NT qy_ = qy - queryy;
	NT qz_ = qx_ * qx_ + qy_ * qy_;
	NT rx_ = rx - queryx;
	NT ry_ = ry - queryy;
	NT rz_ = rx_ * rx_ + ry_ * ry_;

	return OMC::sign(pz_ * (qx_ * ry_ - qy_ * rx_) -
	                 qz_ * (px_ * ry_ - py_ * rx_) +
	                 rz_ * (px_ * qy_ - py_ * qx_));
}

template <typename NT>
Sign InSphere_GNR<NT>::operator()(const PointT &a, const PointT &b,
                                  const PointT &c, const PointT &d,
                                  const PointT &query)
{
	return operator()(a.x(), a.y(), a.z(), b.x(), b.y(), b.z(), c.x(), c.y(),
	                  c.z(), d.x(), d.y(), d.z(), query.x(), query.y(),
	                  query.z());
}

template <typename NT>
Sign InSphere_GNR<NT>::operator()(const NT &ax, const NT &ay, const NT &az,
                                  const NT &bx, const NT &by, const NT &bz,
                                  const NT &cx, const NT &cy, const NT &cz,
                                  const NT &dx, const NT &dy, const NT &dz,
                                  const NT &queryx, const NT &queryy,
                                  const NT &queryz)
{
	// Lift a, b, c, d and query to parabolla, then test the orientation of lifted
	// query with respect to the hyperplane formed by lifted a, b, c and d.

	// More efficiently, lift (a-query), (b-query), (c-query) and (d-query), then
	// test the orientation of origin with respect to the hyperplane.
	NT ax_   = ax - queryx;
	NT ay_   = ay - queryy;
	NT az_   = az - queryz;
	NT bx_   = bx - queryx;
	NT by_   = by - queryy;
	NT bz_   = bz - queryz;
	NT cx_   = cx - queryx;
	NT cy_   = cy - queryy;
	NT cz_   = cz - queryz;
	NT dx_   = dx - queryx;
	NT dy_   = dy - queryy;
	NT dz_   = dz - queryz;
	NT ab    = ax_ * by_ - bx_ * ay_;
	NT bc    = bx_ * cy_ - cx_ * by_;
	NT cd    = cx_ * dy_ - dx_ * cy_;
	NT da    = dx_ * ay_ - ax_ * dy_;
	NT ac    = ax_ * cy_ - cx_ * ay_;
	NT bd    = bx_ * dy_ - dx_ * by_;
	NT abc   = az_ * bc + cz_ * ab - bz_ * ac;
	NT bcd   = bz_ * cd + dz_ * bc - cz_ * bd;
	NT cda   = cz_ * da + az_ * cd + dz_ * ac;
	NT dab   = dz_ * ab + bz_ * da + az_ * bd;
	NT alift = ax_ * ax_ + ay_ * ay_ + az_ * az_;
	NT blift = bx_ * bx_ + by_ * by_ + bz_ * bz_;
	NT clift = cx_ * cx_ + cy_ * cy_ + cz_ * cz_;
	NT dlift = dx_ * dx_ + dy_ * dy_ + dz_ * dz_;

	NT det = clift * dab - dlift * abc + alift * bcd - blift * cda;
	return OMC::sign(det);
}

} // namespace OMC