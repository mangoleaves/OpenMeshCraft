#pragma once

#include "Primitive2.h"

#include "OpenMeshCraft/NumberTypes/NumberUtils.h"

namespace OMC {

/**
 * @brief Axis-Aligned Bounding Box in 2D that contains two points as min/max
 * bounds.
 * @tparam NT The number type.
 * @tparam VecT The vector type.
 * @tparam PointT The point type.
 */
template <typename _NT, typename _VecT, typename _PointT>
class BoundingBox2T : public Primitive2<_NT>
{
public:
	using NT     = _NT;
	using VecT   = _VecT;
	using PointT = _PointT;
	using BboxT  = BoundingBox2T<NT, VecT, PointT>;

public:
	/**
	 * @brief Construct a new BoundingBox2T object.
	 *  By default, its bounds are set to numeric limits.
	 */
	BoundingBox2T()
	  : m_min_bound(NumericLimits<double>::max())
	  , m_max_bound(NumericLimits<double>::lowest())
	{
	}

	BoundingBox2T(const BoundingBox2T &)            = default;
	BoundingBox2T(BoundingBox2T &&)                 = default;
	BoundingBox2T &operator=(const BoundingBox2T &) = default;
	BoundingBox2T &operator=(BoundingBox2T &&)      = default;

	/**
	 * @brief Construct a new BoundingBox2T object from given bounds \p minB and
	 * \p maxB .
	 * @param[in] minB minimal bound
	 * @param[in] maxB maximal bound
	 */
	BoundingBox2T(const PointT &minB, const PointT &maxB)
	  : m_min_bound(minB)
	  , m_max_bound(maxB)
	{
	}

	/**
	 * @brief Construct a new BoundingBox2T object from given point \p point .
	 * The point are set as two bounds.
	 * @param[in] point given point
	 */
	BoundingBox2T(const PointT &point)
	  : m_min_bound(point)
	  , m_max_bound(point)
	{
	}

	/**
	 * @brief Construct a new BoundingBox2T object from given points \p points .
	 * The result is the minimal bounding box that bounds \p points .
	 * @param[in] points a set of points, stored in STL container.
	 */
	template <template <typename, typename> class STL_Container>
	BoundingBox2T(const STL_Container<PointT, std::allocator<PointT>> &points)
	  : BoundingBox2T()
	{
		for (const auto &p : points)
			*this += p;
	}

	/**
	 * @brief Check if two bounding boxes are same.
	 * If their min/max bounds are same respectively, they are same.
	 * @param b the other bounding box.
	 * @return true if they are same, false if they are different.
	 */
	inline bool operator==(const BboxT &b) const;

	/// @name Access min/max bounds.
	/// @{
	PointT       &min_bound() { return m_min_bound; }
	PointT       &max_bound() { return m_max_bound; }
	const PointT &min_bound() const { return m_min_bound; }
	const PointT &max_bound() const { return m_max_bound; }
	NT           &min_coord(size_t i) { return m_min_bound[i]; }
	NT           &max_coord(size_t i) { return m_max_bound[i]; }
	const NT     &min_coord(size_t i) const { return m_min_bound[i]; }
	const NT     &max_coord(size_t i) const { return m_max_bound[i]; }
	/// @}

	/**
	 * @brief Calculate a bounding box that is the union of \p this and \p b .
	 * @param[in] b the given box.
	 * @return BoundingBox2T the result of union.
	 */
	inline BboxT operator+(const PointT &p) const;

	/**
	 * @brief Set this bounding box \p this to the union of \p this and \p b .
	 * @param[in] b the given box.
	 * @return BoundingBox2T \p this .
	 */
	inline BboxT &operator+=(const PointT &p);

	/**
	 * @brief Calculate a bounding box that is the union of \p this and \p b .
	 * @param[in] p the given point.
	 * @return BoundingBox2T the result of union.
	 */
	inline BboxT operator+(const BboxT &b) const;

	/**
	 * @brief Set this bounding box \p this to the union of \p this and \p b .
	 * @param[in] p the given point.
	 * @return BoundingBox2T \p this .
	 */
	inline BboxT &operator+=(const BboxT &b);

	/**
	 * @brief Enlarge the bounding box by \p offset .
	 * The min bound is subtracted by \p offset and the max bound is added by \p
	 * offset.
	 * @param[in] offset
	 * @note \p offset is allowed to contain negative values.
	 */
	inline void enlarge(const VecT &offset);

	/**
	 * @brief Get the longest axis of this bounding box.
	 * @return The dimension where longest axis is in.
	 * 0, 1, 2 -> x, y, z.
	 */
	inline size_t longest_axis() const;

private:
	PointT m_min_bound, m_max_bound;
};

} // namespace OMC

#ifdef OMC_HAS_IMPL
	#include "BoundingBox2T.inl"
#endif