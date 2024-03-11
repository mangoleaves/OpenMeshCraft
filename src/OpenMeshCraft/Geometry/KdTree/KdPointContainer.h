#pragma once

namespace OMC {

template <typename KdTraits>
class KdPointContainer
{
public:
	using NT = typename KdTraits::NT;

	using KdPointPtr      = typename KdTraits::KdPointPtr;
	using KdPointPtrsIter = typename KdTraits::KdPointPtrsIter;
	using KdBox           = typename KdTraits::KdBox;

public:
	KdPointContainer() {}

	KdPointContainer(KdPointPtrsIter begin, KdPointPtrsIter end);

	inline size_t          size() const { return std::distance(m_begin, m_end); }
	inline KdPointPtrsIter begin() const { return m_begin; }
	inline KdPointPtrsIter end() const { return m_end; }
	inline const KdBox    &bbox() const { return m_bbox; }
	inline const KdBox    &tbox() const { return m_tbox; }

	/**
	 * @brief Split points in \p this to two containers \p this and \p c_low.
	 * Points are split by given separate dimension and separate value.
	 * @param c_low The other point container storing points.
	 * @param sep_dim separate dimension.
	 * @param sep_val separate value.
	 * @param sliding whether sliding points.
	 */
	void split(KdPointContainer &c_low, size_t &sep_dim, NT &sep_val,
	           bool sliding = false);

private:
	KdPointPtrsIter m_begin, m_end;
	KdBox           m_bbox, m_tbox;
	size_t          m_build_dim;
};

} // namespace OMC

#ifdef OMC_HAS_IMPL
	#include "KdPointContainer.inl"
#endif