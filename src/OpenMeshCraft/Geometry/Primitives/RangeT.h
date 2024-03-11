#pragma once

namespace OMC {

/**
 * @brief A range (aka interval) in 1D.
 * @tparam NT number type
 */
template <class NT>
class RangeT
{
public:
	RangeT() {}
	RangeT(const NT &low, const NT &high)
	  : m_low(low)
	  , m_high(high)
	{
	}

	NT   low() const { return m_low; }
	void set_low(const NT &low) { m_low = low; }

	NT   high() const { return m_high; }
	void set_high(const NT &high) { m_high = high; }

	int contains(const NT &value) const;
	int contains(const RangeT<NT> &range) const;

	int operator==(const RangeT<NT> &range) const;

private:
	NT m_low;
	NT m_high;
};

template <class NT>
int RangeT<NT>::contains(const NT &value) const
{
	return ((value >= m_low) && (value <= m_high));
}

template <class NT>
int RangeT<NT>::contains(const RangeT<NT> &range) const
{
	return ((range.low() >= m_low) && (range.High() <= m_high));
}

template <class NT>
int RangeT<NT>::operator==(const RangeT<NT> &range) const
{
	return ((range.low() == m_low) && (range.High() == m_high));
}

} // namespace OMC