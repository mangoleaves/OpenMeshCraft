#pragma once

#include "Primitive3.h"

#include <cmath>
#include <type_traits>

namespace OMC {

/**
 * @brief Implicit normal cone in 3D, normal cone's apex is at the origin and
 * goes infinitely far along its center axis. We store the direction of center
 * axis and the apex angle.
 * @tparam NT The number type.
 * @tparam VecT The vector type.
 * @note center axis is supposed to be normalized.
 * @note apex angle is represented in radian
 */
template <typename _NT, typename _VecT>
class NormalCone3T : public Primitive3<_NT>
{
public:
	using NT          = _NT;
	using VecT        = _VecT;
	using NormalConeT = NormalCone3T<NT, VecT>;

	static_assert(std::is_floating_point_v<NT>);

public:
	/**
	 * @brief Construct a NormalCone3T object by default.
	 */
	NormalCone3T() {}

	/**
	 * @brief Construct a NormalCone3T object from given \p center_axis and \p
	 * apex_angle.
	 */
	NormalCone3T(const VecT &center_axis, const NT &apex_angle)
	  : m_center_axis(center_axis)
	  , m_apex_angle(apex_angle)
	{
	}

	/// @note center axis is supposed to be normalized.
	VecT       &center_axis() { return m_center_axis; }
	/// @note center axis is supposed to be normalized.
	const VecT &center_axis() const { return m_center_axis; }

	/// @note apex angle is represented in radian
	NT       &apex_angle() { return m_apex_angle; }
	/// @note apex angle is represented in radian
	const NT &apex_angle() const { return m_apex_angle; }

	void set_apex_angle_degree(NT angle) { m_apex_angle = angle / 180. * M_PI; }
	NT   apex_angle_degree() const { return m_apex_angle / M_PI * 180.; }

	/**
	 * @brief Calculate the maximal angle difference between this cone and a
	 * normal.
	 */
	NT max_angle_delta(const VecT &normal) const;

	/**
	 * @brief Calculate the maximal angle difference between two cones.
	 * @param cone another normal cone
	 */
	NT max_angle_delta(const NormalCone3T &cone) const;

	/**
	 * @brief merge this normal cone with a normal.
	 * @param normal a unit normal.
	 */
	void merge(const VecT &normal);

	/**
	 * @brief merge two normal cones to one.
	 * @param cone another normal cone.
	 */
	void merge(const NormalCone3T &cone);

private:
	//        /
	//      /
	//    /
	//  /   apex_angle
	// ------------------ center_axis
	// \    apex_angle
	//   \ 
  //     \ 
  //       \ 

	VecT m_center_axis;
	/// @note apex angle is represented in radian
	NT   m_apex_angle;
};

} // namespace OMC

#ifdef OMC_HAS_IMPL
	#include "NormalCone3T.inl"
#endif