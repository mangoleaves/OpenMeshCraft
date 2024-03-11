#pragma once

#include "NormalCone3T.h"

namespace OMC {

template <typename NT, typename VecT>
auto NormalCone3T<NT, VecT>::max_angle_delta(const VecT &normal) const -> NT
{
	// angle delta between axis and normal
	NT cos_angle_delta = m_center_axis.dot(normal);
	NT angle_delta =
	  cos_angle_delta >= 1.0
	    ? 0.0
	    : (cos_angle_delta <= -1.0 ? M_PI : std::acos(cos_angle_delta));
	return angle_delta + m_apex_angle;
}

template <typename NT, typename VecT>
auto NormalCone3T<NT, VecT>::max_angle_delta(const NormalCone3T &cone) const
  -> NT
{
	// angle delta between two axes
	NT cos_angle_delta = m_center_axis.dot(cone.m_center_axis);
	NT angle_delta =
	  cos_angle_delta >= NT(1.0)
	    ? NT(0.0)
	    : (cos_angle_delta <= NT(-1.0) ? NT(M_PI) : std::acos(cos_angle_delta));
	NT side_angle0 = std::max(m_apex_angle - angle_delta, cone.m_apex_angle);
	NT side_angle1 = std::max(cone.m_apex_angle - angle_delta, m_apex_angle);

	return angle_delta + side_angle0 + side_angle1;
}

template <typename NT, typename VecT>
void NormalCone3T<NT, VecT>::merge(const VecT &normal)
{
	NT cos_angle_delta = m_center_axis.dot(normal);
	if (std::fabs(cos_angle_delta) < 0.99999)
	{
		// calculate new apex angle.
		NT angle_delta = std::acos(cos_angle_delta);
		NT min_angle   = -m_apex_angle;
		NT max_angle   = std::max(m_apex_angle, angle_delta);
		m_apex_angle   = (max_angle - min_angle) * NT(0.5);
		// angle between m_center_axis and new axis
		NT axis_angle  = (max_angle + min_angle) * NT(0.5);
		// new axis by SLERP
		m_center_axis  = ((m_center_axis * std::sin(angle_delta - axis_angle) +
                      normal * std::sin(axis_angle)) /
                     std::sin(angle_delta));
	}
	else
	{
		// axes are almost parallel
		if (cos_angle_delta < 0.0)
			// axes point in opposite directions
			m_apex_angle = NT(M_PI);
	}
}

template <typename NT, typename VecT>
void NormalCone3T<NT, VecT>::merge(const NormalCone3T &cone)
{
	NT cos_angle_delta = m_center_axis.dot(cone.m_center_axis);
	if (std::fabs(cos_angle_delta) < 0.99999)
	{
		// calculate new apex angle.
		NT angle_delta = std::acos(cos_angle_delta);
		NT min_angle   = std::min(-m_apex_angle, angle_delta - cone.m_apex_angle);
		NT max_angle   = std::max(m_apex_angle, angle_delta + cone.m_apex_angle);
		m_apex_angle   = (max_angle - min_angle) * NT(0.5);
		// angle between m_center_axis and new axis
		NT axis_angle  = (min_angle + max_angle) * NT(0.5);
		// new axis by SLERP (interpolate quaternion)
		m_center_axis  = ((m_center_axis * std::sin(angle_delta - axis_angle) +
                      cone.m_center_axis * std::sin(axis_angle)) /
                     std::sin(angle_delta));
	}
	else
	{
		// axes are almost parallel
		if (cos_angle_delta > 0.0)
			// axes point in same direction
			m_apex_angle = std::max(m_apex_angle, cone.m_apex_angle);
		else
			// axes point in opposite directions
			m_apex_angle = NT(M_PI);
	}
}

} // namespace OMC