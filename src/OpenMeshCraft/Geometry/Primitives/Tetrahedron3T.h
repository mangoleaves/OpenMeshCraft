#pragma once

#include "Primitive3.h"

namespace OMC {

/**
 * @brief Tetrahedron in 3D that contains four points as its vertices.
 * @tparam NT The number type.
 * @tparam PointT The point type.
 */
template <typename _NT, typename _PointT>
class Tetrahedron3T : public Primitive3<_NT>
{
public:
  using NT           = _NT;
  using PointT       = _PointT;
  using TetrahedronT = Tetrahedron3T<NT, PointT>;

public:
  /**
   * @brief Construct tetrahedron by default.
   */
  Tetrahedron3T() {}

  /**
   * @brief Construct tetrahedron from four points.
   */
  Tetrahedron3T(const PointT &v0, const PointT &v1, const PointT &v2, const PointT &v3)
    : m_v0(v0)
    , m_v1(v1)
    , m_v2(v2)
    , m_v3(v3)
  {
  }

  /// @name Access data of tetrahedron.
  /// @{
  /**
   * @brief Get the i-th vertex.
   * @warning Is "pointer cast + offset" safe?
   */
  PointT       &operator[](size_t i) { return (&m_v0)[i]; }
  const PointT &operator[](size_t i) const { return (&m_v0)[i]; }

  PointT       &v0() { return m_v0; }
  const PointT &v0() const { return m_v0; }
  PointT       &v1() { return m_v1; }
  const PointT &v1() const { return m_v1; }
  PointT       &v2() { return m_v2; }
  const PointT &v2() const { return m_v2; }
  PointT       &v3() { return m_v3; }
  const PointT &v3() const { return m_v3; }
  /// @}

protected:
  PointT m_v0;
  PointT m_v1;
  PointT m_v2;
  PointT m_v3;
};

} // namespace OMC