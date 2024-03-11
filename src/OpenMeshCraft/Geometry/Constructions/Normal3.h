#pragma once

namespace OMC {

template <typename Kernel>
class ConstructNormal3K
{
public:
	using VecT    = typename Kernel::Vec3;
	using EPointT = typename Kernel::EPoint3;
	using GPointT = typename Kernel::GPoint3;

	using TriangleT = typename Kernel::Triangle3;

	using ToEP = typename Kernel::ToEP;

public:
	template <typename GPT,
	          typename = std::enable_if_t<std::is_same_v<GPT, GPointT> &&
	                                      !std::is_same_v<GPT, EPointT>>>
	VecT operator()(const GPT &v0, const GPT &v1, const GPT &v2) const
	{
		return (ToEP()(v1) - ToEP()(v0)).cross(ToEP()(v2) - ToEP()(v0));
	}

	VecT operator()(const EPointT &v0, const EPointT &v1, const EPointT &v2) const
	{
		return (v1 - v0).cross(v2 - v0);
	}

	VecT operator()(const TriangleT &tri) const
	{
		return operator()(tri.v0(), tri.v1(), tri.v2());
	}
};

} // namespace OMC