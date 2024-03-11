#pragma once

#include "ExplicitPoint2T.h"
#include "GenericPoint2T.h"
#include "GlobalCachedValues.h"

#include <memory>

namespace OMC {

/// @brief Implicit exact 2D point defined by the intersection of two lines.
template <typename IT_, typename ET_>
class ImplicitPoint2T_SSI : public GenericPoint2T<IT_, ET_>
{
public: /* types *************************************************************/
	using FT = double;
	using IT = IT_;
	using ET = ET_;

	using EP = ExplicitPoint2T<IT, ET>;
	using IP = ImplicitPoint2T_SSI<IT, ET>;
	using GP = GenericPoint2T<IT, ET>;

	using PointType = typename GP::PointType;

	using GCV = GlobalCachedValues<IT, ET, OnePointCachedValues2<IT, ET>>;

public: /* functions about types *********************************************/
	void get_Explicit(EP &e) const;
	EP   to_Explicit() const;

public: /* Constructors ******************************************************/
	ImplicitPoint2T_SSI() noexcept;
	/// @brief init SSI point with segment(l11, l12) and segment(l21, l22)
	ImplicitPoint2T_SSI(const EP &l11, const EP &l12, const EP &l21,
	                    const EP &l22) noexcept;
	virtual ~ImplicitPoint2T_SSI() noexcept;

	ImplicitPoint2T_SSI(const IP &rhs) noexcept;
	ImplicitPoint2T_SSI(IP &&rhs) noexcept;

	IP &operator=(const IP &rhs);
	IP &operator=(IP &&rhs);

public: /* Members ***********************************************************/
	const EP &L1_1() const { return *l1_1; }
	const EP &L1_2() const { return *l1_2; }
	const EP &L2_1() const { return *l2_1; }
	const EP &L2_2() const { return *l2_2; }

public: /* Lambdas ***********************************************************/
	bool getFilteredLambda(FT &lx, FT &ly, FT &d, FT &mv) const;
	bool getIntervalLambda(IT &lx, IT &ly, IT &d) const;
	void getExactLambda(ET &lx, ET &ly, ET &d) const;
	void getExpansionLambda(FT **lx, int &lx_len, FT **ly, int &ly_len, FT **d,
	                        int &d_len) const;

	static GCV &gcv() { return global_cached_values; }

private:
	const EP *l1_1, *l1_2; // the first segment
	const EP *l2_1, *l2_2; // the second segment

	static GCV global_cached_values;
};

/// static member variables
template <typename IT_, typename ET_>
typename ImplicitPoint2T_SSI<IT_, ET_>::GCV
  ImplicitPoint2T_SSI<IT_, ET_>::global_cached_values =
    ImplicitPoint2T_SSI<IT_, ET_>::GCV();

} // namespace OMC

#ifdef OMC_HAS_IMPL
	#include "ImplicitPoint2T.inl"
#endif