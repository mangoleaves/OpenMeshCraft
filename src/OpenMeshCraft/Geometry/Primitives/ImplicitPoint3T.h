#pragma once

#include "ExplicitPoint3T.h"
#include "GenericPoint3T.h"
#include "GlobalCachedValues.h"

#include <memory>

namespace OMC {

/// @brief Implicit point defined by the intersection of two coplanar segments
/// in 3D.
template <typename IT_, typename ET_>
class ImplicitPoint3T_SSI : public GenericPoint3T<IT_, ET_>
{
public: /* Types *************************************************************/
	using FT = double;
	using IT = IT_;
	using ET = ET_;

	using EP = ExplicitPoint3T<IT, ET>;
	using IP = ImplicitPoint3T_SSI<IT, ET>;
	using GP = GenericPoint3T<IT, ET>;

	using PointType = typename GP::PointType;

	using GCV = GlobalCachedValues<IT, ET, OnePointCachedValues3<IT, ET>>;

public: /* functions about types *********************************************/
	// Calculates an explicit approximation of the implicit point.
	void get_Explicit(EP &e) const;
	EP   to_Explicit() const;

public: /* Constructors ******************************************************/
	/// @brief init SSI point with segment(a,b) and segment(p,q)
	ImplicitPoint3T_SSI() noexcept;
	ImplicitPoint3T_SSI(const EP &_a, const EP &_b, const EP &_p, const EP &_q,
	                    int _plane) noexcept;
	virtual ~ImplicitPoint3T_SSI() noexcept;

	ImplicitPoint3T_SSI(const IP &rhs) noexcept;
	ImplicitPoint3T_SSI(IP &&rhs) noexcept;

	IP &operator=(const IP &rhs);
	IP &operator=(IP &&rhs);

public: /* Members ***********************************************************/
	const EP &A() const { return *ia; }
	const EP &B() const { return *ib; }
	const EP &P() const { return *ip; }
	const EP &Q() const { return *iq; }

public: /* Lambdas ***********************************************************/
	bool getFilteredLambda(FT &lx, FT &ly, FT &lz, FT &d, FT &mv) const;
	bool getIntervalLambda(IT &lx, IT &ly, IT &lz, IT &d) const;
	void getExactLambda(ET &lx, ET &ly, ET &lz, ET &d) const;
	void getExpansionLambda(FT **lx, int &lx_len, FT **ly, int &ly_len, FT **lz,
	                        int &lz_len, FT **d, int &d_len) const;

	static GCV &gcv() { return global_cached_values; }

private:
	const EP *ia, *ib; // The first segment AB
	const EP *ip, *iq; // The second segment PQ
	// The 2D plane where two segments are not degenerate.
	// Its value means the index of max component of plane's normal.
	// plane == 0 -> YZ plane;
	// plane == 1 -> ZX plane;
	// plane == 2 -> XY plane.
	int       plane;

	static GCV global_cached_values;
};

/// @brief Implicit point defined by the intersection of a line and a plane
template <typename IT_, typename ET_>
class ImplicitPoint3T_LPI : public GenericPoint3T<IT_, ET_>
{
public: /* Types *************************************************************/
	using FT = double;
	using IT = IT_;
	using ET = ET_;

	using EP = ExplicitPoint3T<IT, ET>;
	using IP = ImplicitPoint3T_LPI<IT, ET>;
	using GP = GenericPoint3T<IT, ET>;

	using PointType = typename GP::PointType;

	using GCV = GlobalCachedValues<IT, ET, OnePointCachedValues3<IT, ET>>;

public: /* functions about types *********************************************/
	// Calculates an explicit approximation of the implicit point.
	void get_Explicit(EP &e) const;
	EP   to_Explicit() const;

public: /* Constructors ******************************************************/
	/// @brief init LPI point with line(p,q) and plane(r,s,t)
	ImplicitPoint3T_LPI() noexcept;
	ImplicitPoint3T_LPI(const EP &_p, const EP &_q, const EP &_r, const EP &_s,
	                    const EP &_t) noexcept;
	virtual ~ImplicitPoint3T_LPI() noexcept;

	ImplicitPoint3T_LPI(const IP &rhs) noexcept;
	ImplicitPoint3T_LPI(IP &&rhs) noexcept;

	IP &operator=(const IP &rhs);
	IP &operator=(IP &&rhs);

public: /* Members ***********************************************************/
	const EP &P() const { return *ip; }
	const EP &Q() const { return *iq; }
	const EP &R() const { return *ir; }
	const EP &S() const { return *is; }
	const EP &T() const { return *it; }

public: /* Lambdas ***********************************************************/
	bool getFilteredLambda(FT &lx, FT &ly, FT &lz, FT &d, FT &mv) const;
	bool getIntervalLambda(IT &lx, IT &ly, IT &lz, IT &d) const;
	void getExactLambda(ET &lx, ET &ly, ET &lz, ET &d) const;
	void getExpansionLambda(FT **lx, int &lx_len, FT **ly, int &ly_len, FT **lz,
	                        int &lz_len, FT **d, int &d_len) const;

	static GCV &gcv() { return global_cached_values; }

private:
	const EP *ip, *iq;      // The line
	const EP *ir, *is, *it; // The plane

	static GCV global_cached_values;
};

/// @brief Implicit point defined by the intersection of three planes
template <typename IT_, typename ET_>
class ImplicitPoint3T_TPI : public GenericPoint3T<IT_, ET_>
{
public: /* Types *************************************************************/
	using FT = double;
	using IT = IT_;
	using ET = ET_;

	using EP = ExplicitPoint3T<IT, ET>;
	using IP = ImplicitPoint3T_TPI<IT, ET>;
	using GP = GenericPoint3T<IT, ET>;

	using PointType = typename GP::PointType;

	using GCV = GlobalCachedValues<IT, ET, OnePointCachedValues3<IT, ET>>;

public: /* functions about types *********************************************/
	// Calculates an explicit approximation of the implicit point.
	void get_Explicit(EP &e) const;
	EP   to_Explicit() const;

public: /* Constructors ******************************************************/
	ImplicitPoint3T_TPI() noexcept;
	ImplicitPoint3T_TPI(const EP &_v1, const EP &_v2, const EP &_v3,
	                    const EP &_w1, const EP &_w2, const EP &_w3,
	                    const EP &_u1, const EP &_u2, const EP &_u3) noexcept;
	virtual ~ImplicitPoint3T_TPI() noexcept;

	ImplicitPoint3T_TPI(const IP &rhs) noexcept;
	ImplicitPoint3T_TPI(IP &&rhs) noexcept;

	IP &operator=(const IP &rhs);
	IP &operator=(IP &&rhs);

public: /* Members ***********************************************************/
	const EP &V1() const { return *iv1; }
	const EP &V2() const { return *iv2; }
	const EP &V3() const { return *iv3; }
	const EP &W1() const { return *iw1; }
	const EP &W2() const { return *iw2; }
	const EP &W3() const { return *iw3; }
	const EP &U1() const { return *iu1; }
	const EP &U2() const { return *iu2; }
	const EP &U3() const { return *iu3; }

public: /* Lambdas ***********************************************************/
	bool getFilteredLambda(FT &lx, FT &ly, FT &lz, FT &d, FT &mv) const;
	bool getIntervalLambda(IT &lx, IT &ly, IT &lz, IT &d) const;
	void getExactLambda(ET &lx, ET &ly, ET &lz, ET &d) const;
	void getExpansionLambda(FT **lx, int &lx_len, FT **ly, int &ly_len, FT **lz,
	                        int &lz_len, FT **d, int &d_len) const;

	static GCV &gcv() { return global_cached_values; }

private:
	const EP *iv1, *iv2, *iv3; // Plane 1
	const EP *iw1, *iw2, *iw3; // Plane 2
	const EP *iu1, *iu2, *iu3; // Plane 3

	static GCV global_cached_values;
};

/// @brief Implicit point defined as a linear combination of two points.
/// start with p and end with q: (1-t)p + tq  or  p + t(q-p).
template <typename IT_, typename ET_>
class ImplicitPoint3T_LNC : public GenericPoint3T<IT_, ET_>
{
public: /* Types *************************************************************/
	using FT = double;
	using IT = IT_;
	using ET = ET_;

	using EP = ExplicitPoint3T<IT, ET>;
	using IP = ImplicitPoint3T_LNC<IT, ET>;
	using GP = GenericPoint3T<IT, ET>;

	using PointType = typename GP::PointType;

	using GCV = GlobalCachedValues<IT, ET, OnePointCachedValues3<IT, ET>>;

public: /* functions about types *********************************************/
	// Calculates an explicit approximation of the implicit point.
	void get_Explicit(EP &e) const;
	EP   to_Explicit() const;

public: /* Constructors ******************************************************/
	ImplicitPoint3T_LNC() noexcept;
	ImplicitPoint3T_LNC(const EP &_p, const EP &_q, const double _t) noexcept;

	virtual ~ImplicitPoint3T_LNC() noexcept;

	ImplicitPoint3T_LNC(const ImplicitPoint3T_LNC &rhs) noexcept;
	ImplicitPoint3T_LNC(ImplicitPoint3T_LNC &&rhs) noexcept;

	IP &operator=(const IP &rhs);
	IP &operator=(IP &&rhs);

	const EP     &P() const { return *ip; }
	const EP     &Q() const { return *iq; }
	const double &T() const { return t; }

	bool getFilteredLambda(FT &lx, FT &ly, FT &lz, FT &d, FT &mv) const;
	bool getIntervalLambda(IT &lx, IT &ly, IT &lz, IT &d) const;
	void getExactLambda(ET &lx, ET &ly, ET &lz, ET &d) const;
	void getExpansionLambda(FT **lx, int &lx_len, FT **ly, int &ly_len, FT **lz,
	                        int &lz_len, FT **d, int &d_len) const;

	static GCV &gcv() { return global_cached_values; }

private:
	const EP *ip, *iq; // The two points
	double    t;       // The parameter (0 = ip, 1 = iq)

	static GCV global_cached_values;
};

// static member variables
template <typename IT_, typename ET_>
typename ImplicitPoint3T_SSI<IT_, ET_>::GCV
  ImplicitPoint3T_SSI<IT_, ET_>::global_cached_values =
    ImplicitPoint3T_SSI<IT_, ET_>::GCV();

template <typename IT_, typename ET_>
typename ImplicitPoint3T_LPI<IT_, ET_>::GCV
  ImplicitPoint3T_LPI<IT_, ET_>::global_cached_values =
    ImplicitPoint3T_LPI<IT_, ET_>::GCV();

template <typename IT_, typename ET_>
typename ImplicitPoint3T_TPI<IT_, ET_>::GCV
  ImplicitPoint3T_TPI<IT_, ET_>::global_cached_values =
    ImplicitPoint3T_TPI<IT_, ET_>::GCV();

template <typename IT_, typename ET_>
typename ImplicitPoint3T_LNC<IT_, ET_>::GCV
  ImplicitPoint3T_LNC<IT_, ET_>::global_cached_values =
    ImplicitPoint3T_LNC<IT_, ET_>::GCV();

inline void normalizeLambda3D(double *lx, int &lxl, double *ly, int &lyl,
                              double *lz, int &lzl, double *d, int &dl)
{
	expansionObject o;
	double          maxd, maxsd, ad, aad;
	maxsd = o.To_Double(lxl, lx);
	maxd  = fabs(maxsd);
	if ((aad = fabs((ad = o.To_Double(lyl, ly)))) > maxd)
	{
		maxd  = aad;
		maxsd = ad;
	}
	if ((aad = fabs((ad = o.To_Double(lzl, lz)))) > maxd)
	{
		maxd  = aad;
		maxsd = ad;
	}
	if ((aad = fabs((ad = o.To_Double(dl, d)))) > maxd)
	{
		maxd  = aad;
		maxsd = ad;
	}

	int e;
	frexp(maxsd, &e);
	const double m = ldexp(2, -e);

	o.ExactScale(lxl, lx, m);
	o.ExactScale(lyl, ly, m);
	o.ExactScale(lzl, lz, m);
	o.ExactScale(dl, d, m);
}

} // namespace OMC

#ifdef OMC_HAS_IMPL
	#include "ImplicitPoint3T_SSI.inl"
	#include "ImplicitPoint3T_LPI.inl"
	#include "ImplicitPoint3T_TPI.inl"
	#include "ImplicitPoint3T_LNC.inl"
#endif