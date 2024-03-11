#include "BoostMultiprecision.h"
#include "gmp.h"
#include "mpfr.h"

namespace OMC {

using Pair = std::pair<double, double>;

Pair BoostMpToInterval::operator()(const BoostRational &x)
{
#if MPFR_VERSION_MAJOR >= 3
	mpfr_exp_t emin = mpfr_get_emin();
	mpfr_set_emin(-1073);
	MPFR_DECL_INIT(y, 53); /* Assume IEEE-754 */
	int r    = mpfr_set_q(y, x.backend().data(), MPFR_RNDA);
	r        = mpfr_subnormalize(y, r, MPFR_RNDA); /* Round subnormals */
	double i = mpfr_get_d(y, MPFR_RNDA);           /* EXACT but can overflow */
	mpfr_set_emin(emin); /* Restore old value, users may care */
	// With mpfr_set_emax(1024) we could drop the is_finite test
	if (r == 0 && std::isfinite(i))
		return std::pair<double, double>(i, i);
	else
	{
		double s = nextafter(i, 0);
		if (i < 0)
			return std::pair<double, double>(i, s);
		else
			return std::pair<double, double>(s, i);
	}
#else
	mpfr_t y;
	mpfr_init2(y, 53); /* Assume IEEE-754 */
	mpfr_set_q(y, x.backend().data(), GMP_RNDD);
	double i = mpfr_get_d(y, GMP_RNDD); /* EXACT but can overflow */
	mpfr_set_q(y, x.backend().data(), GMP_RNDU);
	double s = mpfr_get_d(y, GMP_RNDU); /* EXACT but can overflow */
	mpfr_clear(y);
	return std::pair<double, double>(i, s);
#endif
}

std::pair<double, double> BoostMpToInterval::operator()(const BoostFloat &x)
{
#if MPFR_VERSION_MAJOR >= 3
	MPFR_DECL_INIT(y, 53); /* Assume IEEE-754 */
	int    r = mpfr_set_f(y, x.backend().data(), MPFR_RNDA);
	double i = mpfr_get_d(y, MPFR_RNDA);

	if (r == 0 && std::isfinite(i))
		return std::pair<double, double>(i, i);
	else
	{
		double s = std::nextafter(i, 0);
		if (i < 0)
			return std::pair<double, double>(i, s);
		else
			return std::pair<double, double>(s, i);
	}
#else
	mpfr_t y;
	mpfr_init2(y, 53); /* Assume IEEE-754 */
	mpfr_set_f(y, x.backend().data(), GMP_RNDD);
	double i = mpfr_get_d(y, GMP_RNDD); /* EXACT but can overflow */
	mpfr_set_f(y, x.backend().data(), GMP_RNDU);
	double s = mpfr_get_d(y, GMP_RNDU); /* EXACT but can overflow */
	mpfr_clear(y);
	return std::pair<double, double>(i, s);
#endif
}

Pair BoostMpToInterval::operator()(const BoostMpzInt &x)
{
#if MPFR_VERSION_MAJOR >= 3
	MPFR_DECL_INIT(y, 53); /* Assume IEEE-754 */
	int    r = mpfr_set_z(y, x.backend().data(), MPFR_RNDA);
	double i = mpfr_get_d(y, MPFR_RNDA); /* EXACT but can overflow */
	if (r == 0 && std::isfinite(i))
		return std::pair<double, double>(i, i);
	else
	{
		double s = std::nextafter(i, 0);
		if (i < 0)
			return std::pair<double, double>(i, s);
		else
			return std::pair<double, double>(s, i);
	}
#else
	mpfr_t y;
	mpfr_init2(y, 53); /* Assume IEEE-754 */
	mpfr_set_z(y, x.backend().data(), GMP_RNDD);
	double i = mpfr_get_d(y, GMP_RNDD); /* EXACT but can overflow */
	mpfr_set_z(y, x.backend().data(), GMP_RNDU);
	double s = mpfr_get_d(y, GMP_RNDU); /* EXACT but can overflow */
	mpfr_clear(y);
	return std::pair<double, double>(i, s);
#endif
}

} // namespace OMC