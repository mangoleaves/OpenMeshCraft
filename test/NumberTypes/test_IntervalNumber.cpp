#include "OpenMeshCraft/NumberTypes/IntervalNumber.h"

#include "test_utils.h"

#if defined(OMC_SSE2) && !defined(__SSE2__)
#define __SSE2__
#endif
#include "CGAL/Interval_nt.h"

#include <type_traits>

using namespace OMC;

TEST(IntervalNumber, Construct)
{
	IntervalNumber it0(1.);
	IntervalNumber it1(1., 2.);
	// EXPECT_THROW(it1 = IntervalNumber(-1.1, 0.9), std::runtime_error);
}

TEST(IntervalNumber, OperatorsComparedWithCGAL)
{
	auto interval_number_equal =
	  [](const IntervalNumber<std::true_type> &it, const CGAL::Interval_nt<false> &_it)
	{
		if (!(it.inf() == _it.inf() && it.sup() == _it.sup()))
		{
			std::cout.precision(20);
			std::cout << "it " << it.inf() << ", " << it.sup() << std::endl;
			std::cout << "_it " << _it.inf() << ", " << _it.sup() << std::endl;
			return false;
		}
		else
			return true;
	};

	srand(0);
	FPU_set_round(FE_UPWARD);
	for (size_t i = 0; i < 10000000; i++)
	{
		double d1 = rand() * M_PI, d2 = rand() * M_PI;

		IntervalNumber<std::true_type> it1(d1), it2(d2);
		IntervalNumber<std::true_type> it3  = it1 + it2;
		IntervalNumber<std::true_type> it4  = it1 - it2;
		IntervalNumber<std::true_type> it5  = it1 * it2;
		IntervalNumber<std::true_type> it6  = it3 * it4;
		IntervalNumber<std::true_type> it7  = it6 - it5;
		IntervalNumber<std::true_type> it8  = it7 + it6;
		IntervalNumber<std::true_type> it9  = it8 * it3;
		IntervalNumber<std::true_type> it10 = it9 / it3;
		IntervalNumber<std::true_type> it11 = it10 / it7;

		CGAL::Interval_nt<false> _it1(d1), _it2(d2);
		CGAL::Interval_nt<false> _it3  = _it1 + _it2;
		CGAL::Interval_nt<false> _it4  = _it1 - _it2;
		CGAL::Interval_nt<false> _it5  = _it1 * _it2;
		CGAL::Interval_nt<false> _it6  = _it3 * _it4;
		CGAL::Interval_nt<false> _it7  = _it6 - _it5;
		CGAL::Interval_nt<false> _it8  = _it7 + _it6;
		CGAL::Interval_nt<false> _it9  = _it8 * _it3;
		CGAL::Interval_nt<false> _it10 = _it9 / _it3;
		CGAL::Interval_nt<false> _it11 = _it10 / _it7;

		ASSERT_TRUE(interval_number_equal(it1, _it1));
		ASSERT_TRUE(interval_number_equal(it2, _it2));
		ASSERT_TRUE(interval_number_equal(it3, _it3));
		ASSERT_TRUE(interval_number_equal(it4, _it4));
		ASSERT_TRUE(interval_number_equal(it5, _it5));
		ASSERT_TRUE(interval_number_equal(it6, _it6));
		ASSERT_TRUE(interval_number_equal(it7, _it7));
		ASSERT_TRUE(interval_number_equal(it8, _it8));
		ASSERT_TRUE(interval_number_equal(it9, _it9));
		ASSERT_TRUE(interval_number_equal(it10, _it10));
		ASSERT_TRUE(interval_number_equal(it11, _it11));
	}
	FPU_set_round(FE_TONEAREST);
}