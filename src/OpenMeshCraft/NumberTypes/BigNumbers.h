#pragma once

#include "OpenMeshCraft/Utils/Exception.h"
#include "NumberUtils.h"

#include <cstdlib>
#include <iostream>
#include <string>

namespace OMC {

/////////////////////////////////////////////////////////////////////
//
// 	   B I G   N A T U R A L
//
/////////////////////////////////////////////////////////////////////

// A bignatural is an arbitrarily large non-negative integer.
// It is made of a sequence of digits in base 2^32.
// Leading zero-digits are not allowed.
// The value 'zero' is represented by an empty digit sequence.

class bignatural
{
	uint32_t *digits;     // Ptr to the digits
	uint32_t  m_size;     // Actual number of digits
	uint32_t  m_capacity; // Current vector capacity

	inline static uint32_t *BN_ALLOC(uint32_t num_bytes)
	{
		return (uint32_t *)malloc(num_bytes);
	}
	inline static void BN_FREE(uint32_t *ptr) { free(ptr); }

	void init(const bignatural &m);
	void init(const uint64_t m);

public:
	// Creates a 'zero'
	bignatural()
	  : digits(NULL)
	  , m_size(0)
	  , m_capacity(0)
	{
	}

	~bignatural() { BN_FREE(digits); }

	bignatural(const bignatural &m) { init(m); }

	bignatural(bignatural &&m) noexcept
	  : digits(m.digits)
	  , m_size(m.m_size)
	  , m_capacity(m.m_capacity)
	{
		m.digits = nullptr;
	}

	// Creates from uint64_t
	bignatural(uint64_t m) { init(m); }

	// If the number fits a uint64_t convert and return true
	bool toUint64(uint64_t &n) const;

	// If the number fits a uint32_t convert and return true
	bool toUint32(uint32_t &n) const;

	bignatural &operator=(const bignatural &m);

	bignatural &operator=(const uint64_t m);

	inline const uint32_t &back() const { return digits[m_size - 1]; }

	inline const uint32_t &operator[](int i) const { return digits[i]; }

	inline uint32_t size() const { return m_size; }

	inline bool empty() const { return m_size == 0; }

	// Left-shift by n bits and possibly add limbs as necessary
	void operator<<=(uint32_t n);

	// Right-shift by n bits
	void operator>>=(uint32_t n);

	bool operator==(const bignatural &b) const;

	bool operator!=(const bignatural &b) const;

	bool operator>=(const bignatural &b) const;

	bool operator>(const bignatural &b) const;

	bool operator<=(const bignatural &b) const { return b >= *this; }

	bool operator<(const bignatural &b) const { return b > *this; }

	bignatural operator+(const bignatural &b) const;

	// Assume that b is smaller than or equal to this number!
	bignatural operator-(const bignatural &b) const;

	bignatural operator*(const bignatural &b) const;

	// Short division algorithm
	bignatural divide_by(const uint32_t D, uint32_t &remainder) const;

	uint32_t getNumSignificantBits() const;

	bool getBit(uint32_t b) const;

	// Long division
	bignatural divide_by(const bignatural &divisor, bignatural &remainder) const;

	// Greatest common divisor (Euclidean algorithm)
	bignatural GCD(const bignatural &D) const;

	// String representation in decimal form
	std::string get_dec_str() const;

	// String representation in binary form
	std::string get_str() const;

	// Count number of zeroes on the right (least significant binary digits)
	uint32_t countEndingZeroes() const;

public:
	inline uint32_t &back() { return digits[m_size - 1]; }

	inline void pop_back() { m_size--; }

	inline uint32_t &operator[](int i) { return digits[i]; }

	inline void push_back(uint32_t b)
	{
		if (m_size == m_capacity)
			increaseCapacity((m_capacity | 1) << 2);
		digits[m_size++] = b;
	}

	inline void push_bit_back(uint32_t b)
	{
		if (m_size)
		{
			operator<<=(1);
			back() |= b;
		}
		else if (b)
			push_back(1);
	}

	inline void reserve(uint32_t n)
	{
		if (n > m_capacity)
			increaseCapacity(n);
	}

	inline void resize(uint32_t n)
	{
		reserve(n);
		m_size = n;
	}

	inline void fill(uint32_t v)
	{
		for (uint32_t i = 0; i < m_size; i++)
			digits[(int)i] = v;
	}

	inline void pop_front()
	{
		for (uint32_t i = 1; i < m_size; i++)
			digits[i - 1] = digits[i];
		pop_back();
	}

	// Count number of zeroes on the left (most significant digits)
	uint32_t countLeadingZeroes() const;

	// Count number of zeroes on the right of the last 1 in the least significant
	// limb Assumes that number is not zero and last limb is not zero!
	uint32_t countEndingZeroesLSL() const
	{
		const uint32_t d = back();
		uint32_t       i = 31;
		while (!(d << i))
			i--;
		return 31 - i;
	}

	inline void pack()
	{
		uint32_t i = 0;
		while (i < m_size && digits[i] == 0)
			i++;
		if (i)
		{
			uint32_t *dold = digits + i;
			uint32_t *dnew = digits;
			uint32_t *dend = digits + m_size;
			while (dold < dend)
				*dnew++ = *dold++;
			m_size -= i;
		}
	}

	// a and b must NOT be this number!
	void toSum(const bignatural &a, const bignatural &b);

	// a and b must NOT be this number!
	// Assume that b is smaller or equal than a!
	void toDiff(const bignatural &a, const bignatural &b);

	// a and b must NOT be this number!
	void toProd(const bignatural &a, const bignatural &b);

private:
	// Multiplies by a single limb, left shift, and add to accumulator. Does not
	// pack!
	void addmul(uint32_t b, uint32_t left_shifts, bignatural &result) const;

	void increaseCapacity(uint32_t new_capacity);

	friend class bigfloat;
};

inline std::ostream &operator<<(std::ostream &os, const bignatural &p)
{
	os << p.get_dec_str();
	return os;
}

inline void bignatural::init(const bignatural &m)
{
	m_size     = m.m_size;
	m_capacity = m.m_capacity;
	if (m_capacity)
	{
		digits = (uint32_t *)BN_ALLOC(sizeof(uint32_t) * m_capacity);
		memcpy(digits, m.digits, sizeof(uint32_t) * m_size);
	}
	else
		digits = NULL;
}

inline void bignatural::init(const uint64_t m)
{
	if (m == 0)
	{
		m_size = m_capacity = 0;
		digits              = NULL;
	}
	else if (m <= UINT32_MAX)
	{
		m_size = m_capacity = 1;
		digits              = (uint32_t *)BN_ALLOC(sizeof(uint32_t));
		digits[0]           = (uint32_t)m;
	}
	else
	{
		m_size = m_capacity = 2;
		digits              = (uint32_t *)BN_ALLOC(sizeof(uint32_t) * 2);
		digits[0]           = (uint32_t)(m >> 32);
		digits[1]           = (uint32_t)(m & UINT32_MAX);
	}
}

inline bool bignatural::toUint64(uint64_t &n) const
{
	if (m_size == 0)
		n = 0;
	else if (m_size == 1)
		n = digits[0];
	else if (m_size == 2)
	{
		n = (((uint64_t)digits[0]) << 32) + digits[1];
	}
	else
		return false;

	return true;
}

inline bool bignatural::toUint32(uint32_t &n) const
{
	if (m_size == 0)
		n = 0;
	else if (m_size == 1)
		n = digits[0];
	else
		return false;

	return true;
}

inline bignatural &bignatural::operator=(const bignatural &m)
{
	if (digits != m.digits)
	{
		BN_FREE(digits);
		init(m);
	}

	return *this;
}

inline bignatural &bignatural::operator=(const uint64_t m)
{
	BN_FREE(digits);
	init(m);
	return *this;
}

inline void bignatural::operator<<=(uint32_t n)
{
	uint32_t s  = n & 0x0000001f;
	uint32_t lz = countLeadingZeroes();
	if (lz < s)
	{ // Need a further limb
		push_back(0);
		s = 32 - s;
		for (int i = (int)m_size - 1; i > 0; i--)
		{
			digits[i] >>= s;
			digits[i] |= (digits[i - 1] << (32 - s));
		}
		digits[0] >>= s;
	}
	else if (s)
	{ // Leading zeroes are enough
		for (int i = 0; i < (int)m_size - 1; i++)
		{
			digits[i] <<= s;
			digits[i] |= (digits[i + 1] >> (32 - s));
		}
		back() <<= s;
	}

	while (n >= 32)
	{
		push_back(0);
		n -= 32;
	}
}

inline void bignatural::operator>>=(uint32_t n)
{
	while (n >= 32)
	{
		pop_back();
		n -= 32;
	}
	if (!n)
		return;

	for (uint32_t i = m_size; i > 1; i--)
	{
		digits[i - 1] >>= n;
		digits[i - 1] |= (digits[i - 2] << (32 - n));
	}
	digits[0] >>= n;
	if (digits[0] == 0)
		pop_front();
}

inline bool bignatural::operator==(const bignatural &b) const
{
	if (size() != b.size())
		return false;
	for (uint32_t i = 0; i < size(); i++)
		if (digits[i] == b.digits[i])
			return false;
	return true;
}

inline bool bignatural::operator!=(const bignatural &b) const
{
	if (size() != b.size())
		return true;
	for (uint32_t i = 0; i < size(); i++)
		if (digits[i] == b.digits[i])
			return true;
	return false;
}

inline bool bignatural::operator>=(const bignatural &b) const
{
	const int s = (size() > b.size()) - (size() < b.size());
	if (s)
		return (s > 0);
	uint32_t i;
	for (i = 0; i < size() && digits[i] == b.digits[i]; i++)
		;
	return (i == size() || digits[i] > b.digits[i]);
}

inline bool bignatural::operator>(const bignatural &b) const
{
	const int s = (size() > b.size()) - (size() < b.size());
	if (s)
		return (s > 0);
	uint32_t i;
	for (i = 0; i < size() && digits[i] == b.digits[i]; i++)
		;
	return (i != size() && digits[i] > b.digits[i]);
}

inline bignatural bignatural::operator+(const bignatural &b) const
{
	bignatural result;
	result.toSum(*this, b);
	return result;
}

// Assume that b is smaller than or equal to this number!
inline bignatural bignatural::operator-(const bignatural &b) const
{
	bignatural result;
	result.toDiff(*this, b);
	return result;
}

inline bignatural bignatural::operator*(const bignatural &b) const
{
	bignatural result;
	result.toProd(*this, b);
	return result;
}

// Short division algorithm
inline bignatural bignatural::divide_by(const uint32_t D,
                                        uint32_t      &remainder) const
{
	/*if (D == 0)
	  OMC_ASSERT(false, "Division by zero\n");*/
	if (m_size == 0)
		return 0;

	// If both dividend fits into 64 bits, use hardware division
	uint64_t n;
	if (toUint64(n))
	{
		remainder = n % D;
		return n / D;
	}

	bignatural Q;
	uint32_t   next_digit = 0;
	uint64_t   dividend   = digits[next_digit++];
	for (;;)
	{
		uint64_t tmp_div = dividend / D;
		if (!Q.empty() || tmp_div)
			Q.push_back((uint32_t)tmp_div);
		dividend -= (tmp_div * D);
		if (next_digit < m_size)
		{
			dividend <<= 32;
			dividend += digits[next_digit++];
		}
		else
			break;
	}
	remainder = (uint32_t)dividend;

	return Q;
}

inline uint32_t bignatural::getNumSignificantBits() const
{
	if (!m_size)
		return 0;
	int nsb = 31;
	while (!(digits[0] & (1 << nsb)))
		nsb--;
	nsb++;
	return nsb + (m_size - 1) * 32;
}

inline bool bignatural::getBit(uint32_t b) const
{
	const uint32_t dig = (m_size - (b >> 5)) - 1;
	const uint32_t bit = b & 31;
	return (digits[dig] & (1 << bit));
}

// Long division
inline bignatural bignatural::divide_by(const bignatural &divisor,
                                        bignatural       &remainder) const
{
	/*if (divisor.empty())
	  OMC_ASSERT(false, "Division by zero\n");*/
	if (empty())
		return 0;

	// If divisor fits into 32 bits, revert to short division
	uint32_t d32, rem;
	if (divisor.toUint32(d32))
	{
		bignatural q = divide_by(d32, rem);
		remainder    = rem;
		return q;
	}

	// If both dividend and divisor fit into 64 bits, use hardware division
	uint64_t n, d;
	if (toUint64(n) && divisor.toUint64(d))
	{
		remainder = n % d;
		return n / d;
	}

	// If divisor is greater than dividend...
	if (divisor > *this)
	{
		remainder = *this;
		return 0;
	}

	// Use binary (per-bit) long division
	const bignatural &dividend = *this;

	bignatural quotient, loc_dividend;
	uint32_t   next_dividend_bit = dividend.getNumSignificantBits();

	do
	{
		loc_dividend.push_bit_back(dividend.getBit(--next_dividend_bit));
		if (loc_dividend >= divisor)
		{
			loc_dividend = loc_dividend - divisor;
			quotient.push_bit_back(1);
		}
		else if (!quotient.empty())
			quotient.push_bit_back(0);
	} while (next_dividend_bit);

	remainder = loc_dividend;

	return quotient;
}

// Greatest common divisor (Euclidean algorithm)
inline bignatural bignatural::GCD(const bignatural &D) const
{
	bignatural A = *this;
	bignatural B = D;
	bignatural R;
	while (!A.empty() && !B.empty())
	{
		A.divide_by(B, R);
		A = B;
		B = R;
	}
	if (A.empty())
		return B;
	else
		return A;
}

// String representation in decimal form
inline std::string bignatural::get_dec_str() const
{
	std::string st;
	bignatural  N = *this;
	uint32_t    R;
	if (N.empty())
		return "0";
	while (!N.empty())
	{
		N = N.divide_by(10, R);
		st += (std::string("0") + std::to_string(R));
	}
	std::reverse(st.begin(), st.end());

	return st;
}

// String representation in binary form
inline std::string bignatural::get_str() const
{
	std::string st;
	char        s[33];
	s[32] = 0;
	for (uint32_t j = 0; j < m_size; j++)
	{
		for (int i = 0; i < 32; i++)
			s[i] = (digits[j] & (((uint32_t)1) << (31 - i))) ? '1' : '0';
		st += s;
	}
	return st;
}

// Count number of zeroes on the right (least significant binary digits)
inline uint32_t bignatural::countEndingZeroes() const
{
	if (m_size == 0)
		return 0;
	uint32_t i    = m_size - 1;
	uint32_t shft = 0;
	while (!digits[i])
	{
		i--;
		shft += 32;
	}

	const uint32_t d = digits[i];
	uint32_t       j = 31;
	while (!(d << j))
		j--;
	return shft + 31 - j;

	// uint32_t s = UINT32_MAX;
	// uint32_t m = digits[i];
	// while ((s & m) == m) {
	//	s <<= 1;
	//	shft++;
	// }
	// return shft - 1;
}

inline uint32_t bignatural::countLeadingZeroes() const
{
	uint32_t       s    = UINT32_MAX;
	const uint32_t m    = digits[0];
	uint32_t       shft = 0;
	while ((s & m) == m)
	{
		s >>= 1;
		shft++;
	}
	return shft - 1;
}

inline void bignatural::toSum(const bignatural &a, const bignatural &b)
{
	if (a.m_size == 0)
		operator=(b);
	else if (b.m_size == 0)
		operator=(a);
	else
	{
		const uint32_t a_s   = a.m_size;
		const uint32_t b_s   = b.m_size;
		uint64_t       carry = 0;
		uint32_t      *dig_a = a.digits + a_s;
		uint32_t      *dig_b = b.digits + b_s;

		if (a_s == b_s)
		{
			resize(a_s + 1);
			uint32_t *dig_r = digits + a_s + 1;
			do
			{
				const uint64_t da = *(--dig_a);
				const uint64_t db = *(--dig_b);
				const uint64_t sm = da + db + carry;
				*(--dig_r)        = sm & UINT32_MAX;
				carry             = (sm >> 32);
			} while (dig_a != a.digits);
		}
		else if (a_s < b_s)
		{
			resize(b_s + 1);
			uint32_t *dig_r = digits + b_s + 1;
			do
			{
				const uint64_t da = *(--dig_a);
				const uint64_t db = *(--dig_b);
				const uint64_t sm = da + db + carry;
				*(--dig_r)        = sm & UINT32_MAX;
				carry             = (sm >> 32);
			} while (dig_a != a.digits);
			do
			{
				const uint64_t db = *(--dig_b);
				const uint64_t sm = db + carry;
				*(--dig_r)        = sm & UINT32_MAX;
				carry             = (sm >> 32);
			} while (dig_b != b.digits);
		}
		else
		{
			resize(a_s + 1);
			uint32_t *dig_r = digits + a_s + 1;
			do
			{
				const uint64_t da = *(--dig_a);
				const uint64_t db = *(--dig_b);
				const uint64_t sm = da + db + carry;
				*(--dig_r)        = sm & UINT32_MAX;
				carry             = (sm >> 32);
			} while (dig_b != b.digits);
			do
			{
				const uint64_t da = *(--dig_a);
				const uint64_t sm = da + carry;
				*(--dig_r)        = sm & UINT32_MAX;
				carry             = (sm >> 32);
			} while (dig_a != a.digits);
		}
		if (carry)
			digits[0] = (uint32_t)carry;
		else
		{
			uint32_t *dold = digits + 1;
			uint32_t *dnew = digits;
			uint32_t *dend = digits + m_size;
			while (dold < dend)
				*dnew++ = *dold++;
			m_size--;
		}
	}
}

// a and b must NOT be this number!
// Assume that b is smaller or equal than a!
inline void bignatural::toDiff(const bignatural &a, const bignatural &b)
{
	if (b.m_size == 0)
		operator=(a);
	else
	{
		const uint32_t a_s      = a.m_size;
		const uint32_t b_s      = b.m_size;
		uint32_t       res_size = a_s;
		if (b_s > res_size)
			res_size = b_s;
		resize(res_size);

		uint64_t debt = 0;
		for (uint32_t i = 1; i <= res_size; i++)
		{
			const uint64_t da = (i <= a_s) ? (a.digits[(int)(a_s - i)]) : (0);
			const uint64_t db =
			  ((i <= b_s) ? (b.digits[(int)(b_s - i)]) : (0)) + debt;
			debt = !(da >= db);
			if (debt)
				digits[(int)(res_size - i)] =
				  (uint32_t)((da + (((uint64_t)1) << 32)) - db);
			else
				digits[(int)(res_size - i)] = (uint32_t)(da - db);
		}
		pack();
	}
}

// a and b must NOT be this number!
inline void bignatural::toProd(const bignatural &a, const bignatural &b)
{
	if (a.empty())
		operator=(a);
	else if (b.empty())
		operator=(b);
	else
	{
		resize(a.m_size + b.m_size);
		fill(0);

		uint32_t ls = 0;
		for (uint32_t i = b.m_size; i > 0; i--)
			a.addmul(b[(int)(i - 1)], ls++, *this);

		pack();
	}
}

inline void bignatural::addmul(uint32_t b, uint32_t left_shifts,
                               bignatural &result) const
{
	uint64_t carry = 0;
	int      d     = (int)(result.m_size - m_size - left_shifts);
	for (uint32_t i = m_size; i > 0; i--)
	{
		uint64_t pm =
		  ((uint64_t)digits[(int)(i - 1)]) * b + carry + result[(int)i + d - 1];
		result[(int)i + d - 1] = (uint32_t)pm;
		carry                  = pm >> 32;
	}
	result[d - 1] = (uint32_t)carry;
}

inline void bignatural::increaseCapacity(uint32_t new_capacity)
{
	m_capacity      = new_capacity;
	uint32_t *tmp_d = (uint32_t *)BN_ALLOC(sizeof(uint32_t) * m_capacity);
	memcpy(tmp_d, digits, sizeof(uint32_t) * m_size);
	BN_FREE(digits);
	digits = tmp_d;
}

/////////////////////////////////////////////////////////////////////
//
// 	   B I G   F L O A T
//
/////////////////////////////////////////////////////////////////////

// A bigfloat is a floting point number with arbitrarily large mantissa.
// In principle, we could have made the exponent arbitrarily large too,
// but in practice this appears to be useless.
// Exponents are in the range [-INT32_MAX, INT32_MAX]
//
// A bigfloat f evaluates to f = sign * mantissa * 2^exponent
//
// mantissa is a bignatural whose least significant bit is 1.
// Number is zero if mantissa is empty.

class bigfloat
{
	bignatural mantissa; // .back() is less significant. Use 32-bit limbs to avoid
	                     // overflows using 64-bits
	int32_t exponent; // In principle we might still have under/overflows, but not
	                  // in practice
	int32_t sign;     // Redundant but keeps alignment

public:
	// Default constructor creates a zero-valued bigfloat
	bigfloat()
	  : exponent(0)
	  , sign(0)
	{
	}

	// Lossless conversion from double
	bigfloat(const double d);

	// Truncated approximation
	double get_d() const;

	bigfloat operator+(const bigfloat &b) const;

	bigfloat operator-(const bigfloat &b) const;

	bigfloat operator*(const bigfloat &b) const;

	void invert() { sign = -sign; }

	bigfloat inverse() const
	{
		bigfloat r = *this;
		r.invert();
		return r;
	}

	inline int sgn() const { return sign; }

	std::string get_str() const;

	const bignatural &getMantissa() const { return mantissa; }
	int32_t           getExponent() const { return exponent; }

private:
	// Right-shift as long as the least significant bit is zero
	void pack();

	// Left-shift the mantissa by n bits and reduce the exponent accordingly
	void leftShift(uint32_t n)
	{
		mantissa <<= n;
		exponent -= n;
	}
};

inline bigfloat operator-(const bigfloat &f) { return f.inverse(); }

inline bigfloat operator*(double d, const bigfloat &f)
{
	return f * bigfloat(d);
}

inline bigfloat::bigfloat(const double d)
{
	union
	{
		uint64_t ui;
		double   d;
	} u;
	u.d = d;

	sign = (d > 0) - (d < 0);

	if (sign)
	{
		const uint64_t m = (u.ui & 0x000fffffffffffff) + 0x0010000000000000;
		mantissa.push_back(m >> 32);
		mantissa.push_back(m & 0x00000000ffffffff);
		u.ui <<= 1;
		u.ui >>= 53;
		exponent = ((int32_t)u.ui) - 1075; // Exp

		pack();
	}
	else
		exponent = 0;
}

inline double bigfloat::get_d() const
{
	union
	{
		uint64_t ui;
		double   d;
	} u;
	u.ui = 0;
	if (mantissa.empty())
		return 0.0;

	uint64_t m;
	int32_t  e;
	uint32_t shft;

	if (mantissa.size() == 1)
	{
		m    = ((uint64_t)mantissa[0]);
		shft = mantissa.countLeadingZeroes() + 21;
		m <<= shft;
		e = exponent - shft;
	}
	else
	{
		m    = (((uint64_t)mantissa[0]) << 32) | ((uint64_t)mantissa[1]);
		e    = exponent + 32 * ((uint32_t)mantissa.size() - 2);
		shft = mantissa.countLeadingZeroes();

		if (shft < 11)
		{
			m >>= (11 - shft);
			e += (11 - shft);
		}
		if (shft > 11)
		{
			m <<= (shft - 11);
			e -= (shft - 11);
			if (mantissa.size() > 2)
				m |= (mantissa[2] >> (43 - shft));
		}
	}
	m &= (~0x0010000000000000); // Remove implicit digit
	e += 52;

	if (e < (-1022))
		return 0.0;
	if (e > 1023)
		return sign * INFINITY;

	if (sign < 0)
		u.ui |= 0x8000000000000000;           // Set sign
	u.ui |= (((uint64_t)(e + 1023)) << 52); // Set exponent
	u.ui |= m;                              // Set mantissa

	return u.d;
}

inline bigfloat bigfloat::operator+(const bigfloat &b) const
{
	if (mantissa.empty())
		return b;
	if (b.mantissa.empty())
		return *this;

	if (exponent == b.exponent)
	{
		bigfloat result;

		if (sign == b.sign)
		{
			result.mantissa.toSum(mantissa, b.mantissa);
			result.sign = sign;
		}
		else if (b.mantissa >= mantissa)
		{
			result.mantissa.toDiff(b.mantissa, mantissa);
			result.sign = b.sign;
		}
		else
		{
			result.mantissa.toDiff(mantissa, b.mantissa);
			result.sign = sign;
		}

		result.exponent = exponent;
		result.pack();
		return result;
	}
	else if (exponent > b.exponent)
	{
		bigfloat op(*this);
		op.leftShift(exponent - b.exponent);
		return op + b;
	}
	else
	{ // exponent < b.exponent
		bigfloat op(b);
		op.leftShift(b.exponent - exponent);
		return op + *this;
	}
}

inline bigfloat bigfloat::operator-(const bigfloat &b) const
{
	if (mantissa.empty())
		return b.inverse();
	if (b.mantissa.empty())
		return *this;

	if (exponent == b.exponent)
	{
		bigfloat result;

		if (sign != b.sign)
		{
			result.mantissa.toSum(mantissa, b.mantissa);
			result.sign = sign;
		}
		else if (b.mantissa >= mantissa)
		{
			result.mantissa.toDiff(b.mantissa, mantissa);
			result.sign = -sign;
		}
		else
		{
			result.mantissa.toDiff(mantissa, b.mantissa);
			result.sign = sign;
		}

		result.exponent = exponent;
		result.pack();
		return result;
	}
	else if (exponent > b.exponent)
	{
		bigfloat op(*this);
		op.leftShift(exponent - b.exponent);
		return op - b;
	}
	else
	{ // exponent < b.exponent
		bigfloat op(b);
		op.leftShift(b.exponent - exponent);
		return *this - op;
	}
}

inline bigfloat bigfloat::operator*(const bigfloat &b) const
{
	if (mantissa.empty() || b.mantissa.empty())
		return 0;

	// Left-shift operator with highest exponent
	if (exponent == b.exponent)
	{
		bigfloat result;
		result.mantissa.toProd(mantissa, b.mantissa);
		result.exponent = exponent;
		result.sign     = sign * b.sign;
		result.leftShift(result.exponent - exponent);
		result.exponent *= 2;

		result.pack();
		return result;
	}
	else if (exponent > b.exponent)
	{
		bigfloat op(*this);
		op.leftShift(exponent - b.exponent);
		return op * b;
	} // exponent < b.exponent
	else
	{
		bigfloat op(b);
		op.leftShift(b.exponent - exponent);
		return op * *this;
	}
}

inline std::string bigfloat::get_str() const
{
	std::string s;
	if (sign == 0)
		s += "0";
	if (sign < 0)
		s += "-";
	s += mantissa.get_str();
	s += " * 2^";
	s += std::to_string(exponent);
	return s;
}

inline void bigfloat::pack()
{
	if (mantissa.empty())
	{
		sign = exponent = 0;
		return;
	}

	while (mantissa.back() == 0)
	{
		mantissa.pop_back();
		exponent += 32;
	}

	const uint32_t s = mantissa.countEndingZeroesLSL();
	if (s)
	{
		for (int i = (int)mantissa.size() - 1; i > 0; i--)
		{
			mantissa[i] >>= s;
			mantissa[i] |= (mantissa[i - 1] << (32 - s));
		}
		mantissa[0] >>= s;
		exponent += s;
	}

	mantissa.pack();
}

template <>
class UnaryOperators<bigfloat>
{
public:
	using NT = bigfloat;

public:
	static Sign sign(const NT &x) { return static_cast<Sign>(x.sgn()); }

	static NT abs(const NT &x) { return x.sgn() == -1 ? x.inverse() : x; }

	static NT negate(const NT &x) { return x.inverse(); }
};

class BigNumbersToDouble
{
public:
	double operator()(const bigfloat &x) { return x.get_d(); }
};

template <>
inline double to_double(const bigfloat &x)
{
	return BigNumbersToDouble()(x);
}


} // namespace OMC