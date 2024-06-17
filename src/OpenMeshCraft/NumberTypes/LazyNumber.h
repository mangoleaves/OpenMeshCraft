#pragma once

#include "BoostMultiprecision.h"
#include "IntervalNumber.h"

#include <atomic>
#include <memory>
#include <mutex>
#include <optional>
#include <thread>
#include <tuple>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

namespace OMC {

/**
 * @brief
 * @tparam HasApprox Set to true to enable an approximate floating-point number.
 * @tparam Protected True if rounding mode of interval number is protected
 * outside LazyNumber.
 * @tparam ThreadSafe When ThreadSafe is true, we guarantee that writing /
 * updating is thread-safe. Otherwise we don't guarantee thread-safe and have
 * higher efficiency. So, disable ThreadSafe when there is definitely no data
 * race.
 */
template <typename HasApprox, typename Protected, typename ThreadSafe>
class LazyExpressionBase
{
public:
	using AT  = IntervalNumber<Protected>;
	using ET  = BoostRational;
	using E2A = BoostMpToInterval;

public:
	LazyExpressionBase() = default;

	LazyExpressionBase(const LazyExpressionBase &)            = delete;
	LazyExpressionBase &operator=(const LazyExpressionBase &) = delete;

	virtual ~LazyExpressionBase() noexcept {}
};

/// @brief Has no approx type version
template <typename Protected, typename ThreadSafe>
class LazyExpressionBase<std::false_type, Protected, ThreadSafe>
{
public:
	using AT  = IntervalNumber<Protected>;
	using ET  = BoostRational;
	using E2A = BoostMpToInterval;

public:
	LazyExpressionBase() = default;

	LazyExpressionBase(const LazyExpressionBase &)            = delete;
	LazyExpressionBase &operator=(const LazyExpressionBase &) = delete;

	template <typename AT_RR>
	LazyExpressionBase(OMC_UNUSED AT_RR &&_at)
	{ // do nothing
	}
};

/// @brief Has approx type version
template <typename Protected, typename ThreadSafe>
class LazyExpressionBase<std::true_type, Protected, ThreadSafe>
{
public:
	using AT  = IntervalNumber<Protected>;
	using ET  = BoostRational;
	using E2A = BoostMpToInterval;

public:
	LazyExpressionBase() = default;

	LazyExpressionBase(const LazyExpressionBase &)            = delete;
	LazyExpressionBase &operator=(const LazyExpressionBase &) = delete;

	template <typename AT_RR>
	LazyExpressionBase(AT_RR &&_at)
	{
		at = std::forward<AT_RR>(_at);
	}

public:
	const AT &approx() const { return at; }

	void set_approx(const AT &_at) const { at = _at; }
	void set_approx(AT &&_at) const { at = std::move(_at); }
	void set_approx(ET *_pet) const { at = AT(E2A()(*_pet)); }

private:
	mutable AT at;
};

template <typename HasApprox, typename Protected, typename ThreadSafe>
class LazyExpression
  : public LazyExpressionBase<HasApprox, Protected, ThreadSafe>
{
public:
	LazyExpression() = default;
	virtual ~LazyExpression() noexcept {}
};

template <typename HasApprox, typename Protected>
class LazyExpression<HasApprox, Protected, /*ThreadSafe*/ std::false_type>
  : public LazyExpressionBase<HasApprox, Protected, std::false_type>
{
public:
	using Base = LazyExpressionBase<HasApprox, Protected, std::false_type>;
	using AT   = typename Base::AT;
	using ET   = typename Base::ET;

public: /* constructions and destructions */
	/// @brief Initialize the lazy expression with an approx number.
	template <typename AT_RR>
	LazyExpression(AT_RR &&_at)
	  : Base(std::forward<AT_RR>(_at))
	{
	}

	/// @brief Initialize with an approx number and an exact number
	template <typename AT_RR, typename ET_RR>
	LazyExpression(AT_RR &&_at, ET_RR &&_et)
	  : Base(std::forward<AT_RR>(_at))
	  , pet(new ET(std::forward<ET_RR>(_et)))
	{
	}

	virtual ~LazyExpression() noexcept {}

public: /* Data access */
	const ET &exact() const
	{
		if (still_lazy())
			this->update_exact();
		return *pet;
	}

	void set_exact(ET *_pet) const { pet = _pet; }

	bool still_lazy() const { return pet == nullptr; }

public: /* Data update */
	virtual void update_exact() const = 0;

private:
	mutable std::unique_ptr<ET> pet;
};

template <typename HasApprox, typename Protected>
class LazyExpression<HasApprox, Protected, /*ThreadSafe*/ std::true_type>
  : public LazyExpressionBase<HasApprox, Protected, std::true_type>
{
public:
	using Base = LazyExpressionBase<HasApprox, Protected, std::true_type>;
	using AT   = typename Base::AT;
	using ET   = typename Base::ET;

public: /* constructions and destructions */
	/// @brief Initialize the lazy expression with an approx number.
	template <typename AT_RR>
	LazyExpression(AT_RR &&_at)
	  : Base(std::forward<AT_RR>(_at))
	{
	}

	/// @brief Initialize with an approx number and an exact number
	template <typename AT_RR, typename ET_RR>
	LazyExpression(AT_RR &&_at, ET_RR &&_et)
	  : Base(std::forward<AT_RR>(_at))
	  , pet(new ET(std::forward<ET_RR>(_et)))
	{
	}

	virtual ~LazyExpression() noexcept { release_pet_if_exists(); }

public: /* Data access */
	const ET &exact() const
	{
		std::call_once(once, [this]() { this->update_exact(); });
		return *pet.load(std::memory_order_relaxed);
	}

	void set_exact(ET *_pet) const { pet.store(_pet, std::memory_order_release); }

	bool still_lazy() const
	{
		return pet.load(std::memory_order_relaxed) == nullptr;
	}

public: /* Data update */
	virtual void update_exact() const = 0;

private:
	/// @brief Delete the pointer to exact number if it exists.
	void release_pet_if_exists()
	{
		ET *p = pet.load(std::memory_order_relaxed);
		if (p != nullptr)
		{
			std::atomic_thread_fence(std::memory_order_acquire);
			delete p;
		}
	}

private:
	mutable std::atomic<ET *> pet = nullptr;
	mutable std::once_flag    once;
};

template <typename CT, typename HasApprox, typename Protected,
          typename ThreadSafe>
class LazyConstant final
  : public LazyExpression<HasApprox, Protected, ThreadSafe>
{
public:
	using LazyExpr       = LazyExpression<HasApprox, Protected, ThreadSafe>;
	using SharedLazyExpr = std::shared_ptr<LazyExpr>;
	using AT             = typename LazyExpr::AT;
	using ET             = typename LazyExpr::ET;

public:
	template <typename T>
	LazyConstant(T &&_c)
	  : LazyExpr(std::forward<T>(_c))
	  , c(std::forward<T>(_c))
	{
	}

	virtual ~LazyConstant() noexcept {}

	virtual void update_exact() const { this->set_exact(new ET(c)); }

private:
	CT c;
};

template <typename HasApprox, typename Protected, typename ThreadSafe>
class LazyConstant<BoostRational, HasApprox, Protected, ThreadSafe> final
  : public LazyExpression<HasApprox, Protected, ThreadSafe>
{
public:
	using LazyExpr       = LazyExpression<HasApprox, Protected, ThreadSafe>;
	using SharedLazyExpr = std::shared_ptr<LazyExpr>;
	using AT             = typename LazyExpr::AT;
	using ET             = typename LazyExpr::ET;

public:
	LazyConstant(const BoostRational &_n)
	  : LazyExpr(HasApprox::value ? BoostMpToInterval()(_n) : AT(), _n)
	{
	}

	LazyConstant(BoostRational &&_n)
	  : LazyExpr(HasApprox::value ? BoostMpToInterval()(_n) : AT(), std::move(_n))
	{
	}

	virtual ~LazyConstant() noexcept {}

	virtual void update_exact() const { /*do nothing*/ }
};

template <template <typename> class UnaryOp, typename HasApprox,
          typename Protected, typename ThreadSafe>
class LazyUnaryExpression final
  : public LazyExpression<HasApprox, Protected, ThreadSafe>
{
public:
	using LazyExpr       = LazyExpression<HasApprox, Protected, ThreadSafe>;
	using SharedLazyExpr = std::shared_ptr<LazyExpr>;
	using AT             = typename LazyExpr::AT;
	using ET             = typename LazyExpr::ET;

public:
	LazyUnaryExpression(const SharedLazyExpr &_op1)
	  : LazyExpr(calc_approx(_op1))
	  , op1(_op1)
	{
	}

	virtual ~LazyUnaryExpression() noexcept {}

	void prune_dag() const { op1.reset(); }

	void update_exact() const
	{
		ET *p = new ET(UnaryOp<ET>()(op1->exact()));
		if constexpr (HasApprox::value)
		{
			if (!this->approx().is_point())
				this->set_approx(p);
		}
		this->set_exact(p);
		this->prune_dag();
	}

private:
	AT calc_approx(const SharedLazyExpr &_op1) const
	{
		if constexpr (HasApprox::value)
			return UnaryOp<AT>()(_op1->approx());
		else
			return AT();
	}

private:
	mutable SharedLazyExpr op1;
};

template <template <typename> class BinaryOp, typename HasApprox,
          typename Protected, typename ThreadSafe>
class LazyBinaryExpression final
  : public LazyExpression<HasApprox, Protected, ThreadSafe>
{
public:
	using LazyExpr       = LazyExpression<HasApprox, Protected, ThreadSafe>;
	using SharedLazyExpr = std::shared_ptr<LazyExpr>;
	using AT             = typename LazyExpr::AT;
	using ET             = typename LazyExpr::ET;

public:
	LazyBinaryExpression(const SharedLazyExpr &_op1, const SharedLazyExpr &_op2)
	  : LazyExpr(calc_approx(_op1, _op2))
	  , op1(_op1)
	  , op2(_op2)
	{
	}

	virtual ~LazyBinaryExpression() noexcept {}

	void prune_dag() const
	{
		op1.reset();
		op2.reset();
	}

	void update_exact() const
	{
		ET *p = new ET(BinaryOp()(op1->exact(), op2->exact()));
		if constexpr (HasApprox::value)
		{
			if (!this->approx().is_point())
				this->set_approx(p);
		}
		this->set_exact(p);
		this->prune_dag();
	}

private:
	AT calc_approx(const SharedLazyExpr &_op1, const SharedLazyExpr &_op2) const
	{
		if constexpr (HasApprox::value)
			return BinaryOp<AT>()(_op1->approx(), _op2->approx());
		else
			return AT();
	}

private:
	mutable SharedLazyExpr op1, op2;
};

/**
 * @brief Lazy evaluated exact number type.
 * @tparam HasApprox Set to true to enable an approximate floating-point number.
 * @tparam Protected Protect the round mode when doing interval number
 * computation.
 * @tparam ThreadSafe Guarantee evaluating the number is thread-safe, otherwise
 * it isn't.
 */
template <typename HasApprox, typename Protected, typename ThreadSafe>
class LazyNumber
{
public:
	using LazyExpr       = LazyExpression<HasApprox, Protected, ThreadSafe>;
	using SharedLazyExpr = std::shared_ptr<LazyExpr>;
	using AT             = typename LazyExpr::AT;
	using ET             = typename LazyExpr::ET;

	using Protector = FPU_RoundingProtector<std::false_type>;

	// clang-format off
	/***** Unary operations for lazy expression *****/
	template<typename CT>
	using LazyConst = LazyConstant<CT, HasApprox, Protected, ThreadSafe>;

	/***** Unary operations for lazy expression *****/
	using LazyUnaryAbs = LazyUnaryExpression<OMC::absolute, HasApprox, Protected, ThreadSafe>;
	using LazyUnaryNeg = LazyUnaryExpression<std::negate, HasApprox, Protected, ThreadSafe>;

	/***** Binary operations for lazy expression *****/
	using LazyBinaryAdd = LazyBinaryExpression<std::plus, HasApprox, Protected, ThreadSafe>;
	using LazyBinarySub = LazyBinaryExpression<std::minus, HasApprox, Protected, ThreadSafe>;
	using LazyBinaryMul = LazyBinaryExpression<std::multiplies, HasApprox, Protected, ThreadSafe>;
	using LazyBinaryDiv = LazyBinaryExpression<std::divides, HasApprox, Protected, ThreadSafe>;
	// clang-format on

public: /* Constructors ******************************************************/
	LazyNumber() = default;

	template <typename T, typename = std::enable_if_t<(std::is_arithmetic_v<T> ||
	                                                   std::is_enum_v<T>) &&
	                                                  !std::is_same_v<T, ET>>>
	LazyNumber(T _n)
	  : expr(std::make_shared<LazyConst<T>>(_n))
	{
	}

	LazyNumber(const ET &_e)
	  : expr(std::make_shared<LazyConst<ET>>(_e))
	{
	}

	LazyNumber(ET &&_e)
	  : expr(std::make_shared<LazyConst<ET>>(std::move(_e)))
	{
	}

	LazyNumber(LazyExpr *_pexpr)
	  : expr(_pexpr)
	{
	}

	LazyNumber(SharedLazyExpr _expr)
	  : expr(_expr)
	{
	}

	LazyNumber(const LazyNumber &_n)
	  : expr(_n.expr)
	{
	}

	LazyNumber(LazyNumber &&_n)
	  : expr(std::move(_n.expr))
	{
	}

	LazyNumber &operator=(const LazyNumber &_n)
	{
		expr = _n.expr;
		return *this;
	}

	LazyNumber &operator=(LazyNumber &&_n)
	{
		expr = std::move(_n.expr);
		return *this;
	}

	void reset() { expr.reset(); }

public: /* Data access ******************************************************/
	const AT &approx() const
	{
		static_assert(HasApprox::value);
		return expr->approx();
	}
	const ET             &exact() const { return expr->exact(); }
	const SharedLazyExpr &expression() const { return expr; }

public: /* Unary operators **************************************************/
	LazyNumber operator+() const { return LazyNumber(expr); }
	LazyNumber operator-() const
	{
		return SharedLazyExpr(new LazyUnaryNeg(expr));
	}

public: /* Binary operators *************************************************/
	friend LazyNumber operator+(const LazyNumber &lhs, const LazyNumber &rhs)
	{
		return SharedLazyExpr(new LazyBinaryAdd(lhs.expr, rhs.expr));
	}
	friend LazyNumber operator-(const LazyNumber &lhs, const LazyNumber &rhs)
	{
		return SharedLazyExpr(new LazyBinarySub(lhs.expr, rhs.expr));
	}
	friend LazyNumber operator*(const LazyNumber &lhs, const LazyNumber &rhs)
	{
		return SharedLazyExpr(new LazyBinaryMul(lhs.expr, rhs.expr));
	}
	friend LazyNumber operator/(const LazyNumber &lhs, const LazyNumber &rhs)
	{
		return SharedLazyExpr(new LazyBinaryDiv(lhs.expr, rhs.expr));
	}

	LazyNumber &operator+=(const LazyNumber &rhs)
	{
		expr = SharedLazyExpr(new LazyBinaryAdd(expr, rhs.expr));
		return *this;
	}
	LazyNumber &operator-=(const LazyNumber &rhs)
	{
		expr = SharedLazyExpr(new LazyBinarySub(expr, rhs.expr));
		return *this;
	}
	LazyNumber &operator*=(const LazyNumber &rhs)
	{
		expr = SharedLazyExpr(new LazyBinaryMul(expr, rhs.expr));
		return *this;
	}
	LazyNumber &operator/=(const LazyNumber &rhs)
	{
		expr = SharedLazyExpr(new LazyBinaryDiv(expr, rhs.expr));
		return *this;
	}

public: /* Comprators ********************************************************/
	friend bool operator<(const LazyNumber &lhs, const LazyNumber &rhs)
	{
		if constexpr (HasApprox::value)
		{
			Certainty res = lhs.approx() < rhs.approx();
			if (is_certain(res))
				return get_certain(res);
		}
		return lhs.exact() < rhs.exact();
	}
	friend bool operator<=(const LazyNumber &lhs, const LazyNumber &rhs)
	{
		if constexpr (HasApprox::value)
		{
			Certainty res = lhs.approx() <= rhs.approx();
			if (is_certain(res))
				return get_certain(res);
		}
		return lhs.exact() <= rhs.exact();
	}
	friend bool operator==(const LazyNumber &lhs, const LazyNumber &rhs)
	{
		if constexpr (HasApprox::value)
		{
			Certainty res = lhs.approx() == rhs.approx();
			if (is_certain(res))
				return get_certain(res);
		}
		return lhs.exact() == rhs.exact();
	}
	friend bool operator!=(const LazyNumber &lhs, const LazyNumber &rhs)
	{
		if constexpr (HasApprox::value)
		{
			Certainty res = lhs.approx() != rhs.approx();
			if (is_certain(res))
				return get_certain(res);
		}
		return lhs.exact() != rhs.exact();
	}
	friend bool operator>=(const LazyNumber &lhs, const LazyNumber &rhs)
	{
		if constexpr (HasApprox::value)
		{
			Certainty res = lhs.approx() >= rhs.approx();
			if (is_certain(res))
				return get_certain(res);
		}
		return lhs.exact() >= rhs.exact();
	}
	friend bool operator>(const LazyNumber &lhs, const LazyNumber &rhs)
	{
		if constexpr (HasApprox::value)
		{
			Certainty res = lhs.approx() > rhs.approx();
			if (is_certain(res))
				return get_certain(res);
		}
		return lhs.exact() > rhs.exact();
	}

	template <typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
	friend bool operator<(const LazyNumber &lhs, T rhs)
	{
		if constexpr (HasApprox::value)
		{
			Certainty res = lhs.approx() < rhs;
			if (is_certain(res))
				return get_certain(res);
		}
		return lhs.exact() < rhs;
	}
	template <typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
	friend bool operator<=(const LazyNumber &lhs, T rhs)
	{
		if constexpr (HasApprox::value)
		{
			Certainty res = lhs.approx() <= rhs;
			if (is_certain(res))
				return get_certain(res);
		}
		return lhs.exact() <= rhs;
	}
	template <typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
	friend bool operator==(const LazyNumber &lhs, T rhs)
	{
		if constexpr (HasApprox::value)
		{
			Certainty res = lhs.approx() == rhs;
			if (is_certain(res))
				return get_certain(res);
		}
		return lhs.exact() == rhs;
	}
	template <typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
	friend bool operator!=(const LazyNumber &lhs, T rhs)
	{
		if constexpr (HasApprox::value)
		{
			Certainty res = lhs.approx() != rhs;
			if (is_certain(res))
				return get_certain(res);
		}
		return lhs.exact() != rhs;
	}
	template <typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
	friend bool operator>=(const LazyNumber &lhs, T rhs)
	{
		if constexpr (HasApprox::value)
		{
			Certainty res = lhs.approx() >= rhs;
			if (is_certain(res))
				return get_certain(res);
		}
		return lhs.exact() >= rhs;
	}
	template <typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
	friend bool operator>(const LazyNumber &lhs, T rhs)
	{
		if constexpr (HasApprox::value)
		{
			Certainty res = lhs.approx() > rhs;
			if (is_certain(res))
				return get_certain(res);
		}
		return lhs.exact() > rhs;
	}

	template <typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
	friend bool operator<(T lhs, const LazyNumber &rhs)
	{
		return rhs > lhs;
	}
	template <typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
	friend bool operator<=(T lhs, const LazyNumber &rhs)
	{
		return rhs >= lhs;
	}
	template <typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
	friend bool operator==(T lhs, const LazyNumber &rhs)
	{
		return rhs == lhs;
	}
	template <typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
	friend bool operator!=(T lhs, const LazyNumber &rhs)
	{
		return rhs != lhs;
	}
	template <typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
	friend bool operator>=(T lhs, const LazyNumber &rhs)
	{
		return rhs <= lhs;
	}
	template <typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
	friend bool operator>(T lhs, const LazyNumber &rhs)
	{
		return rhs < lhs;
	}

private:
	SharedLazyExpr expr;
};

template <typename HasApprox, typename Protected, typename ThreadSafe>
class UnaryOperators<LazyNumber<HasApprox, Protected, ThreadSafe>>
{
public:
	using NT = LazyNumber<HasApprox, Protected, ThreadSafe>;

public:
	static Sign sign(const NT &x)
	{
		if (x > 0)
			return Sign::POSITIVE;
		else if (x < 0)
			return Sign::NEGATIVE;
		else
			return Sign::ZERO;
	}

	static NT abs(const NT &x)
	{
		return NT(std::make_shared<typename NT::LazyUnaryAbs>(x.expression()));
	}

	static NT negate(const NT &x) { return -x; }
};

template <typename HasApprox, typename Protected, typename ThreadSafe>
inline double to_double(const LazyNumber<HasApprox, Protected, ThreadSafe> &n)
{
	return OMC::to_double(n.exact());
}

} // namespace OMC