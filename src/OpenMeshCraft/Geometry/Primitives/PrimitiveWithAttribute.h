#pragma once

#include <utility>

namespace OMC {

/**
 * @brief PrimitiveWithAttribute is a class containing any primitive type
 * and can be attached with any attribute.
 * It can be used in BVH and other similar places.
 * @tparam PrimitiveT Type of primitive.
 * @tparam AttributeT Type of attribute.
 * @note Carefully design names to avoid overwritting.
 */
template <typename PrimitiveT, typename AttributeT>
class PrimitiveWithAttribute : public PrimitiveT
{
public:
	using PT    = PrimitiveT;
	using AT    = AttributeT;
	using ThisT = PrimitiveWithAttribute<PT, AT>;

public:
	/**
	 * @brief Construct Primitive_WA by default.
	 */
	PrimitiveWithAttribute() {}

	/**
	 * @brief Construct Primitive_WA from given primitive.
	 */
	PrimitiveWithAttribute(const PT &primitive)
	  : PrimitiveT(primitive)
	{
	}

	/**
	 * @brief Construct Primitive_WA from given primitive and attribute.
	 */
	PrimitiveWithAttribute(const PT &primitive, const AttributeT &attribute)
	  : PrimitiveT(primitive)
	  , m_attribute(attribute)
	{
	}

	/**
	 * @brief Copy constructor.
	 */
	PrimitiveWithAttribute(const ThisT &rhs)
	  : PrimitiveT(rhs.primitive())
	  , m_attribute(rhs.m_attribute)
	{
	}
	ThisT &operator=(const ThisT &rhs)
	{
		primitive() = rhs.primitive();
		m_attribute = rhs.m_attribute;
		return *this;
	}

	/**
	 * @brief Move constructor.
	 */
	PrimitiveWithAttribute(ThisT &&rhs)
	  : PrimitiveT(std::move(std::move(rhs.primitive())))
	  , m_attribute(std::move(rhs.m_attribute))
	{
	}
	ThisT &operator=(ThisT &&rhs)
	{
		primitive() = std::move(rhs.primitive());
		m_attribute = std::move(rhs.m_attribute);
		return *this;
	}

	/// @name Data access
	/// @{
	PT       &primitive() { return *reinterpret_cast<PT *>(this); }
	const PT &primitive() const { return *reinterpret_cast<const PT *>(this); }
	AT       &attribute() { return m_attribute; }
	const AT &attribute() const { return m_attribute; }
	/// @}

private:
	AT m_attribute;
};

} // namespace OMC