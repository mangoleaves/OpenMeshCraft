#pragma once

#include "Exception.h"

#include <cstddef>
#include <cstdlib>
#include <memory>
#include <type_traits>

namespace OMC {
/// @brief std::vector will initialize all elements when resize() is called.
/// It waste several seconds when I will parallely initialize all elements after
/// resizing. So, I write this simple C style vector to avoid such case.
/// @note Use this C style vector to store simple types.
/// I don't know what will happen for complex types. :(
template <typename T>
class CStyleVector
{
	static_assert(std::is_trivially_copyable<T>::value);

public: /* Constructor and Destructor *****************************************/
	CStyleVector()
	  : m_data(nullptr)
	  , m_size(0)
	  , m_capacity(0)
	{
	}

	CStyleVector(const CStyleVector &) = delete;

	CStyleVector(CStyleVector &&other)
	{
		m_data           = other.m_data;
		m_size           = other.m_size;
		m_capacity       = other.m_capacity;
		other.m_data     = nullptr;
		other.m_size     = 0;
		other.m_capacity = 0;
	}

	~CStyleVector()
	{
		if (m_data)
			free(m_data);
	}

public: /* Acess and Modify **************************************************/
	const T &operator[](size_t idx) const { return m_data[idx]; }

	T &operator[](size_t idx) { return m_data[idx]; }

	T       *data() { return m_data; }
	const T *data() const { return m_data; }

	T       &front() { return m_data[0]; }
	const T &front() const { return m_data[0]; }

	T       &back() { return m_data[m_size - 1]; }
	const T &back() const { return m_data[m_size - 1]; }

	size_t size() const { return m_size; }

	size_t capacity() const { return m_capacity; }

	bool empty() const { return m_size == 0; }

	void reserve(size_t new_capacity)
	{
		if (new_capacity <= m_capacity)
			return;

		// out of memory exception may be throwed by malloc, just let it throw.
		T *tmp_data = (T *)malloc(new_capacity * sizeof(T));
		// it's better to use copy constructor of T.
		// but we only use this CStyleVector for simple types, just memcpy them.
		if (m_size != 0) // always keep data
			memcpy(tmp_data, m_data, sizeof(T) * m_size);
		// m_data may not be nullptr even if m_size != 0.
		if (m_data)
			free(m_data);
		// set m_data and m_size
		m_data     = tmp_data;
		m_capacity = new_capacity;
		// m_size is unchanged.
	}

	void resize(size_t new_size, bool keep_data = true)
	{
		if (new_size == m_size)
			return;
		else if (new_size < m_capacity)
			m_size = new_size;
		else // insufficient capacity
		{
			// out of memory exception may be throwed by malloc, just let it throw.
			T *tmp_data = (T *)malloc(new_size * sizeof(T));
			// it's better to use copy constructor of T.
			// but we only use this CStyleVector for simple types, just memcpy them.
			if (keep_data && m_size != 0)
				memcpy(tmp_data, m_data, sizeof(T) * m_size);
			// m_data may not be nullptr even if m_size != 0.
			if (m_data)
				free(m_data);
			// set m_data and m_size
			m_data     = tmp_data;
			m_size     = new_size;
			m_capacity = new_size;
		}
	}

	void clear() { m_size = 0; }

	void shrink_to_fit()
	{
		if (m_size == m_capacity)
			return;
	}

public: /* Iterator **********************************************************/
	      // These iterators are very simple and unsafe, use them carefully.
	      // clang-format off
	class iterator
	{
	public:
		using iterator_category = std::random_access_iterator_tag;
		using difference_type   = std::ptrdiff_t;
		using value_type        = T;
		using pointer           = T*;
		using reference         = T&;
	public:
		iterator(pointer _cur_elem)
		  : cur_elem(_cur_elem) {}
		iterator &operator++() { cur_elem++; return *this; }
		iterator &operator--() { cur_elem--; return *this; }
		iterator operator++(int) { iterator cpy(*this); cur_elem++; return cpy; }
		iterator operator--(int) { iterator cpy(*this); cur_elem--; return cpy; }
		iterator &operator+=(size_t step) { cur_elem += step; return *this; }
		iterator operator+(size_t step) const { iterator cpy(*this); cpy.cur_elem += step; return cpy; }
		iterator &operator-=(size_t step) { cur_elem -= step; return *this; }
		iterator operator-(size_t step) const { iterator cpy(*this); cpy.cur_elem -= step; return cpy; }
		size_t operator-(const iterator& other) const { return cur_elem - other.cur_elem; }
		T &operator*() const { return *cur_elem; }
		T &operator*() { return *cur_elem; }
		T *operator->() const { return cur_elem; }
		T *operator->() { return cur_elem; }
		bool operator==(const iterator &rhs) const { return cur_elem == rhs.cur_elem; }
		bool operator!=(const iterator &rhs) const { return cur_elem != rhs.cur_elem; }
	public:
		T* cur_elem;
	};

	iterator begin() {return iterator(m_data);}
	iterator end() {return iterator(m_data + m_size);}

	class const_iterator
	{
	public:
		using iterator_category = std::random_access_iterator_tag;
		using difference_type   = std::ptrdiff_t;
		using value_type        = T;
		using pointer           = const T*;
		using reference         = const T&;
	public:
		const_iterator(pointer _cur_elem)
		  : cur_elem(_cur_elem) {}
		const_iterator &operator++() { cur_elem++; return *this; }
		const_iterator &operator--() { cur_elem--; return *this; }
		const_iterator operator++(int) { const_iterator cpy(*this); cur_elem++; return cpy; }
		const_iterator operator--(int) { const_iterator cpy(*this); cur_elem--; return cpy; }
		const_iterator &operator+=(size_t step) { cur_elem += step; return *this; }
		const_iterator operator+(size_t step) { const_iterator cpy(*this); cpy.cur_elem += step; return cpy; }
		const_iterator &operator-=(size_t step) { cur_elem -= step; return *this; }
		const_iterator operator-(size_t step) { const_iterator cpy(*this); cpy.cur_elem -= step; return cpy; }
		size_t operator-(const const_iterator& other) const { return cur_elem - other.cur_elem; }
		T &operator*() const { return *cur_elem; }
		T &operator*() { return *cur_elem; }
		T *operator->() const { return cur_elem; }
		T *operator->() { return cur_elem; }
		bool operator==(const const_iterator &rhs) const { return cur_elem == rhs.cur_elem; }
		bool operator!=(const const_iterator &rhs) const { return cur_elem != rhs.cur_elem; }
	public:
		T* cur_elem;
	};
	
	const_iterator cbegin() const {return const_iterator(m_data);}
	const_iterator cend() const {return const_iterator(m_data + m_size);}
	const_iterator begin() const {return const_iterator(m_data);}
	const_iterator end() const {return const_iterator(m_data + m_size);}
	      // clang-format on
protected:
	T     *m_data;
	size_t m_size;
	size_t m_capacity;
};
} // namespace OMC