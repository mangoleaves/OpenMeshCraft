#include "StringUtils.h"

namespace OMC {

/**
 * @brief Split string with given delimitor.
 * @param str Input string
 * @param delim delimitor
 * @return std::vector<std::string> splitting result.
 */
std::vector<std::string> split_string(const std::string &str, char delim)
{
	std::size_t              previous = 0;
	std::size_t              current  = str.find(delim);
	std::vector<std::string> elems;
	while (current != std::string::npos)
	{
		if (current > previous)
		{
			elems.push_back(str.substr(previous, current - previous));
		}
		previous = current + 1;
		current  = str.find(delim, previous);
	}
	if (previous != str.size())
	{
		elems.push_back(str.substr(previous));
	}
	return elems;
};

/**
 * @brief Trim the spaces before and after string.
 * @param[inout] str The string to trim.
 */
void trim_string(std::string &str)
{
	size_t start = str.find_first_not_of(" \t\r\n");
	size_t end   = str.find_last_not_of(" \t\r\n");

	if ((std::string::npos == start) || (std::string::npos == end))
		str = "";
	else
		str = str.substr(start, end - start + 1);
}

bool starts_with(const std::string &big, const std::string &small)
{
	if (&big == &small)
		return true;
	const typename std::string::size_type big_size   = big.size();
	const typename std::string::size_type small_size = small.size();

	const bool valid_       = (big_size >= small_size);
	const bool starts_with_ = (big.compare(0, small_size, small) == 0);
	return valid_ && starts_with_;
}

bool ends_with(const std::string &big, const std::string &small)
{
	if (&big == &small)
		return true;
	const typename std::string::size_type big_size   = big.size();
	const typename std::string::size_type small_size = small.size();

	const bool valid_ = (big_size >= small_size);
	const bool ends_with_ =
	  (big.compare(big_size - small_size, small_size, small) == 0);
	return valid_ && ends_with_;
}

std::string replace_first(const std::string &str, const std::string &orig,
                          const std::string &rep)
{
	size_t off = str.find(orig);
	if (off == size_t(-1))
		return str;
	// get the substr before replacing position
	std::string res = str.substr(0, off);
	// append replacing string
	res += rep;
	// append the substr after replacing position
	res += str.substr(off + orig.size());
	return res;
}

std::string replace_last(const std::string &str, const std::string &orig,
                         const std::string &rep)
{
	size_t      off = str.rfind(orig);
	// get the substr before replacing position
	std::string res = str.substr(0, off);
	// append replacing string
	res += rep;
	// append the substr after replacing position
	res += str.substr(off + orig.size());
	return res;
}

std::string replace_all(const std::string &str, const std::string &orig,
                        const std::string &rep)
{
	size_t      last_off = 0;
	size_t      off      = 0;
	std::string res;
	while ((off = str.find(orig)) != std::string::npos)
	{
		// get the substr before replacing position
		res += str.substr(last_off, off);
		// append replacing string
		res += rep;
		// move the next one
		off += orig.size();
		last_off = off;
	}
	// append the substr after replacing position
	res += str.substr(off);
	return res;
}

} // namespace OMC