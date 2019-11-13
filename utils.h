#pragma once

#include <iomanip>
#include <iostream>

namespace utils {
	inline void ncout(double n, std::ostream& out = std::cout, int spacing = 8, std::string extra = "", bool nline = false) {
		int s = 1;
		int ni = static_cast<int>(n);
		
		if (ni < 0) ++s;

		do {
			// reduce n;
			ni /= 10;
			++s;
		} while (ni != 0);

		//std::cout << "@" << s << "@";

		out << std::setw(spacing - s) << " ";

		out << n << " " << extra;
		if (nline)
			out << std::endl;
	}

	inline int ndigits(int n) {
		int s = 0;
		do {
			++s;
			n /= 10;
		} while (n != 0);

		return s;
	}

	inline std::vector<std::vector<double>> power_set(std::vector<double> in) {
		std::vector<std::vector<double>> nset;
		const int pow_set_size = static_cast<const int>(pow(2, in.size()));
		for (int counter = 0; counter < pow_set_size; ++counter) {
			std::vector<double> vset;
			for (int j = 0; j < in.size(); j++) {
				if (counter & (1 << j))
					vset.push_back(in[j]);
			}
			nset.push_back(vset);
		}

		return nset;
	}

	inline int fact(int n)
	{
		int result = 1;
		while (n > 1) {
			result *= n--;
		}
		return result;
	}

	template <typename T>
	void permutation(std::vector<T> v)
	{
		std::sort(v.begin(), v.end());
		do {
			std::copy(v.begin(), v.end(), std::ostream_iterator<T>(std::cout, " "));
			std::cout << std::endl;
		} while (std::next_permutation(v.begin(), v.end()));
	}

	// FROM: http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.353.930&rep=rep1&type=pdf
	// next_mapping
	template < class BidirectionalIterator, class T >
	bool next_mapping(BidirectionalIterator first, BidirectionalIterator last, T first_value, T last_value) {
		if (last == first) {
			return false;
		}
		do {
			if (++(*(--last)) != last_value) {
				return true;
			}
			*last = first_value;
		} while (last != first);
		return false;
	}	template < class BidirectionalIterator, class T >
	bool prev_mapping(BidirectionalIterator first, BidirectionalIterator last, T first_value, T last_value)
	{
		if (last == first) {
			return false;
		}
		--last_value;
		do {
			if (*(--last) != first_value) {
				--(*last);
				return true;
			}
			*last = last_value;
		} while (last != first);
		return true;
	}

	template <typename T>
	std::vector<std::vector<T>> permutation_mapk(std::vector<T>& v, std::size_t count) {
		std::vector<std::vector<T>> pmk;
		std::vector<double> ffs(v.begin(), v.end());
		std::vector<std::vector<double>::const_iterator> v_iter(count, ffs.begin());
		do {
			std::vector<double> hey;
			for(auto pmk2 : v_iter) {
				hey.push_back(*pmk2);
			}
			pmk.push_back(hey);
		} while (next_mapping(v_iter.begin(), v_iter.end(), ffs.begin(), ffs.end()));

		return pmk;
	}

	template <typename T>
	std::vector<std::vector<T>> combination_k(const std::vector<T>& v, std::size_t count)
	{
		std::vector<std::vector<T>> cnk;
		assert(count <= v.size());
		std::vector<bool> bitset(v.size() - count, 0);
		bitset.resize(v.size(), 1);

		do {
			std::vector<T> vs;
			for (std::size_t i = 0; i != v.size(); ++i) {
				if (bitset[i]) {
					vs.push_back(v[i]);
				}
			}
			cnk.push_back(vs);
		} while (std::next_permutation(bitset.begin(), bitset.end()));

		return cnk;
	}

	inline bool increase(std::vector<bool>& bs)
	{
		for (std::size_t i = 0; i != bs.size(); ++i) {
			bs[i] = !bs[i];
			if (bs[i] == true) {
				return true;
			}
		}
		return false; // overflow
	}

	template <typename T>
	void powerset(const std::vector<T>& v)
	{
		std::vector<bool> bitset(v.size());

		do {
			for (std::size_t i = 0; i != v.size(); ++i) {
				if (bitset[i]) {
					std::cout << v[i] << " ";
				}
			}
			std::cout << std::endl;
		} while (increase(bitset));
	}
}
