#pragma once

#include "defines.h"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/math/special_functions/round.hpp>
#include <boost/math/common_factor.hpp>
#include <boost/rational.hpp>

#include <Eigen/Dense>

#include <sqlite3.h>
#include <chrono>
#include <numeric>
#include <algorithm>
#include <iterator>
#include <map>

namespace math {
	// return lcm of 2 numbers
	inline int lcm(int a, int b) {
		int temp = boost::math::gcd(a, b);

		return temp ? (a / temp * b) : 0;
	}

	// David Eppstein / UC Irvine / 8 Aug 1993
	// find rational approximation to given real number
	// TODO: there are some quirks that need fixing here but superlow priority
	inline std::pair<int, int> rational(double x, long maxden) {
		long m[2][2];
		double startx;
		long ai;

		startx = x;

		/* initialize matrix */
		m[0][0] = m[1][1] = 1;
		m[0][1] = m[1][0] = 0;

		/* loop finding terms until denom gets too big */
		while (m[1][0] * (ai = (long)x) + m[1][1] <= maxden) {
			long t;
			t = m[0][0] * ai + m[0][1];
			m[0][1] = m[0][0];
			m[0][0] = t;
			t = m[1][0] * ai + m[1][1];
			m[1][1] = m[1][0];
			m[1][0] = t;
			if (x == (double)ai) break;     // AF: division by zero
			x = 1 / (x - (double)ai);
			if (x > (double)0x7FFFFFFF) break;  // AF: representation failure
		}

		/* now remaining x is between 0 and 1/ai */
		/* approx as either 0 or 1/m where m is max that will fit in maxden */
		/* first try zero */
		std::pair<int, int> fp = { m[0][0], m[1][0] };
		double er = abs(startx - ((double)m[0][0] / (double)m[1][0]));
		/* now try other possibility */
		ai = (maxden - m[1][1]) / m[1][0];
		m[0][0] = m[0][0] * ai + m[0][1];
		m[1][0] = m[1][0] * ai + m[1][1];

		std::pair<int, int>  lp = { m[0][0], m[1][0] };
		double ner = abs(startx - ((double)m[0][0] / (double)m[1][0]));

		return fp;
		if (ner > er || ner == er) {
			return fp;
		}
		else {
			return lp;
		}
	}

	// split a vector of ints into a pair of vectors, TODO: move this to utils or something
	inline std::pair<std::vector<int>, std::vector<int>> rational_split(std::vector<std::pair<int, int>> rs) {
		std::vector<int> n, d;

		for (auto r : rs) {
			n.push_back(r.first);
			d.push_back(r.second);
		}

		return { n, d };
	}

	// precision-based round
	template<typename T>
	inline T round(T& x, int precision = 3) {
		T n;
		std::stringstream ss;
		ss << std::setprecision(precision) << x;
		ss >> n;
		return n;
	}

	// generate combinations of nk with repetitions
	template<typename V, typename Callable>
	void for_each_combination(V &v, size_t gp_sz, Callable f) {
		V gp(gp_sz);
		auto total_n = std::pow(v.size(), gp.size());
		for (auto i = 0; i < total_n; ++i) {
			auto n = i;
			for (auto j = 0ul; j < gp.size(); ++j) {
				gp[gp.size() - j - 1] = v[n % v.size()];
				n /= v.size();
			}
			f(gp);
		}
	}

	// returns n! given n
	/*inline unsigned long long factorial(int n)
	{
		unsigned long long result = 1;
		while (n > 1) {
			result *= n--;
			std::cout << "N: " << n << "Fact: " << result << std::endl;
		}
		return result;
	}*/
}