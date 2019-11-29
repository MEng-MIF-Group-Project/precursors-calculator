#pragma once

#include "defines.h"

#include "element.h"


#include <algorithm>
#include <vector>

class Reagent
{
	std::vector<Element> _self;
	double _mass = 0;

public:
	Reagent();
	Reagent(std::vector<std::string> name, std::vector<double> atom_count);

	const std::vector<Element>& operator()() const;
	const Element& operator()(int idx) const;

	void insert(const Element& element);
	void remove(int idx);

	std::string str() {
		std::string s;

		std::sort(_self.begin(), _self.end(), [](const auto& e1, const auto& e2) {
			double a, b;
			std::stringstream(e1().el_neg) >> a;
			std::stringstream(e2().el_neg) >> b;
	
			return (a < b);
		});
		
		for (auto &el : _self) {
			if (el().q != 0) {
				s += el().n;
				// change this if you want to see quantity on 1
				if (el().q > 1) {
					s += std::to_string(el().q);
				}

				s += "-";
			}
		}

		return s.substr(0, s.size() - 1);
	}

	double mass(double qualifier = 1.0f) const {
		return qualifier * _mass;
	}
};