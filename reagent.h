#pragma once

#include "defines.h"

#include "element.h"

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
		for (auto &el : _self) {
			if (el().q != 0) {
				s += el().n;
				// change this if you want to see quantity on 1
				if (el().q > 1) {
					s += std::to_string(el().q);
				}
			}
		}

		return s;
	}

	double mass() const {
		return _mass;
	}
};