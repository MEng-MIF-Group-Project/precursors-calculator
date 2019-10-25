#pragma once

#include <unordered_map>

#include "defines.h"
#include "element.h"


class ElementDB
{
	static std::unordered_map<std::string, Element> _self;
	//static std::set<std::string> _keys;

public:
	
	/*ElementDB(const std::string& f = "data/element_data.txt") {
		std::fstream fin(f);

		std::string l;
		int i = 0;
		boost::char_separator<char> s0(" ");
		while (std::getline(fin, l)) {
			if (l[0] == '#')
				continue;

			boost::tokenizer<boost::char_separator<char>> t(l, s0);

			auto it = t.begin();
			std::string e = *it;
			_keys.emplace(e);

			// keep advancing token iterator from pos 1
			ElementData ed;
			ed.n = e;
			// new stuff, needs re-organizing later
			// only care about mass for now
			ed.symbol = e;
			ed.name = *it;
			std::stringstream(*it) >> ed.z;
			std::stringstream(*it) >> ed.mass;
			ed.r_cov = *it;
			ed.e_affinity = *it;
			ed.p_eig = *it;
			ed.s_eig = *it;
			ed.abundance = *it;
			ed.el_neg = *it;
			ed.ion_pot = *it;
			ed.dipol = *it;
		
			_self.emplace(e, Element(ed));
		}
	}
	*/

	static void init(const std::string& f = "data/element_data.txt") {
		std::fstream fin(f);

		std::string l;
		int i = 0;
		boost::char_separator<char> s0(" ");
		while (std::getline(fin, l)) {
			if (l[0] == '#')
				continue;

			boost::tokenizer<boost::char_separator<char>> t(l, s0);

			auto it = t.begin();
			//_keys.emplace(e);

			// keep advancing token iterator from pos 1
			ElementData ed;
			ed.n = *it;
			ed.q = 1;
			// new stuff, needs re-organizing later
			// only care about mass for now
			ed.symbol = *it; ++it;
			ed.name = *it; ++it;
			std::stringstream(*it) >> ed.z; ++it;
			std::stringstream(*it) >> ed.mass; ++it;
			ed.r_cov = *it; ++it;
			ed.e_affinity = *it; ++it;
			ed.p_eig = *it; ++it;
			ed.s_eig = *it; ++it;
			ed.abundance = *it; ++it;
			ed.el_neg = *it; ++it;
			ed.ion_pot = *it; ++it;
			ed.dipol = *it; ++it;

			_self.emplace(ed.n, Element(ed));
		}
	}

	//std::unordered_map<std::string, Element>& operator()() const {
	//	return _self;
	//}

	static std::unordered_map<std::string, Element>& map() {
		return _self;
	}
};
