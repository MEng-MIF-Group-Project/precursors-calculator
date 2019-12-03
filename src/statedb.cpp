#include "statedb.h"

StateDB::StateDB(const std::string& f) {
	std::fstream fin(f);

	std::string l;
	int i = 0;
	boost::char_separator<char> s0(" ");
	while (std::getline(fin, l)) {
		boost::tokenizer<boost::char_separator<char>> t(l, s0);

		auto it = t.begin();
		std::string e = *it;
		++it;
		_keys.emplace(e);

		while (it != t.end()) {
			std::istringstream(*it) >> i;
			_self.emplace(e, i);

			++it;
		}
	}
}

const std::multimap<std::string, int>& StateDB::operator()() const {
	return _self;
}

state_db_map_range StateDB::operator()(const std::string& k) {
	return _self.equal_range(k);
}

std::vector<std::vector<int>> StateDB::populate(std::vector<std::string> e) {
	std::vector<int> charge_count;
	charge_count.reserve(e.size());

	for (const auto& n : e) {
		charge_count.push_back(_self.count(n));
	}

	int mc = 0;
	for (auto cc : charge_count) {
		if (cc > mc)
			mc = cc;
	}

	std::vector<int> charges;
	charges.reserve(charge_count.size());

	for (int i = 0; i < charge_count.size(); i++)
		charges.push_back(_self.find(e[i])->second);

	std::vector<std::vector<int>> last;
	last.push_back(charges);
	for (int i = charge_count.size() - 1; i >= 0; i--) {
		std::vector<int> range;

		auto r = _self.equal_range(e[i]);
		for (auto it = r.first; it != r.second; ++it) {
			range.push_back(it->second);
		}

		std::vector<std::vector<int>> current = last;
		for (auto& j : range) {
			for (int k = 0; k < last.size(); k++) {
				std::vector<int> nv = current[k];

				nv[i] = j;

				bool s = false;
				for (int o = 0; o < last.size(); o++) {
					int sc = 0;
					for (int m = 0; m < current[o].size(); m++) {
						if (nv[m] == current[o][m])
							sc++;
					}

					if (sc == current[o].size()) {
						s = true;
					}
				}

				if (!s) {
					current.push_back(nv);
				}
			}
		}

		last = current;
	}

	return last;
}

void StateDB::print() {
	std::cout << "---------------- PRINTING OXIDATION STATES ----------------" << std::endl;
	for (auto ei : _keys) {
		std::cout << ei;
		std::cout << std::setw(4 - ei.size()) << " ";
		std::cout << ">> ";
		auto range = _self.equal_range(ei);
		for (auto it = range.first; it != range.second; ++it) {
			if (it->second >= 0) std::cout << "+";
			std::cout << it->second << " ";
		}
		std::cout << std::endl;
	}
	std::cout << "--------------------------- END --------------------------" << std::endl;
	std::cout << std::endl;
}
