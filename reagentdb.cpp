#include "reagentdb.h"

const std::vector<Reagent>& ReagentDB::operator()() const {
	return _self;
}

const Reagent& ReagentDB::operator()(int idx) const {
	return _self[idx];
}

void ReagentDB::insert(const Reagent& reagent) {
	_self.push_back(reagent);
}

void ReagentDB::remove(int idx) {
	_self.erase(_self.begin() + idx);
}

void ReagentDB::trim(const Reagent& compound) {
	auto it = _self.begin();
	while (it != _self.end()) {
		bool useful = false;
		for (auto re : (*it)()) {
			for (auto e : compound()) {
				if (e().n == re().n)
					useful = true;
			}
		}
		if (!useful) {
			it = _self.erase(it);
		}
		else {
			++it;
		}
	}
}

std::vector<int> ReagentDB::quanitfy(std::vector<std::string> names) {
	std::vector<int> temp;
	if (!names.empty())
		temp.reserve(names.size());

	for (auto n : names) {
		int atom_count = 0;
		for (auto r : _self) {
			for (auto ri : r()) {
				if (ri().n == n)
					atom_count += ri().q;
			}
		}

		temp.push_back(atom_count);
	}

	return temp;
}

int ReagentDB::validate(const Reagent& compound) {
	int index = 0;
	for (auto e : compound()) {
		bool possible_element = false;
		bool possible_atomcnt = false;
		int ec = 0;
		for (auto ri : _self) {
			for (auto re : ri()) {
				if (e().n == re().n) {
					possible_element = true;
					ec += re().q;
				}
			}
		}
		if (ec >= e().q) {
			possible_atomcnt = true;
		}
		if (!possible_element || !possible_atomcnt) {
			if (!possible_element) {
				//std::cout << "Reconstruction impossible, missing: " << e.first << std::endl;
				return 1 + index;
			}
			if (!possible_atomcnt && possible_element) {
				//std::cout << "Reconstruction impossible, not enough " << e.first << "; Present " << ec << " out of " << e.second << " atoms" << std::endl;
				return -1 - index;
			}
			//NLINE

			//return false;
		}
		++index;
	}

	//std::cout << "Reconstruction possible, all elements present" << std::endl;
	//NLINE
	return 0;
}
