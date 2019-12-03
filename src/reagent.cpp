#include "reagent.h"
#include "elementdb.h"

#include <iostream>

Reagent::Reagent() = default;

Reagent::Reagent(std::vector<std::string> name, std::vector<double> atom_count) {
	for (auto i = 0; i < name.size(); i++) {
		_self.emplace_back(name[i], atom_count[i]);
		_mass += _self.back()().mass * _self.back()().q;
		//std::cout << ElementDB::map()[name[i]]().mass << " " << std::endl;
	}
}

const std::vector<Element>& Reagent::operator()() const {
	return _self;
}

const Element& Reagent::operator()(int idx) const {
	return _self[idx];
}

void Reagent::insert(const Element& element) {
	_self.push_back(element);
	_mass += element().mass * element().q;
}

void Reagent::remove(int idx) {
	_self.erase(_self.begin() + idx);
}
