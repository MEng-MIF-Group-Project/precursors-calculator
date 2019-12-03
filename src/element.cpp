#include "element.h"
#include "elementdb.h"

Element::Element() = default;

Element::Element(const ElementData& d) {
	_self = d;
}

Element::Element(std::string s, double i, int e) {
	_self = ElementDB::map()[s]();

	_self.n = s;
	_self.q = i;
	_self.s = e;
}

const ElementData& Element::operator()() const {
	return _self;
}
