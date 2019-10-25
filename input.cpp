#include "input.h"

Input::Input(std::string f) {
	parse(f);
}

void Input::parse(std::string f) {
	// Line split
	std::vector<std::string> ld;
	std::ifstream file(f);

	std::string line;
	while (std::getline(file, line)) {
		ld.push_back(line);
	}

	// Line Parsing
	// 0. Reserved
	// 1. List of Elements
	// 2. Available reagents
	boost::char_separator<char> s0(" ");
	boost::char_separator<char> s1(":");
	boost::tokenizer<boost::char_separator<char>> t1(ld[1], s0);

	// Retrieve name of elements
	for (const auto& t : t1) {
		boost::tokenizer<boost::char_separator<char>> t2_0(t, s1);
		for (const auto& m : t2_0) {
			int c = 0;
			const auto r = m.find('#');
			if (r == std::string::npos) {
				_self.els.push_back(m);
				_self.amounts.push_back(1);
			}
			else {
				const std::string n = m.substr(0, r);

				std::stringstream ss;
				ss << m.substr(r + 1);
				ss >> c;

				_self.els.push_back(n);
				_self.amounts.push_back(c);
			}
		}
	}

	boost::tokenizer<boost::char_separator<char>> t2(ld[2], s0);
	for (const auto& t : t2) {
		boost::tokenizer<boost::char_separator<char>> t2_0(t, s1);
		Reagent reagent;
		for (const auto& m : t2_0) {
			int c = 0;
			const auto r = m.find('#');
			if (r == std::string::npos) {
				reagent.insert(Element(m, 1));
			}
			else {
				const std::string n = m.substr(0, r);

				std::stringstream ss;
				ss << m.substr(r + 1);
				ss >> c;

				reagent.insert(Element(n, c));
			}
		}
		_self.rdb.insert(reagent);
	}
}

const Input::IOdata& Input::operator()() const {
	return _self;
}
