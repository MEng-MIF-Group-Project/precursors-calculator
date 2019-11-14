#pragma once

#include "defines.h"

typedef std::pair<std::multimap<std::string, int>::iterator, std::multimap<std::string, int>::iterator> state_db_map_range;

class StateDB
{
	std::multimap<std::string, int> _self;
	std::set<std::string> _keys;
public:
	StateDB(const std::string& f = "data/oxidation_states.txt");;

	const std::multimap<std::string, int>& operator()() const;
	state_db_map_range operator()(const std::string& k);

	std::vector<std::vector<int>> populate(std::vector<std::string> e);

	void print();
};