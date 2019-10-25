#pragma once

#include "reagent.h"

class ReagentDB
{
	std::vector<Reagent> _self;
public:
	const std::vector<Reagent>& operator()() const;
	const Reagent& operator()(int idx) const;

	void insert(const Reagent& reagent);
	void remove(int idx);

	void trim(const Reagent& compound);
	std::vector<int> quanitfy(std::vector<std::string> names);
	int validate(const Reagent& compound);
};