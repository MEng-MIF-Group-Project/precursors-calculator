#pragma once

#include "solution.h"

#include "input.h"
#include "reagentdb.h"

#include <map>

class SolutionDB
{
	std::map<std::string, Solution> _self;
public:
	const std::map<std::string, Solution>& operator()() const;
	Solution operator()(std::string key);

	void insert(Solution& solution);
	void insert(Eigen::VectorXd vec);
	void remove(std::string key);

	void trim(int score);

	void exportcsv(Reagent r, ReagentDB rdb, Input::IOdata io);
};