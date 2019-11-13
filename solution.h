#pragma once

#include "defines.h"
#include "math.h"

class Solution
{
	std::vector<double> _self;
	std::vector<std::pair<int, int>> _selfr;
	std::string _key;
	int _score;

public:
	Solution() = default;
	Solution(Eigen::VectorXd vec);

	const std::vector<double>& operator()() const;
	const double& operator()(int idx) const;

	void insert(const double& value, int maxden = 99);
	void remove(int idx);

	std::string key();
	std::vector<std::pair<int, int>> rational();
	int score();
	
	void validate(int precision = 3, int maxdenum = 99);
};