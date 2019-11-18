#include "solution.h"

Solution::Solution(Eigen::VectorXd vec)
{
	for (int i = 0; i < vec.size(); ++i) {
		insert(vec(i));
	}

	validate();

	// Compute key
	std::stringstream ss;
	int nz = vec.size();
	for (auto s : _self) {
		if (s == 0)
			nz--;
	}
	ss << nz << "!";
	for (int i = 0; i < _self.size() - 1; ++i) {
		ss << _self[i] << "#";
	}
	ss << _self[_self.size() - 1];

	_key = ss.str();
	// End Compute Key

	// Compute score
	auto rv = math::rational_split(_selfr);
	_score = abs(std::accumulate(rv.second.begin(), rv.second.end(), 1, math::lcm));
}

const std::vector<double>& Solution::operator()() const
{
	return _self;
}

const double & Solution::operator()(int idx) const
{
	return _self[idx];
}

void Solution::insert(const double& value, int maxden)
{
	_self.push_back(value);
	_selfr.push_back(math::rational(value, maxden));
}

void Solution::remove(int idx)
{
	_self.erase(_self.begin() + idx);
}

std::string Solution::key()
{
	return _key;
}

std::vector<std::pair<int, int>> Solution::rational()
{
	return _selfr;
}

void Solution::stringify()
{
	for (auto& s : _self) {

	}
}

std::string Solution::str()
{
	return _tag;
}

int Solution::score()
{
	return _score;
}

void Solution::validate(int precision, int maxdenum)
{
	for (auto& s : _self) {
		auto ra = math::rational(s, maxdenum);
		double y = ((double)ra.first / (double)ra.second);
		s = math::round<double>(y, precision);
	}
}
