#include "solutiondb.h"

#include "input.h"

const std::map<std::string, Solution>& SolutionDB::operator()() const
{
	return _self;
}

Solution SolutionDB::operator()(std::string key)
{
	return _self[key];
}

void SolutionDB::insert(Solution& solution)
{
	_self[solution.key()] = solution;
}

void SolutionDB::insert(Eigen::VectorXd vec)
{
	Solution s(vec);
	_self[s.key()] = s;
}

void SolutionDB::remove(std::string key)
{
	_self.erase(key);
}

void SolutionDB::trim(int score)
{
	auto it = _self.begin();
	while (it != _self.end()) {
		if (it->second.score() > score) {
			_self.erase(it++);
		}
		else {
			it++;
		}
	}
}

void SolutionDB::exportcsv(Reagent r, ReagentDB rdb, bool incache)
{
	auto start_time = std::chrono::high_resolution_clock::now();

	std::ofstream g(r.str() + ".csv");

	g << "Key,";
	std::string header = "";
	for (auto r : rdb()) {
		header += r.str() + ",";
	}
	header.erase(header.size() - 1);
	g << header;
	g << ",Score";
	g << std::endl;

	for (auto s : _self) {
		g << s.first << ",";
		for (int i = 0; i < s.second().size() - 1; ++i) {
			g << s.second()[i] << ",";
		}
		g << s.second()[s.second().size() - 1];
		g << "," << s.second.score();
		g << std::endl;
	}
	g.close();

	if (incache == false) {
		std::ofstream icg(INPUTCACHEPATH);
		icg << "Reserved" << std::endl;
		icg << r.str() << std::endl;
		for (auto pr : rdb()) {
			icg << pr.str() << " ";
		}
		icg << std::endl;
	}

	auto current_time = std::chrono::high_resolution_clock::now();
	auto clockval = std::chrono::duration<double>(current_time - start_time).count();
	std::cout << std::endl << ">> CSV operation completed in: " << clockval << "s" << std::endl;
}
