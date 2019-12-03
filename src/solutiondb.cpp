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
	insert(s);
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

void SolutionDB::exportcsv(Reagent r, ReagentDB rdb, Input::IOdata io)
{
	auto start_time = std::chrono::high_resolution_clock::now();
	std::string filename = std::to_string(io.mode);
	filename += r.str() + ".csv";

	std::ofstream g(filename);

	// assume 1 mol = 66 grams, want 0.5 grams instead
	// 66
	double mult = io.dmass / r.mass();

	std::string header = "Key,";
	if (io.mode == 0) {
		header += "Name,Mass,";
	}
	else {
		header += "Name,In_Mass,";
	}
	for (auto r : rdb()) {
		header += r.str() + ",";
		if (io.mode == 1) {
			header += r.str() + "_Mass,";
		}
	}
	header.erase(header.size() - 1);
	g << header;
	g << ",Score";
	g << std::endl;

	for (auto s : _self) {
		std::stringstream line;
		bool bad_q = false;
		line << s.first << ",";
		std::string name = "";

		if (io.mode == 0) {
			auto ra = s.second.rational();
			double mass = 0;
			for (int i = 0; i < ra.size(); ++i) {
				auto r = rdb()[i];
				bool in_use = false;
				int ratio = boost::math::round<int>(static_cast<double>(ra[i].first * s.second.score()) / static_cast<double>(ra[i].second));
				if (ratio > 0) {
					name += r.str();
					mass += r.mass() * ratio;

					if (ratio >= 2) {
						name += std::to_string(ratio);

						if (ratio >= 20) {
							bad_q = true;
						}
					}

					name += "-";
				}
				
			}
			line << name.substr(0, name.size() - 1) << ",";
			line << mass << ",";

			//std::cout << "Name: " << name << std::endl;
		}
		else {
			line << r.str() << ",";
			line << r.mass() * mult << ",";
		}
		if (bad_q == true)
			continue;

		for (int i = 0; i < s.second().size() - 1; ++i) {
			line << s.second()[i] << ",";
			if (io.mode == 1) {
				line << (rdb()[i].mass() * s.second()[i] * mult) << ",";
			}
		}
		line << s.second()[s.second().size() - 1] << ",";
		if (io.mode == 1) {
			line << (rdb()[s.second().size() - 1].mass() * s.second()[s.second().size() - 1]) * mult << ",";
		}
		line << s.second.score();
		line << std::endl;

		if (bad_q == false) {
			g << line.str();
		}
	}
	g.close();

	if (io.options.use_input_cache == false) {
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
