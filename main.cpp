#include "defines.h"
#include "utils.h"
#include "input.h"

#include "elementdb.h"
#include "solutiondb.h"

int main(int argc, char **argv) {
	// start program
	std::chrono::high_resolution_clock::time_point start_time, current_time;
	start_time = std::chrono::high_resolution_clock::now();

	ElementDB::init();
	SolutionDB sdb;

	Input in(argc, argv);

	auto r = in().r;
	auto rdb = in().rdb;

	auto mp = in.matrix();

	Eigen::FullPivLU<Eigen::MatrixXd> lu;

	lu.compute(mp.first);
	auto K = lu.kernel();
	auto S = lu.solve(mp.second);

	in.validate_weights(static_cast<int>(K.cols()));
	auto weights = in.load(WEIGHTSCACHEPATH);

	int sols = 0;

	sdb.insert(S);

	for (auto &ps : weights) {
		Eigen::VectorXd p(S.size());
		//auto p = nS;
		for (int i = 0; i < S.size(); ++i) {
			double x = S(i);
			p(i) = math::round<double>(x);
		}
		for (int i = 0; i < K.cols(); ++i) {
			p += ps[i] * K.col(i);
		}

		//check for negative solutions
		bool ffs = false;
		for (int i = 0; i < rdb().size(); ++i) {
			if (p(i) < 0)
				ffs = true;
		}
		if (ffs == true)
			continue;

		sdb.insert(p);
	}

	std::cout << "Total solutions: " << sdb().size() << std::endl;
	std::vector<int> ranks;
	for (auto s : sdb()) {
		ranks.push_back(s.second.score());
	}
	std::sort(ranks.begin(), ranks.end());
	std::cout << "Best score is: " << ranks[0] << std::endl;

	// trim based on samples
	sdb.trim(ranks[in().samples]);

	current_time = std::chrono::high_resolution_clock::now();
	auto clockval = std::chrono::duration<double>(current_time - start_time).count();
	std::cout << ">> Points computed in: " << clockval << "s" << std::endl << std::endl;

	sdb.exportcsv(r, rdb, in().options.use_input_cache);

	return 0;
}