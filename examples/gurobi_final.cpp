#include "defines.h"

#include "input.h"
#include "modeldb.h"
#include "elementdb.h"
#include "utils.h"

#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

int main() {
    const ElementDB edb;
	ElementDB::init();

	StateDB sdb;
	ModelDB mdb(edb);

	std::cout.precision(4); 
	std::cout.setf(ios::fixed);

	const Input in(INPUT);
	const auto els = in().els;
	const auto amounts = in().amounts;
	const auto rdb = in().rdb;

	const auto coefficients = sdb.populate(els);


	std::stringbuf nbuf;
	std::ostream log(&nbuf);
	std::ofstream g("out.txt");


	std::cout << std::endl;

	/*
	 * First parse of mdb, generates solutions based on possible element states
	 * 1: ensures charge parity, AX + BY = 0, where X and Y are states and A and B are quantities of elements
	 * 2: ensures at least one element is present, A + B >= 1 where A and B are quantities of elements
	 * returns a 3D array containing:
	 * > 3rd layer contains all models solved
	 * > 2nd layer contains all solutions to the current model
	 * > 1st layer contains all variables from the current solution
	 */
	for (const auto c : coefficients) {
		std::vector<std::vector<std::pair<std::string, double>>> lhs, rhs;
		std::vector<char> sense = { '=', '>' };

		std::vector<std::pair<std::string, double>> temp1, temp0;
		for (int i = 0; i < c.size(); ++i) {
			temp1.push_back({ els[i], c[i] });
			temp0.push_back({ els[i], 1 });
		}

		lhs.push_back(temp1);
		lhs.push_back(temp0);
		
		rhs.push_back({ {"", 0} });
		rhs.push_back({ {"", 1} });

		// var type is int as we are looking for quantities of elements (mols i think)
		//mdb.insert(lhs, rhs, sense, 'I', { 0, 10 });
	};

	std::stringbuf test;
	std::ostream hello(&test);
	//mdb.solve(rdb, log);

	//g << log.rdbuf();

	// TODO: maybe make solve clear data as well? ask mike
	// empty model data so mdb can be reused later
	//mdb.clear();

	//clear redundant info from cmd
	system("CLS");

	// TODO : take a loot and split this later
	/*
	const char var = 'A';
	std::vector<double> solution = { 1, 1, 1, 1 };
	std::vector<std::vector<std::pair<string, double>>> lhs, rhs;
	std::vector<char> sense;
			
	// result from previous solve as the new rhs
	for (int i = 0; i < solution.size(); ++i) {
		sense.push_back('=');
			
		rhs.push_back({ {"", solution[i]} });
	}

	const int nv = rdb().size();

	std::vector<std::string> vn;
	for (int i = 0; i < nv; ++i) {
		vn.push_back(std::string(1, var + i));
	}

	// create lhs from available reagents
	// TODO: needs combinations or something for multiple solutions but will probably just swithc to google-ot anyway
	for (int i = 0; i < els.size(); ++i) {
		std::vector<std::pair<string, double>> lhs_row;
		for (int j = 0; j < nv; ++j) {
			double val = 0;
			for (auto &k : rdb(j)()) {
				if (k().n == els[i]) {
					val = k().q;
					// TODO: think about this break
					break;
				}

			}

			lhs_row.push_back({ vn[j], val });
		}
		lhs.push_back(lhs_row);
	}

	// type is S/C as we want continous values
	// range of values is between 0 and 10 inclusive, ask mike about this
	mdb.insert(lhs, rhs, sense, 'C', { 0.00, 10 }, Reagent(els, solution));
	*/
	// combinations of reagents
	std::vector<double> inset;
	for (int i = 0; i < rdb().size(); ++i) {
		// TODO: do pretrimming here to save speed
		inset.push_back(i);
	}
	auto pset = utils::power_set(inset);

	const char var = 'A';
	std::vector<std::string> vn;


	for (int i = 0; i < inset.size(); ++i) {
		vn.push_back(std::string(1, var + i));
	}

	for (auto &ps : pset) {
		std::vector<std::vector<std::pair<string, double>>> lhs, rhs;
		std::vector<char> sense;

		// nr of cols
		
		for (int i = 0; i < els.size(); ++i) {
			sense.push_back('=');

			std::vector<std::pair<std::string, double>> lhs_row;
			for (int j = 0; j < ps.size(); ++j) {
				double val = 0;
				for (auto &k : rdb(static_cast<int>(ps[j]))()) {
					if (k().n == els[i]) {
						val = k().q;
						break;
					}

				}
				
				lhs_row.push_back({ vn[j], val });
			}

			lhs.push_back(lhs_row);
			rhs.push_back({ {"", amounts[i]} });
		}

		// per column
		std::vector<std::string> col_labels;
		for (int j = 0; j < ps.size(); ++j) {
			Reagent r = rdb(static_cast<int>(ps[j]));

			int non_zero = 0;
			for (int i = 0; i < els.size(); ++i) {
				if (lhs[i][j].second == 0)
					non_zero++;
			}

			// remove all zeero colums
			if (non_zero == els.size()) {
				for (int i = 0; i < els.size(); ++i) {
					try {
						if (lhs[i].size() > j)
							lhs[i].erase(lhs[i].begin() + j);
						else {
							lhs[i].erase(lhs[i].begin() + lhs[i].size() - 1);
						}
					} catch (...) {
						
					}
					
				}

				continue;
			}

			col_labels.push_back(r.str());
		}

		mdb.insert(lhs, rhs, sense, 'C', { 0.0005, 10.0 }, col_labels, Reagent(els, amounts));
	}

	//mdb.property.repeat = true;
	mdb.property.use_edb = true;

	mdb.solve(rdb, log);

	std::cout << "OUTPUT >> out.txt" << std::endl;

	g << log.rdbuf();

	std::cout << std::endl;
	std::cout << "DONE" << std::endl;

	return 0;
}