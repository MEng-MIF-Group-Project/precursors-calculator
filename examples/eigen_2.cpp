#include "defines.h"

#include "input.h"
#include "modeldb.h"
#include "elementdb.h"
#include "utils.h"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/math/special_functions/round.hpp>
#include <boost/math/common_factor.hpp>
#include <boost/rational.hpp>

#include <Eigen/Dense>

#include <sqlite3.h>
#include <chrono>
#include <numeric>
#include <algorithm>
#include <iterator>
#include <map>

using namespace Eigen;
using namespace std;

static int sql_callback(void *NotUsed, int argc, char **argv, char **azColName) {
	int i;
	for (i = 0; i < argc; i++) {
		printf("%s = %s\n", azColName[i], argv[i] ? argv[i] : "NULL");
	}
	printf("\n");
	return 0;
}

std::pair<std::pair<std::vector<std::string>, std::vector<double>>, ReagentDB> parse_input(std::string stoichs, std::string precursors) {
	// Line Parsing
	// 0. Reserved
	// 1. List of Elements
	// 2. Available reagents
	boost::char_separator<char> s0(" ");
	boost::char_separator<char> s1(":");
	boost::tokenizer<boost::char_separator<char>> t1(stoichs, s0);

	ReagentDB rdb;
	std::vector<std::string> els;
	std::vector<double> amounts;

	// Retrieve name of elements
	for (const auto& t : t1) {
		boost::tokenizer<boost::char_separator<char>> t2_0(t, s1);
		for (const auto& m : t2_0) {
			int c = 0;
			const auto r = m.find('#');
			if (r == std::string::npos) {
				els.push_back(m);
				amounts.push_back(1);
			}
			else {
				const std::string n = m.substr(0, r);

				std::stringstream ss;
				ss << m.substr(r + 1);
				ss >> c;

				els.push_back(n);
				amounts.push_back(c);
			}
		}
	}

	boost::tokenizer<boost::char_separator<char>> t2(precursors, s0);
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
		rdb.insert(reagent);
	}
	
	
	return {{els, amounts}, rdb};
}

void makedb(ReagentDB rdb, std::vector<std::vector<double>> solutions, std::string ies, bool use_input_cache = true, std::string cmd_input_stoics = "", std::string cmd_input_precursors = "") {
	auto start_time = std::chrono::high_resolution_clock::now();
	// at this point we have the solutions, export to a database
	// if we keep the structure of precursors as columns this would
	// mean a unique table for each unique scenario

	std::string db_name = "precursors";

	// table creation
	// TODO: check for existing table WITH the same precursors as the current test before dropping it
	std::string sql = "DROP TABLE IF EXISTS " + ies + ";\n";
	sql += "CREATE TABLE ";
	sql += ies;
	sql += "(\n";
	int id = 0;
	for (auto rp : rdb()) {
		sql += "\t'" + rp.str() + "'\tREAL NOT NULL,\n";
		//sql += "\t'g/Mol" + std::to_string(id++) + "'\tREAL NOT NULL,\n";
	}
	sql.erase(sql.size() - 2);
	sql += "\n);\n\n";

	// values insertion
	sql += "INSERT INTO\n\t" + ies + "(";
	id = 0;
	for (auto rp : rdb()) {
		sql += "'" + rp.str() + "',";
		//sql += "'g/Mol" + std::to_string(id++) + "',";
	}
	sql.erase(sql.size() - 1);
	sql += ")\n";
	sql += "VALUES\n";
	for (auto sol : solutions) {
		sql += "\t(";
		int id = 0;
		for (auto sv : sol) {
			sql += std::to_string(sv) + ", ";
			//sql += std::to_string(rdb()[id].mass() * sv) + ", ";
		}
		sql.erase(sql.size() - 2);
		sql += "),\n";
	}
	sql.erase(sql.size() - 2);
	sql += ";\n\n";

	//table_mass calculations
	sql += "DROP TABLE IF EXISTS " + ies + "_mass;\n";
	sql += "CREATE TABLE ";
	sql += ies;
	sql += "_mass(\n";
	//int id = 0;
	for (auto rp : rdb()) {
		sql += "\t'" + rp.str() + "'\tREAL NOT NULL,\n";
		//sql += "\t'g/Mol" + std::to_string(id++) + "'\tREAL NOT NULL,\n";
	}
	sql.erase(sql.size() - 2);
	sql += "\n);\n\n";

	// values insertion
	sql += "INSERT INTO\n\t" + ies + "_mass(";
	id = 0;
	for (auto rp : rdb()) {
		sql += "'" + rp.str() + "',";
		//sql += "'g/Mol" + std::to_string(id++) + "',";
	}
	sql.erase(sql.size() - 1);
	sql += ")\n";
	sql += "VALUES\n";
	for (auto sol : solutions) {
		sql += "\t(";
		int id = 0;
		for (auto sv : sol) {
			//sql += std::to_string(sv) + ", ";
			sql += std::to_string(rdb()[id++].mass() * sv) + ", ";
		}
		sql.erase(sql.size() - 2);
		sql += "),\n";
	}
	sql.erase(sql.size() - 2);
	sql += ";\n\n";

	// table creation done
	ofstream g(ies + ".sql");
	g << sql;
	g.close();

	// opening db
	sqlite3 *db;
	char *zErrMsg = 0;
	int rc;

	rc = sqlite3_open("precursors.db", &db);

	if (rc) {
		fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(db));
	}
	else {
		fprintf(stderr, "Opened database successfully\n");
	}

	rc = sqlite3_exec(db, sql.c_str(), sql_callback, 0, &zErrMsg);
	if (rc != SQLITE_OK) {
		fprintf(stderr, "SQL error: %s\n", zErrMsg);
		sqlite3_free(zErrMsg);
	}
	else {
		fprintf(stdout, "Records created successfully\n");
	}

	//closing db
	sqlite3_close(db);

	//caching input for safety purposes
	if (use_input_cache == false) {
		ofstream icg("last_input_cache.txt");
		icg << "Reserved" << endl;
		icg << cmd_input_stoics << endl;
		icg << cmd_input_precursors << endl;
	}

	auto current_time = std::chrono::high_resolution_clock::now();
	auto clockval = std::chrono::duration<double>(current_time - start_time).count();
	std::cout << endl << ">> Database operation completed in: " << clockval << "s" << endl;
}

int lcm2(int a, int b)
{
	int temp = boost::math::gcd(a, b);

	return temp ? (a / temp * b) : 0;
}

/*
** find rational approximation to given real number
** David Eppstein / UC Irvine / 8 Aug 1993
*/
std::pair<int, int> new_ra(double x, long maxden) {
	long m[2][2];
	double startx;
	long ai;

	startx = x;

	/* initialize matrix */
	m[0][0] = m[1][1] = 1;
	m[0][1] = m[1][0] = 0;

	/* loop finding terms until denom gets too big */
	while (m[1][0] * (ai = (long)x) + m[1][1] <= maxden) {
		long t;
		t = m[0][0] * ai + m[0][1];
		m[0][1] = m[0][0];
		m[0][0] = t;
		t = m[1][0] * ai + m[1][1];
		m[1][1] = m[1][0];
		m[1][0] = t;
		if (x == (double)ai) break;     // AF: division by zero
		x = 1 / (x - (double)ai);
		if (x > (double)0x7FFFFFFF) break;  // AF: representation failure
	}

	/* now remaining x is between 0 and 1/ai */
	/* approx as either 0 or 1/m where m is max that will fit in maxden */
	/* first try zero */
	std::pair<int, int> fp = { m[0][0], m[1][0] };
	double er = abs(startx - ((double)m[0][0] / (double)m[1][0]));
	/* now try other possibility */
	ai = (maxden - m[1][1]) / m[1][0];
	m[0][0] = m[0][0] * ai + m[0][1];
	m[1][0] = m[1][0] * ai + m[1][1];

	std::pair<int, int>  lp = { m[0][0], m[1][0] };
	double ner = abs(startx - ((double)m[0][0] / (double)m[1][0]));
	
	return fp;
	if (ner > er || ner == er) {
		return fp;
	}
	else {
		return lp;
	}
}

std::vector<std::pair<int, int>> to_rat(std::vector<double> solution, int maxdenum = 99) {
	std::vector<std::pair<int, int>> rs;
	for (auto s : solution) {
		if (s == 0)
			continue;

		auto p = new_ra(s, maxdenum);

		rs.push_back({ p.first, p.second });
	}
	return rs;
}

std::pair<std::vector<int>, std::vector<int>> split_rat(std::vector<std::pair<int, int>> rs) {
	std::vector<int> nums;
	std::vector<int> denums;

	for (auto r : rs) {
		nums.push_back(r.first);
		denums.push_back(r.second);
	}

	return { nums, denums };
}

int score(std::vector<double> solution, int maxdenum = 99) {
	auto ras = to_rat(solution);
	auto sra = split_rat(ras);
	int lcm = abs(std::accumulate(sra.second.begin(), sra.second.end(), 1, lcm2)); // running on denumerators
	/*
	int sum = 0;
	for (auto s : solution) {
		sum += boost::math::iround<double>(s * lcm);
	}
	*/
	if (lcm) {
		//std::cout << " Score: " << lcm << " for S: ";
		//for (int i = 0; i < denoms.size(); ++i)
		//	std::cout << nums[i] << "/" << denoms[i] << ", ";
		//std::cout << std::endl;
	}

	//std::cout << "Score: " << lcm << std::endl << std::endl;

	return lcm;
}

void makecsv(ReagentDB rdb, std::map<std::string, std::vector<double>> solmap, std::string ies, bool use_input_cache = true, std::string cmd_input_stoics = "", std::string cmd_input_precursors = "") {
	//1st row, reagents
	//every other row, solution
	auto start_time = std::chrono::high_resolution_clock::now();

	ofstream g(ies + ".csv");
	
	g << "Key,";
	std::string header = "";
	for (auto r : rdb()) {
		header += r.str() + ",";
	}
	header.erase(header.size() - 1);
	g << header;
	g << ",Score";
	g << std::endl;

	for (auto s : solmap) {
		g << s.first << ",";
		for (int i = 0; i < s.second.size() - 1; ++i) {
			g << s.second[i] << ",";
		}
		g << s.second[s.second.size() - 1];
		g << "," << score(s.second);
		g << std::endl;
	}
	g.close();

	if (use_input_cache == false) {
		ofstream icg("last_input_cache.txt");
		icg << "Reserved" << endl;
		icg << cmd_input_stoics << endl;
		icg << cmd_input_precursors << endl;
	}

	auto current_time = std::chrono::high_resolution_clock::now();
	auto clockval = std::chrono::duration<double>(current_time - start_time).count();
	std::cout << endl << ">> CSV operation completed in: " << clockval << "s" << endl;
} 

std::pair<MatrixXd, VectorXd> handle_input(StateDB sdb, ReagentDB rdb, std::vector<std::string> els, std::vector<double> amounts) {
	const auto coefficients = sdb.populate(els);
	const int nv = rdb().size();

	const char var = 'A';
	std::vector<std::string> vn;
	for (int i = 0; i < nv; ++i) {
		vn.push_back(std::string(1, var + i));
	}

	//std::vector<double> solution = { 1, 1, 1, 1 };
	std::vector<std::vector<std::pair<string, double>>> lhs, rhs;

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
		rhs.push_back({ {"RHS", amounts[i]} });
	}

	//std::cout << "Input: ";
	//for (auto via : amounts) {
	//	std::cout << via << " ";
	//}
	//std::cout << endl;

	MatrixXd A(lhs.size(), lhs[0].size());
	VectorXd B(lhs.size());

	for (int i = 0; i < els.size(); ++i) {
		for (int j = 0; j < nv; ++j) {
			A(i, j) = lhs[i][j].second;
		}

		B(i) = amounts[i];
	}

	return { A, B };
}

double dround(double x, int p = 3) {
	std::stringstream ss;
	ss << std::setprecision(p) << x;
	double meme;
	ss >> meme;
	return meme;
}

double validate_as_rational_approx(double x, int precision = 3, int maxdenum = 99) {
	auto ra = new_ra(x, maxdenum);
	double y = ((double)ra.first / (double)ra.second);

	return dround(y, precision);
}

std::string sol_as_string(std::vector<double> solution) {
	std::stringstream ss;
	int nz = solution.size();
	for (auto s : solution) {
		if (s == 0)
			nz--;
	}
	ss << nz << "!";
	for (int i = 0; i < solution.size() - 1; ++i) {
		ss << solution[i] << "#";
	}
	ss << solution[solution.size() - 1];

	return ss.str();
}

int main(int argc, char **argv) {
	// start program
	std::chrono::high_resolution_clock::time_point start_time, current_time;
	start_time = std::chrono::high_resolution_clock::now();

	double margin = 0.001;

	/*
	EXPERIMENTAL, use map to store solutions, with key being a string of values
	RES
	RES
	*/

	std::map<std::string, std::vector<double>> solmap;



	//config file parsing
	std::string confline;
	ifstream fconfig("config.cfg");
	std::vector<std::pair<std::string, double>> confvals;
	while (std::getline(fconfig, confline)) {
		std::vector<std::string> newstrs;
		boost::split(newstrs, confline, boost::is_any_of("="));
		double val;
		std::stringstream ss;
		ss << newstrs[1];
		ss >> val;
		confvals.push_back({newstrs[0], val});
	}

	for (auto& cval : confvals) {
		if (cval.first == "margin") {
			margin = cval.second;
		}
	}

	// command options
	bool use_input_cache = false;
	bool recache_margin_weights = false;
	bool csv = true;
	std::string cmd_input_stoics, cmd_input_precursors;
	try {
		boost::program_options::options_description desc("Allowed options");
		desc.add_options()
			("help", "Show command and syntax usage")
			("stoichs", boost::program_options::value<std::string>(),"Stoichiometry of the desired elements")
			("precursors", boost::program_options::value<std::string>(), "Present precursors to be used for calculations")
			("cache", boost::program_options::value<bool>(), "Wheter to use cached input or not")
			("margin", boost::program_options::value<double>(), "Margin used for precision and weight generation")
			("csv", boost::program_options::value<bool>(), "Wheter to use SQL output or CSV")
			("debug", "Run with debug flags on");
		boost::program_options::variables_map vm;
		boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
		boost::program_options::notify(vm);

		if (vm.count("help")) {
			cout << desc << endl;
			return 0;
		}

		// use input file then
		if (vm.count("cache")) {
			use_input_cache = vm["cache"].as<bool>();
			cout << "Using cache was set to " << vm["cache"].as<bool>() << endl;
		}

		if (vm.count("margin")) {
			if (margin != vm["margin"].as<double>()) {
				margin = vm["margin"].as<double>();
				recache_margin_weights = true;
			}
		}
		if (vm.count("csv")) {
			csv = vm["csv"].as<bool>();
			cout << "Using csv was set to " << vm["csv"].as<bool>() << endl;
		}

		if ((!vm.count("stoichs") || !vm.count("precursors")) && use_input_cache == false) {
			cout << "Attempting to execute with no input, exiting..." << endl;
			return -1;
		} else if (vm.count("stoichs") && vm.count("precursors")) {
			cmd_input_stoics = vm["stoichs"].as<std::string>();
			cmd_input_precursors = vm["precursors"].as<std::string>();
		}
	} catch(exception& e) {
		cerr << "Error: " << e.what() << endl;
	} catch(...) {
		cerr << "Exception of unknown type!" << endl;
	}
	// now parse the input and make it usable
	std::cout << "Recache margin weights was set to " << recache_margin_weights << endl;
	// end command options
	const ElementDB edb;
	ElementDB::init();

	StateDB sdb;

	std::cout.precision(4);
	std::cout.setf(ios::fixed);

	std::vector<std::string> els;
	std::vector<double> amounts;
	ReagentDB rdb;

	if (use_input_cache == true) {
		//const Input in(INPUT);
		const Input in("last_input_cache.txt");
		els = in().els;
		amounts = in().amounts;
		rdb = in().rdb;
	} else {
		auto parin = parse_input(cmd_input_stoics, cmd_input_precursors);
		els = parin.first.first;
		amounts = parin.first.second;
		rdb = parin.second;
	}


	const auto coefficients = sdb.populate(els);
	const int nv = rdb().size();

	std::cout.precision(3);
	std::cout.setf(ios::fixed);

	bool changed_factor = false;

	std::string ies = "";
	for (auto &e : els) {
		ies += e;
	}

	const char var = 'A';
	std::vector<std::string> vn;
	for (int i = 0; i < nv; ++i) {
		vn.push_back(std::string(1, var + i));
	}

	//std::vector<double> solution = { 1, 1, 1, 1 };
	std::vector<std::vector<std::pair<string, double>>> lhs, rhs;

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
		rhs.push_back({ {"RHS", amounts[i]} });
	}

	//std::cout << "Input: ";
	//for (auto via : amounts) {
	//	std::cout << via << " ";
	//}
	//std::cout << endl;

	MatrixXd A(lhs.size(), lhs[0].size());
	VectorXd B(lhs.size());

	for (int i = 0; i < els.size(); ++i) {
		for (int j = 0; j < nv; ++j) {
			A(i, j) = lhs[i][j].second;
		}

		B(i) = amounts[i];
	}

	//for (int i = 0; i < nv; ++i) {
	//	A(els.size(), i) = 1;
	//	B(els.size()) = 1;
	//}

	// A << 2, 2, 0, 0,
	// 	0, 0, 2, 2,
	// 	0, 1, 0, 3;// ,
	// 	 //1, 1, 1, 1;
	// B << 1, 1, 1;// , 1;

	FullPivLU<MatrixXd> lu;

	//cout << "A is: " << endl << A << endl;
	//cout << "B is: " << endl << B << endl;

	lu.compute(A);
	auto K = lu.kernel();
	VectorXd S = lu.solve(B);

	//cout << "Rank: " << lu.rank() << endl;
	//cout << "Solution:\n" << S << endl;
	//cout << "Kernel is:\n" << K << endl;

	// sol of points is given by P = S + t * K where t is a real number
	std::vector<double> inset;
	for (double i = 0; i < K.cols(); ++i) {
		inset.push_back(i);
	}

	// each kernel colum = column in array
	// each line in array = solution
	// issue is size as it becomes num_sols^n = big
	// maybe do a map? no as we'd need to also add a hashing function, but it wouldnt be too costy? mhm

	// Method 1 - Array
	vector<vector<double>> coeffs;

	// TODO: make it dynamic
	std::vector<std::vector<double>> prk;

	if (boost::filesystem::exists("combination_weights_cache.txt")) {
		ifstream cwcif("combination_weights_cache.txt");
		if (cwcif.peek() == cwcif.eof()) {
			cwcif.close();
			recache_margin_weights = true;
			std::cout << "Found empty margin weights file" << endl;
		}
	} else {
		recache_margin_weights = true;
		std::cout << "Couldn't find margin weights file" << endl;
	}
	if (recache_margin_weights == true) {
		std::cout << "Recaching margin weights... " << endl;
		std::vector<double> fvals;
		for (double i = 0; i <= 1 + margin/ 2; i += margin) {
			fvals.push_back(i);
		}

		//auto nprk = utils::permutation_mapk(fvals, K.cols());
		prk = utils::combination_k(fvals, K.cols());
		//prk = utils::permutation_mapk(fvals, K.cols());
		ofstream g("combination_weights_cache.txt");
		for (auto& cv : prk) {
			for (auto cvv : cv) {
				g << cvv << " ";
			}
			for (int i = cv.size(); i >= 0; --i) {
				g << cv[i] << " ";
			}
			g << "\n";
		}
		g.close();
		std::cout << "Margin weights recached." << endl;

		// adjust config values;
		for (int i = 0; i < confvals.size(); ++i) {
			if (confvals[i].first == "margin") {
				std::fstream ncfg("config.cfg", fstream::in | fstream::out | fstream::binary);
				int cl = 0;
				while (cl < i) {
					ncfg.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
					++cl;
				}

				ncfg.seekp(ncfg.tellp());
				ncfg << "margin=" << margin;
				ncfg.close();

				break;
			}
		}
	}

	ifstream f("combination_weights_cache.txt");
	double n;// , m;
	do {
		if (f.eof()) {
			break;
		}

		std::vector<double> tempvec;
		for (int i = 0; i < K.cols(); ++i) {
			f >> n;
			tempvec.push_back(n);
		}
		prk.push_back(tempvec);
	} while (f.good());
	std::cout << "Loaded weights successfully" << endl;

	int sols = 0;

	// now we have all possible combination weights given the factor
	// so we just multiply the column vectors by their values
	// TODO: multhithread this at program startup as it takes the most time
	// TODO: multhithreading sux, cache it instead
	//auto prk = utils::combination_k(fvals, 2);
	std::vector<std::vector<double>> solutions;
	std::vector<double> first_sol;
	
	for (int i = 0; i < lhs[0].size(); ++i) {
		first_sol.push_back(validate_as_rational_approx(S(i)));
	}
	//std::cout << sol_as_string(first_sol);
	solmap[sol_as_string(first_sol)] = first_sol;
	std::cout << std::endl;
	solutions.push_back(first_sol);

	//auto nS = S;
	Eigen::VectorXd nS(S.size());
	for (int i = 0; i < S.size(); ++i) {
		nS(i) = dround(S(i), utils::ndigits(1.f / margin) / 10);
		
		//solmap
	}

	for (auto &ps : prk) {
		Eigen::VectorXd p(S.size());
		//auto p = nS;
		for (int i = 0; i < S.size(); ++i) {
			p(i) = nS(i);
		}
		for (int i = 0; i < K.cols(); ++i) {
			p += ps[i] * K.col(i);
		}

		//check for negative solutions
		bool ffs = false;
		for (int i = 0; i < lhs[0].size(); ++i) {
			if (p(i) < 0)
				ffs = true;
		}
		if (ffs == true)
			continue;

		//double sum = 0;
		std::vector<double> solution;
		for (int j = 0; j < lhs[0].size(); ++j) {
			//cout << p(j) << " ";
			solution.push_back(validate_as_rational_approx(p(j)));
			//sum += rdb()[j].mass() * p(j);
		}

		solutions.push_back(solution);
		solmap[sol_as_string(solution)] = solution;

		//cout << "M(g/mol): " << sum;
		//std::cout << std::endl;
	}

	std::cout << "Total solutions: " << solutions.size() << endl;
	std::cout << "Unique solutions: " << solmap.size() << endl;
	std::vector<int> ranks;
	for (auto s : solutions) {
		ranks.push_back(score(s));
	}
	std::sort(ranks.begin(), ranks.end());
	std::cout << "Best score is: " << ranks[0] << std::endl;

	current_time = std::chrono::high_resolution_clock::now();
	auto clockval = std::chrono::duration<double>(current_time - start_time).count();
	std::cout << ">> Points computed in: " << clockval << "s" << endl << endl;
	
	//sqlite
	//makedb(rdb, solutions, ies, use_input_cache, cmd_input_stoics, cmd_input_precursors);
	//csv
	makecsv(rdb, solmap, ies, use_input_cache, cmd_input_stoics, cmd_input_precursors);
	//makecsv(rdb, solutions, ies, use_input_cache, cmd_input_stoics, cmd_input_precursors);

	return 0;
}