#include "input.h"
#include "statedb.h"

//#include <random>

Input::Input(int argc, char** argv)
{
#ifndef DEBUG_INPUT
	try {
		boost::program_options::options_description desc("Allowed options");
		boost::program_options::variables_map vm;
		
		desc.add_options()
			("help", "Show command and syntax usage")
			("stoichs", boost::program_options::value<std::string>(), "Stoichiometry of the desired elements")
			("precursors", boost::program_options::value<std::string>(), "Present precursors to be used for calculations")
			("cache", boost::program_options::value<bool>(), "Wheter to use cached input or not")
			("margin", boost::program_options::value<double>(), "Margin used for precision and weight generation")
			("csv", boost::program_options::value<bool>(), "Wheter to use SQL output or CSV")
			("samples", boost::program_options::value<int>(), "Number of point samples to take")
			("mode", boost::program_options::value<int>(), "Switch between different calculators, 0 for stoichs, 1 for precursors")
			("dmass", boost::program_options::value<double>(), "Desired mass of solution")
			("drange", boost::program_options::value<int>(), "Denominator range")
			("mweights", boost::program_options::value<bool>(), "Turn on/off mirrored weights")
			("debug", boost::program_options::value<bool>(), "Run with debug flags on");

		boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
		boost::program_options::notify(vm);


		if (vm.count("help")) {
			std::cout << desc << std::endl;
			exit(0);
		}
		if (vm.count("samples")) {
			_self.samples = vm["samples"].as<int>();
			std::cout << "Number of samples set to " << _self.samples << std::endl; 
		}
		if (vm.count("dmass")) {
			_self.dmass = vm["dmass"].as<double>();
			std::cout << "Desired mass set to " << _self.dmass << std::endl;
		}
		if (vm.count("drange")) {
			_self.drange = vm["drange"].as<int>();
			std::cout << "drange set to " << _self.drange << std::endl;
		}

		if (vm.count("mode")) {
			_self.mode = vm["mode"].as<int>();
		}

		// use input file then
		if (vm.count("cache")) {
			_self.options.use_input_cache = vm["cache"].as<bool>();
			std::cout << "Using cache was set to " << vm["cache"].as<bool>() << std::endl;
		}

		if (vm.count("margin")) {
			_self.margin = vm["margin"].as<double>();
			_self.options.recache_margin_weights = true;
		}

		if (vm.count("csv")) {
			_self.options.csv = vm["csv"].as<bool>();
			std::cout << "Using csv was set to " << vm["csv"].as<bool>() << std::endl;
		}

		if (vm.count("mweights")) {
			_self.mweights = vm["mweights"].as<bool>();
			std::cout << "Mirrored weights was set to " << _self.mweights << std::endl;
		}

		if (vm.count("debug")) {
			_self.debug = vm["debug"].as<bool>();
			std::cout << "Debug mode was set to " << _self.debug << std::endl;
		}

		if (_self.mode == 0) {
			if ((!vm.count("stoichs")) && _self.options.use_input_cache == false) {
				std::cout << "Attempting to execute with no input, exiting..." << std::endl;
				exit(-1);
			}
			else if (vm.count("stoichs")) {
				_self.cmd_input_stoichs = vm["stoichs"].as<std::string>();
			}
		}
		else if (_self.mode == 1) {
			if ((!vm.count("stoichs") || !vm.count("precursors")) && _self.options.use_input_cache == false) {
				std::cout << "Attempting to execute with no input, exiting..." << std::endl;
				exit(-1);
			}
			else if (vm.count("stoichs") && vm.count("precursors")) {
				_self.cmd_input_stoichs = vm["stoichs"].as<std::string>();
				_self.cmd_input_precursors = vm["precursors"].as<std::string>();
			}
		}
	}
	catch (std::exception& e) {
		std::cerr << "Error: " << e.what() << std::endl;
	}
	catch (...) {
		std::cerr << "Exception of unknown type!" << std::endl;
	}
#else if
	_self.samples = 100;
	_self.cmd_input_stoichs = "LiAlSO";
	_self.cmd_input_precursors = "Li2S Al2S3 Al2O3 LiAlO2 Li2O";
	//_self.margin = 0.003;
	_self.dmass = 0.5;
	_self.drange = 20;
	_self.options.recache_margin_weights = true;
	//_self.options.use_input_cache = true;
	_self.mode = 1;
#endif
	std::cout << "Recache margin weights was set to " << _self.options.recache_margin_weights << std::endl;
	if (_self.mode == 0) {
		if (_self.options.use_input_cache == true) {
			parse(INPUTCACHEPATH);
		}
		else {
			parse(_self.cmd_input_stoichs, "");
		}
	}
	else if (_self.mode == 1) {
		if (_self.options.use_input_cache == true) {
			parse(INPUTCACHEPATH);
		}
		else {
			parse(_self.cmd_input_stoichs, _self.cmd_input_precursors);
		}
	}

	std::cout << "Finishing parsing cmline input" << std::endl;
}

Input::Input(std::string f) {
	parse(f);
}

void Input::parse(std::string f) {
	// Line split
	std::vector<std::string> ld;
	std::ifstream file(f);

	std::string line;
	while (std::getline(file, line)) {
		ld.push_back(line);
	}

	// Line Parsing
	// 0. Reserved
	// 1. List of Elements
	// 2. Available reagents
	boost::char_separator<char> s0(" ");
	boost::char_separator<char> s1(":");
	boost::tokenizer<boost::char_separator<char>> t1(ld[1], s0);

	// Retrieve name of elements
	for (const auto& t : t1) {
		boost::tokenizer<boost::char_separator<char>> t2_0(t, s1);
		for (const auto& m : t2_0) {
			int c = 0;
			const auto r = m.find('#');
			if (r == std::string::npos) {
				_self.elements.push_back(m);
				_self.amounts.push_back(1);
				_self.r.insert(Element(m, 1));
			}
			else {
				const std::string n = m.substr(0, r);

				std::stringstream ss;
				ss << m.substr(r + 1);
				ss >> c;

				_self.elements.push_back(n);
				_self.amounts.push_back(c);
				_self.r.insert(Element(n, c));
			}
		}
	}

	boost::tokenizer<boost::char_separator<char>> t2(ld[2], s0);
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
		_self.rdb.insert(reagent);
	}

	if (_self.mode == 0) {
		for (auto e : _self.r()) {
			Reagent re;
			re.insert(e);
			_self.rdb.insert(re);
		}
	}
}

void Input::parse(std::string stoichs, std::string precursors)
{
	boost::char_separator<char> s0(" ");
	boost::char_separator<char> s1(":");
	boost::tokenizer<boost::char_separator<char>> t1(stoichs, s0);
	
	{
		boost::regex re("[A-Z][a-z]?");
		boost::sregex_token_iterator ti(stoichs.begin(), stoichs.end(), re, { -1, 0 });
		std::vector<std::string> tokens;
		std::remove_copy_if(ti, boost::sregex_token_iterator(), std::back_inserter(tokens), [](std::string const &s) { return s.empty(); });
		int last_type = 2;
		for (auto& p : tokens) {
			if (std::isdigit(p[0])) {
				// change this
				std::stringstream ss(p);
				int amount;
				ss >> amount;

				_self.amounts.push_back(amount);
				last_type = 2;
			}
			else {
				if (last_type == 1) {
					_self.amounts.push_back(1);
				}
				_self.elements.push_back(p);
				last_type = 1;
			}
		}
		if (last_type == 1) {
			_self.amounts.push_back(1);
		}

		for (int i = 0; i < _self.elements.size(); ++i) {
			_self.r.insert(Element(_self.elements[i], _self.amounts[i]));
		}
	}
	/*
	// Retrieve name of elements
	for (const auto& t : t1) {
		boost::tokenizer<boost::char_separator<char>> t2_0(t, s1);
		for (const auto& m : t2_0) {
			int c = 0;
			const auto r = m.find('#');
			if (r == std::string::npos) {
				_self.elements.push_back(m);
				_self.amounts.push_back(1);
				_self.r.insert(Element(m, 1));
			}
			else {
				const std::string n = m.substr(0, r);

				std::stringstream ss;
				ss << m.substr(r + 1);
				ss >> c;

				_self.elements.push_back(n);
				_self.amounts.push_back(c);
				_self.r.insert(Element(n, c));
			}
		}
	}
	*/

	boost::tokenizer<boost::char_separator<char>> t2(precursors, s0);
	for (const auto& t : t2) {
		/*/
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
		_self.rdb.insert(reagent);
		*/
		Reagent reagent;
		boost::regex re("[A-Z][a-z]?");
		boost::sregex_token_iterator ti(t.begin(), t.end(), re, { -1, 0 });

		std::vector<std::string> tokens;
		std::remove_copy_if(ti, boost::sregex_token_iterator(), std::back_inserter(tokens), [](std::string const &s) { return s.empty(); });
		int last_type = 2;
		std::vector<std::string> elements;
		std::vector<int> amounts;
		for (auto& p : tokens) {
			if (std::isdigit(p[0])) {
				int amount;
				std::stringstream(p) >> amount;

				amounts.push_back(amount);
				last_type = 2;
			}
			else {
				if (last_type == 1) {
					amounts.push_back(1);
				}
				elements.push_back(p);
				last_type = 1;
			}
		}
		if (last_type == 1) {
			amounts.push_back(1);
		}
		for (int i = 0; i < elements.size(); ++i) { 
			reagent.insert(Element(elements[i], amounts[i]));
		}

		_self.rdb.insert(reagent);
	}

	if (_self.mode == 0) {
		for (auto e : _self.r()) {
			Reagent re;
			re.insert(e);
			_self.rdb.insert(re);
		}

		for (auto rs : _self.rdb()) {
			rs.str();
		}
		
		_self.r.str();
	}
}

void Input::validate_weights(int nulcols, int nulnegs) {
	std::ifstream fconfig(CONFIGPATH);
	std::string confline;
	std::vector<std::pair<std::string, double>> confvals;
	while (std::getline(fconfig, confline)) {
		std::vector<std::string> newstrs;
		boost::split(newstrs, confline, boost::is_any_of("="));
		double val;
		std::stringstream ss;
		ss << newstrs[1];
		ss >> val;
		confvals.push_back({ newstrs[0], val });
	}

	if (boost::filesystem::exists(WEIGHTSCACHEPATH)) {
		std::ifstream cwcif(WEIGHTSCACHEPATH);
		if (cwcif.peek() == static_cast<int>(cwcif.eof())) {
			cwcif.close();
			_self.options.recache_margin_weights = true;
			std::cout << "Found empty margin weights file" << std::endl;
		}
	}
	else {
		_self.options.recache_margin_weights = true;
		std::cout << "Couldn't find margin weights file, generating a new one" << std::endl;
	}

	//experimental - dynamic drange based on samples
	// so amount of weights = (drange+1)^2, which must be > nr of samples
	// we also need to account for negative vectors which might skew the results and 
	// leave us with just a tiny fraction of initial generated weights
	if (_self.debug == true) {
		int wsize = std::pow(_self.drange + 1, 2);
		int adjrange = _self.drange + 1;
		if (wsize < _self.samples) {
			do {
				wsize = std::pow(++adjrange, 2);
			} while (wsize < _self.samples);
			//--adjrange;
		}
		else {
			do {
				wsize = std::pow(--adjrange, 2);
			} while (wsize > _self.samples);
			++adjrange;
		}
		_self.drange = adjrange + nulnegs / nulcols;
		std::cout << "Adjusted drange was set to: " << _self.drange << std::endl;
	}

	// generate new wegihts if needed or load old ones
	if (_self.options.recache_margin_weights == true) {
		std::vector<double> fvals;
		fvals.push_back(1);
		for (int i = 2; i <= _self.drange; ++i) {
			fvals.push_back(1.f / i);
		}
		fvals.push_back(0);

		std::vector<std::vector<double>> i_am_speed;
		math::for_each_combination(fvals, nulcols, [&](std::vector<double> &gp) {
			i_am_speed.push_back(gp);
		});
		_self.weights = i_am_speed;
		
		if (_self.debug) {
			// (k+n-1)! / (k!(n-1)!) comb with rep
			double vf = static_cast<double>(boost::math::factorial<double>(nulcols + fvals.size() - 1))
				/ (boost::math::factorial<double>(nulcols) * boost::math::factorial<double>(fvals.size() - 1));

			std::cout << "Combinations of n: " << fvals.size() << ", k: " << nulcols <<" generated: " << i_am_speed.size() << std::endl;
			std::cout << "Expected (k+n-1)! / (k!(n-1)!) : " << vf << std::endl; 
			// turns out i_am_speed's size is eq to fvals.size() ^ 2 instead of comb formula
		}
		

		std::ofstream g(WEIGHTSCACHEPATH);
		for (auto& cv : _self.weights) {
			for (auto cvv : cv) {
				g << cvv << " ";
			}
			g << '\n';
			if (_self.mweights) {
				for (int i = static_cast<int>(cv.size()) - 1; i >= 0; --i) {
					g << cv[i] << " ";
				}
				g << '\n';
			}
		}
		g.close();

		// adjust config values;
		for (int i = 0; i < confvals.size(); ++i) {
			if (confvals[i].first == "margin") {
				std::fstream ncfg(CONFIGPATH, std::fstream::in | std::fstream::out | std::fstream::binary);
				int cl = 0;
				while (cl < i) {
					ncfg.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
					++cl;
				}

				ncfg.seekp(ncfg.tellp());
				ncfg << "margin=" << _self.margin;
				ncfg.close();

				break;
			}
		}
	}
	else {
		for (auto& cval : confvals) {
			if (cval.first == "margin") {
				_self.margin = cval.second;
				break;
			}
		}
	}
}

std::vector<std::pair<Eigen::MatrixXd, Eigen::VectorXd>> Input::matrix()
{
	if (_self.mode == 0) {
		std::vector<std::pair<Eigen::MatrixXd, Eigen::VectorXd>> temp;
		StateDB sdb;
		const auto coefficients = sdb.populate(_self.elements);

		for (auto c : coefficients) {
			Eigen::MatrixXd A(2, _self.r().size());
			Eigen::VectorXd B(2);
			Eigen::RowVectorXd rv1(_self.r().size()), rv2(_self.r().size());
			//std::cout << "States:\n";
			for (int j = 0; j < _self.r().size(); ++j) {
				rv1(j) = _self.amounts[j];
				rv2(j) = c[j];
				//std::cout << c[j] << " ";
			}
			//std::cout << std::endl << std::endl;

			A << rv1, rv2;
			B << 1, 0;

			//std::cout << std::endl << A << std::endl << std::endl;

			temp.push_back({ A, B });
		}

		//std::cout << "Done processing matrices" << std::endl; 

		return temp;
	}
	else if (_self.mode == 1){
		const char var = 'A';
		std::vector<std::string> vn;
		for (int i = 0; i < _self.rdb().size(); ++i) {
			vn.push_back(std::string(1, var + i));
		}

		std::vector<std::vector<std::pair<std::string, double>>> lhs, rhs;
		for (int i = 0; i < _self.r().size(); ++i) {
			std::vector<std::pair<std::string, double>> lhs_row;
			for (int j = 0; j < _self.rdb().size(); ++j) {
				double val = 0;
				for (auto &k : _self.rdb(j)()) {
					if (k().n == _self.r()[i]().n) {
						val = k().q;
						// TODO: think about this break
						break;
					}
				}

				lhs_row.push_back({ vn[j], val });
			}
			lhs.push_back(lhs_row);
			rhs.push_back({ {"RHS", _self.r()[i]().q} });
		}

		Eigen::MatrixXd A(lhs.size(), lhs[0].size());
		Eigen::VectorXd B(lhs.size());

		for (int i = 0; i < _self.r().size(); ++i) {
			for (int j = 0; j < _self.rdb().size(); ++j) {
				A(i, j) = lhs[i][j].second;
			}

			B(i) = _self.r()[i]().q;
		}

		return { { A, B } };
	}
	
	return {};
}

const Input::IOdata& Input::operator()() const{
	return _self;
}

//template<typename T>
//std::vector<std::vector<T>> Input::load(std::string f)
//{
//	std::vector<std::vector<T>> l;
//
//	std::string line;
//	std::ifstream in(f);
//	while (std::getline(f, line)) {
//		std::vector<T> nl;
//
//		std::stringstream ss(line);
//		T n;
//		while (ss >> n) {
//			nl.push_back(n);
//		}
//		l.push_back(nl);
//	}
//
//	return l;
//}

std::vector<std::vector<double>> Input::load(std::string f)
{
	std::vector<std::vector<double>> l;

	std::string line;
	std::ifstream in(f);
	while (std::getline(in, line)) {
		std::vector<double> nl;

		std::stringstream ss(line);
		double n;
		while (ss >> n) {
			nl.push_back(n);
		}
		l.push_back(nl);
	}

	return l;
}