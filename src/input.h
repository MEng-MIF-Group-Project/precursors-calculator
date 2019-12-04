#pragma once

#include "defines.h"

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string_regex.hpp>
#include <boost/math/special_functions/factorials.hpp>

#include "reagentdb.h"

#include "math.h"
#include "utils.h"

#define INPUTPATH "input.txt"
#define CONFIGPATH "config.cfg"
#define INPUTCACHEPATH "licache.txt"
#define WEIGHTSCACHEPATH "cwcache.txt"

//#define DEBUG_INPUT ON

class Input
{
public:
	struct IOdata
	{
		// member
		std::string cmd_input_stoichs, cmd_input_precursors;
		std::vector<std::vector<double>> weights;
		std::vector<std::string> elements;
		std::vector<double> amounts;
		// 0 for stoichs, 1 for precursors
		int mode = 0;
		Reagent r;
		ReagentDB rdb;

		// property
		double margin = 1.f / 1000;
		int samples = 20;
		int drange = 10; //
		double dmass = 0.5; // 0.5grams
		struct {
			bool use_input_cache = false;
			bool recache_margin_weights = false;
			bool csv = true;
		} options;

		// debug
		bool debug = false;
		bool mweights = false; // whether weights are mirrored, i.e 0.2 and 0.8 also output 0.8 and 0.2
	};
private:
	IOdata _self;

public:
	// initialise input arguments from cli
	Input(int argc, char** argv);
	// initialise input arguments from file, TODO: does this have a purpose?
	Input(std::string f);
	// load input data from a string
	void parse(std::string content);
	// load input data based on cli arguments
	void parse(std::string stoichs, std::string precursors);

	//template<typename T>
	//std::vector<std::vector<T>> load(std::string f);
	std::vector<std::vector<double>> load(std::string f);

	// either generates or loads the weights if they already exist on disk
	void validate_weights(int nulcols, int nulnegs = 0);
	
	// return Eigen based matrices based on input data
	std::vector<std::pair<Eigen::MatrixXd, Eigen::VectorXd>> matrix();

	// returns io data
	const IOdata& operator()() const;
};
