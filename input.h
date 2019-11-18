#pragma once

#include "defines.h"

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>

#include "reagentdb.h"
#include "utils.h"

#define INPUTPATH "input.txt"
#define CONFIGPATH "config.cfg"
#define INPUTCACHEPATH "licache.txt"
#define WEIGHTSCACHEPATH "cwcache.txt"

class Input
{
public:
	struct IOdata
	{
		std::string cmd_input_stoichs, cmd_input_precursors;
		std::vector<std::vector<double>> weights;
		std::vector<std::string> elements;
		std::vector<double> amounts;
		// 0 for stoichs, 1 for precursors
		int mode = 1;
		Reagent r;
		ReagentDB rdb;
		double margin = 1.f / 1000;
		int samples = 20;
		struct {
			bool use_input_cache = false;
			bool recache_margin_weights = false;
			bool csv = true;
		} options;
	};
private:
	IOdata _self;

public:
	Input(int argc, char** argv);
	Input(std::string f);

	void parse(std::string f);
	void parse(std::string stoichs, std::string precursors);

	//template<typename T>
	//std::vector<std::vector<T>> load(std::string f);
	std::vector<std::vector<double>> load(std::string f);

	void validate_weights(int nulcols);
	
	std::vector<std::pair<Eigen::MatrixXd, Eigen::VectorXd>> matrix();

	const IOdata& operator()() const;
};
