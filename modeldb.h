#pragma once

#include "defines.h"

#include "statedb.h"
#include "reagentdb.h"

#include "gurobi_c++.h"
#include <unordered_map>
#include "elementdb.h"

// NOTE: try the saving of models, see if it uses too much memory

class ElementDB;

struct ConstrData
{
	GRBLinExpr lhs = 0, rhs = 0;
	char sense = '=';
	//double rhs = 0;
	std::string name = "";
};

// TODO: clean this whole mess up
struct ModelData
{
	struct Params
	{
		std::string name;
		std::string value;
	};
	std::vector<Params> params;
	std::vector<std::vector<std::pair<std::string, double>>> lhs;
	std::vector<std::vector<std::pair<std::string, double>>> rhs;
	std::vector<char> sense;
	std::vector<ConstrData> constraints;
	std::unordered_map<std::string, GRBVar> vars;
	std::vector<double> coefficients;
	std::vector<std::string> col_names;

	std::string label = "";
	Reagent reagent;

	std::vector<GRBConstr> constrs;
	struct
	{
		std::pair<double, double> var_range = { 0, 20 };
		char type = 'I';
		int num_solutions = 100;
		bool repeat = false;
	} options;
};

class Model
{
private:
	GRBModel* _self;
	ModelData _data;

	friend class ModelDB;
public:
	Model(const GRBEnv& env, const ModelData& md) {
		try {
			_self = new GRBModel(env);
			_data = md;

			// parameters
			//for (auto& p : _data.params) {
			//	_self->set(p.name, p.value);
			//}

			_self->set(GRB_IntParam_OutputFlag, 0);
			_self->set(GRB_IntParam_InfUnbdInfo, 1);
			_self->set(GRB_IntParam_PoolSearchMode, 2);
			_self->set(GRB_IntParam_PoolSolutions, 100);

			// variables & coeff
			// first part, just create empty vars
			for (auto& n : _data.lhs) {
				for (auto& m : n) {
					if (m.first.empty())
						continue;

					if (_data.vars.find(m.first) == _data.vars.end())
						_data.vars[m.first] = _self->addVar(_data.options.var_range.first, _data.options.var_range.second, 1, _data.options.type, m.first);
					//_data.coefficients.push_back(m.second);
				}
			}
			for (auto& n : _data.rhs) {
				for (auto& m : n) {
					if (m.first.empty())
						continue;

					if (_data.vars.find(m.first) == _data.vars.end())
						_data.vars[m.first] = _self->addVar(_data.options.var_range.first, _data.options.var_range.second, 1, _data.options.type, m.first);
					//_data.coefficients.push_back(m.second);
				}
			}

			// expressions
			for (int i = 0; i < _data.lhs.size(); i++) {
				std::vector<GRBVar> lhs_vars, rhs_vars;
				std::vector<double> lhs_coeffs, rhs_coeffs;
				GRBLinExpr lhs = 0, rhs = 0;

				// do lhs expr
				for (auto &n : _data.lhs[i]) {
					if (n.first.empty()) {
						lhs += n.second;
						continue;
					}
						

					lhs_vars.push_back(_data.vars[n.first]);
					lhs_coeffs.push_back(n.second);
				}

				if (!lhs_vars.empty())
					lhs.addTerms(&lhs_coeffs[0], &lhs_vars[0], lhs_vars.size());
				// do rhs expr
				for (auto &n : _data.rhs[i]) {
					if (n.first.empty()) {
						rhs += n.second;
						continue;
					}
						

					rhs_vars.push_back(_data.vars[n.first]);
					rhs_coeffs.push_back(n.second);
				}

				if (!rhs_vars.empty())
					rhs.addTerms(&rhs_coeffs[0], &rhs_vars[0], rhs_vars.size());

				ConstrData cd;
				cd.lhs = lhs;
				cd.rhs = rhs;
				cd.sense = _data.sense[i];

				_data.constraints.push_back(cd);
			}

			// constraints
			for (auto& c : _data.constraints) {
				_self->addConstr(c.lhs, c.sense, c.rhs, c.name);
			}

		} 
		catch (GRBException e) {
			std::cout << "Error code = " << e.getErrorCode() << std::endl;
			std::cout << e.getMessage() << std::endl;
		}
		catch (...)
		{
			std::cout << "Exception during optimization" << std::endl;
		}
	}
	~Model() {
		_self->terminate();
	}
	const Model* operator()() const {
		return this;
	}

	// vector of variables and coefficients as lhs(members)
	void addConstr(GRBLinExpr lhs, GRBLinExpr rhs, char sense = '=') const {
		_self->addConstr(lhs, sense, rhs);
	}

	// compare solutions
	static int compare(std::vector<double> s1, std::vector<double> s2) {
		int l = 0, g = 0, e = 0;
		for (int i = 0; i < s1.size(); ++i) {
			if (s1[i] > s2[i]) ++g;
			if (s1[i] < s2[i]) ++l;
			if (s1[i] == s2[i]) ++e;
		}

		if (l == s1.size())
			return -1;
		if (g == s1.size())
			return -1;
		if (e == s1.size())
			return 0;

		return INT32_MAX;
	}

	// validate signs
	bool validate_state_sign_i() {
		return true;
	}

	// solve
	std::vector<std::vector<double>> solve() {
		_self->optimize();

		auto status = _self->get(GRB_IntAttr_Status);
		if (status != GRB_OPTIMAL) {
			/*std::cout << "The model is infeasible; computing IIS" << std::endl << std::endl;

			_self->computeIIS();
			std::cout << "The following constraints cannot be satisfied" << std::endl;
			auto c = _self->getConstrs();
			for (int i = 0; i < _self->get(GRB_IntAttr_NumConstrs); ++i)
			{
				if (c[i].get(GRB_IntAttr_IISConstr) == 1)
				{
					std::cout << c[i].get(GRB_StringAttr_ConstrName) << std::endl;
				}
			}
			*/
			return {};
		}

		std::vector<std::vector<double>> temp;

		// optimal solution
		std::vector<double> s;
		for (auto& v : _data.vars) {
			s.push_back(v.second.get(GRB_DoubleAttr_X));
		}
		temp.push_back(s);

		int sc = _self->get(GRB_IntAttr_SolCount);
		if (_data.options.num_solutions > 1 && sc > 1 ) {
			for (int i = 0; i < sc; ++i) {
				std::vector<double> iv;
				_self->set(GRB_IntParam_SolutionNumber, i);
				auto nv = _self->getVars();

				for (int j = 0; j < _data.vars.size(); ++j) {
					iv.push_back(nv[j].get(GRB_DoubleAttr_Xn));
				}

				bool search = false;
				for (int q = 0; q < temp.size(); ++q) {
					if (compare(temp[q], iv) == 0)
						search = true;
				}

				if (!search) {
					temp.push_back(iv);
				}
			}
		}

		//_self->terminate();

		return temp;
	}
};

class ModelDB
{
	GRBEnv _self;
	ElementDB _edb;
	std::vector<Model*> _models;

public:
	struct
	{
		bool use_edb = false;
		bool repeat = false;
		double tolerance = 0.001;
	} property;

	ModelDB();
	ModelDB(const ElementDB& edb) {
		_edb = edb;
	}

	static int validate_state_sign_i(std::vector<int> state);

	void insert(Model* m);
	void insert(std::vector<std::vector<std::pair<std::string, double>>> lhs,
		std::vector<std::vector<std::pair<std::string, double>>> rhs, std::vector<char> sense, char type = 'I',
		std::pair<double, double> range = { 0, 20 }, std::vector<std::string> col_names = {}, Reagent re = Reagent());
	void erase(int idx);
	void clear() {
		for (auto &m : _models) {
			delete m;
		}

		_models.clear();
	}

	std::vector<std::vector<std::vector<double>>> solve(ReagentDB rdb, std::ostream& out = std::cout, bool validate = false);

	// TODO: think
	// so data-> model_data needs ConstrData which needs GRBVar to be init
};
