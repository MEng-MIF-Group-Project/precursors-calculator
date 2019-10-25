#include "modeldb.h"
#include "utils.h"
#include <boost/math/common_factor_ct.hpp>
#include <boost/math/tools/precision.hpp>
#include <boost/math/constants/constants.hpp>

ModelDB::ModelDB() = default;

int ModelDB::validate_state_sign_i(std::vector<int> state) {
	int ap = 0, an = 0;
	for (auto i : state) {
		if (i >= 0) ++ap;
		if (i <= 0) ++an;
	}

	if (ap == state.size())
		return 1;
	if (an == state.size())
		return -1;
	
	return 0;
}

void ModelDB::insert(Model* m) {
	_models.push_back(m);
}

void ModelDB::insert(std::vector<std::vector<std::pair<std::string, double>>> lhs,
                     std::vector<std::vector<std::pair<std::string, double>>> rhs, std::vector<char> sense, char type,
                     std::pair<double, double> range, std::vector<std::string> col_names, Reagent re) {
	//params setup
	ModelData md;
	// options phase
	md.options.type = type;
	md.options.var_range = range;
	// default_params
	md.params.push_back({"GRB_IntParam_OutputFlag", "0"});
	md.params.push_back({"GRB_IntParam_InfUnbdInfo", "1"});
	md.params.push_back({"GRB_IntParam_PoolSearchMode", "2"});
	md.params.push_back({"GRB_IntParam_PoolSolutions", "100"});
	//vars first from elements
	//md.terms = terms;
	md.lhs = lhs;
	md.rhs = rhs;
	md.sense = sense;
	// misc stuff

	md.reagent = re;
	md.label = re.str();
	md.col_names = col_names;
	//now create expressions based on oxidation state
	_models.push_back(new Model(_self, md));
}

void ModelDB::erase(int idx) {
	_models.erase(_models.begin() + idx);
}

std::vector<std::vector<std::vector<double>>> ModelDB::solve(ReagentDB rdb, std::ostream& os, bool validate)
{
	std::vector<std::vector<std::vector<double>>> tmodels;
	tmodels.reserve(_models.size());

	std::map<std::string, std::vector<double>> scores;
	std::vector<std::string> streams;
	std::vector<double> scores_vals;

	for (auto &m : _models) {
		// TODO: do somethign with this
		if (!m->validate_state_sign_i())
			continue;

		std::vector<std::string> var_name;
		for (auto& vm : m->_data.vars)
			var_name.push_back(vm.first);
		// only works if var_name == element_name
		auto rq = rdb.quanitfy(var_name);

		auto sv = m->solve();

		if (sv.empty())
			continue;

		std::vector<double> temp_scores;
		double temp_mscore = 0;
		std::string temp_key;
		double temp_tmas = 0;
		for (int i = 0; i < m->_data.col_names.size(); ++i) {
			double mass = 0;
			for (int j = 0; j < rdb().size(); ++j) {
				Reagent r = rdb()[j];
				std::string sr = r.str();
				if (sr == m->_data.col_names[i]) {
					mass = rdb()[j].mass() * sv[0][i];
					temp_tmas += mass;
					double tes = (rdb().size() - j) * mass;// *vi[i];
					temp_scores.push_back(tes);
					temp_mscore += tes;
					//vs.push_back(tes);
					temp_key += std::to_string(j);
					break;
				}
			}
		}
		//out << temp_mscore << " " << boost::math::constants::euler<double>() << " = ";
		temp_mscore *= pow(boost::math::constants::euler<double>(), static_cast<double>(temp_scores.size()));
		//out << temp_mscore << std::endl;

		//deviation
		double dev = std::abs(m->_data.reagent.mass() - temp_tmas);
		double dev_ratio = (dev / m->_data.reagent.mass()) * 100;
		temp_mscore /= (dev_ratio == 0 ? 1 : (dev_ratio < 1 ? 1/dev : dev));

		temp_scores.push_back(temp_mscore);

		auto old_entry = scores.find(temp_key);
		if (old_entry != scores.end()) {
			bool verifier = true;
			for (int i = 0; i < old_entry->second.size(); ++i) {
				if (old_entry->second[i] != temp_scores[i]) {
					verifier = false;
				}
			}
			if (verifier) {
				continue;
			}
		}

		scores[temp_key] = temp_scores;

		std::stringstream out;
		scores_vals.push_back(temp_mscore);
		//out << "Optimising state: ";
		out << m->_data.label << " ";
		out << "(" << m->_data.reagent.mass() << " g/mol)";
		out << std::endl;
		out << std::fixed;
		out << "> Score: " << temp_mscore << std::endl;
		out << "> Deviation: " << dev_ratio << "%" << std::endl;
		/*
		for (auto c : m->_data.col_names) {
			out << std::setw(10-c.size()) << " ";
			out << c;
			out << "\t ";
		}
		out << std::endl;
		
		for (auto l = 0; l < m->_data.lhs.size(); ++l) {
			for (auto r : m->_data.lhs[l]) {
				//std::cout << r.second << " ";
				utils::ncout(r.second, out);
			}
			//std::cout << " | ";
			out << std::setw(4) << " ";
			out << m->_data.sense[l];
			//std::cout << std::setw(4) << " ";

			for (auto r : m->_data.rhs[l]) {
				//std::cout << r.second << " ";
				utils::ncout(r.second, out);
			}

			out << std::endl;
		}
		//std::cout << " >> " << std::endl;
		out << std::endl;
		*/
		std::vector<std::vector<double>> tsols;
		//out << "Num Solutions: " << sv.size() << std::endl;
		int si = 1;
		//tsols.reserve(sv.size());
		for (auto vi : sv) {
			out << "> Reagents: ";// << si;
			for (auto& tcn : m->_data.col_names) {
				out << std::setw(10) << " ";
				out << tcn;
			}
			out << std::endl;
			out << "> Solution: ";// << si;
			out << std::setw(5 - utils::ndigits(si)) << " ";
			int spacing_i = 0;
			std::vector<double> tvals;
			int ti1 = 0;
			for (auto vei : vi) {
				//std::cout << vei << " ";
				if (vei == 0.0)
					vei = 0;
				utils::ncout(vei, out, 8);
				//std::cout << " ";
				if (vei < 10)
					spacing_i += 2;
				else if (vei >= 10 && vei < 100)
					spacing_i += 3;

				tvals.push_back(vei);
			}


			/*
			if (validate) {
				auto rhs = Reagent(var_name, vi);
				auto res = rdb.validate(rhs);//validate_reagents(lhs, rhs);
				out << std::setw(10 - spacing_i) << " ";
				if (res == 0) {
					out << "validated";
				}
				else if (res < 0) {
					out << "not validated, not enough " << var_name[std::abs(res + 1)];
					out << std::setw(3 - var_name[std::abs(res + 1)].size()) << " ";
					out << rq[std::abs(res + 1)] << "/" << vi[std::abs(res + 1)] << " atoms needed";
				}
				else if (res > 0) {
					out << "not validated, missing " << var_name[std::abs(res - 1)];
				}
			}*/
			std::vector<double> vs;
			if (property.use_edb) {
				if (!m->_data.label.empty()) {
					std::string key;
					double tm = 0;
					double score = 0;
					out << std::endl;
					out.precision(4);
					out << std::fixed;
					std::vector<double> mms;
					out << "> Mass (g/mol): ";
					for (int i = 0; i < vi.size(); ++i) {
						double mass = 0;
						for (int j = 0; j < rdb().size(); ++j) {
							Reagent r = rdb()[j];
							std::string sr = r.str();
							if(sr == m->_data.col_names[i]) {
								mass = rdb()[j].mass() * vi[i];
								double tes = (rdb().size() - j);// *vi[i];
								vs.push_back(tes);
								key += std::to_string(j);
								break;
							}
						}
						utils::ncout(mass, out, 8);
						mms.push_back(mass);
						tm += mass;
					}

					for (auto vsv : vs) {
						score += vsv;
					}
					score *= 1 / static_cast<float>(vi.size());
					vs.push_back(score);
					out << std::endl;
					out << "> Total Mass (g/mol): ";
					utils::ncout(tm, out, 8, "", true);

					//out << std::endl;
					
					out << "> Mass (% of total): ";
					double atms = 0;
					for (auto ms : mms) {
						double ams = 0;
						ams = ms / tm;
						out << ams*100 << "%\t\t";
						//utils::ncout(ams, out, 8);
						atms += ams;
					}
					//out << "> Score: " << score << std::endl;
				}				
			}

			out << std::endl;
			++si;
			tsols.push_back(tvals);
		}

		out << std::endl;
		tmodels.push_back(tsols);

		streams.push_back(out.str());
	}

	if (scores_vals.empty()) {
		os << ">> No solutions <<" << std::endl;
		
		return {};
	}

	std::vector<int> indices;
	indices.reserve(scores_vals.size());
	for (int i = 0; i < scores_vals.size(); i++) {
		indices.push_back(i);
	}
	for (int i = 0; i < scores_vals.size() - 1; i++) {
		for (int j = 0; j < scores_vals.size() - i - 1; j++) {
			if (scores_vals[j] < scores_vals[j+1]) {
				std::swap(scores_vals[j], scores_vals[j + 1]);
				std::swap(indices[j], indices[j + 1]);
			}
		}
	}
	for (int i = 0; i < indices.size(); i++) {
		os << streams[indices[i]];
	}

	return tmodels;
}
