#pragma once
#include <vector>

class Matrix
{
	std::vector<std::vector<double>> _self;

public:
	Matrix() = default;
	Matrix(int rows = 1, int cols = 1) {
		// = new std::vector<std::vector<double>>(rows);

		for (auto i = 0; i < rows; ++i) {
			std::vector<double> r;
			for (auto j = 0; j < cols; ++j) {
				r.push_back(0);
			}
			_self.push_back(r);
		}
	}

	const std::vector<std::vector<double>>& operator()() const {
		return _self;
	}

	const std::vector<double>& operator()(int row) const {
		return _self[row];
	}

	const double& operator()(int row, int col) const {
		return _self[row][col];
	}

	void insert_column(std::vector<double> col) {
		int i = 0;
		for (i = 0; i < col.size(); ++i) {
			_self[i].push_back(col[i]);
		}
		for (int j = i; j < _self.size(); ++j) {
			_self[j].push_back(0);
		}
	}

	void insert_row(const std::vector<double>& row) {
		_self.push_back(row);
	}
};
