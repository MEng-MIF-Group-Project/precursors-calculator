#pragma once

// STD Includes
#include <iostream>
#include <iomanip>
#include <vector>

// Eigen Includes
#include <Eigen/Dense>

#include <boost/math/common_factor.hpp>

inline Eigen::VectorXd lu_reduce(Eigen::MatrixXd A, Eigen::VectorXd B) {
	Eigen::FullPivLU<Eigen::MatrixXd> lu(A);

#ifdef DEBUG
	std::cout << "Matrix LU:" << std::endl;
	std::cout << lu.matrixLU() << "\n\n";
	
	MatrixXd L = MatrixXd::Identity(eq, eq);
	L.triangularView<StrictlyLower>() = lu.matrixLU();
	cout << "Matrix L:" << endl;
	cout << L << "\n\n";

	MatrixXd U = lu.matrixLU().triangularView<Upper>();
	cout << "Matrix U:" << endl;
	cout << U << "\n\n";

	MatrixXd iP = lu.permutationP().inverse();
	cout << "Inverse PermutationP:" << endl;
	cout << iP << "\n\n";

	MatrixXd iQ = lu.permutationQ().inverse();
	cout << "Inverse PermutationQ:" << endl;
	cout << iQ << "\n\n";

	MatrixXd R = iP * L * U * iQ;
	cout << "Reconstructed iP * L * U * iQ:" << endl;
	cout << R << "\n\n";

	cout << "Matrix rank:" << endl;
	cout << lu.rank() << "\n\n";
#endif

	Eigen::VectorXd S = lu.solve(B);
	if ((A * S).isApprox(B)) {
		std::cout << "Solution to AX = B:" << std::endl;
		std::cout << S << std::endl;
	}
	else {
		std::cout << "Solution to AX = B does not exist" << "\n\n";
	}

#ifdef DEBUG
	std::cout << "Residual: " << (A * S - B).norm() << "\n\n";

	Eigen::MatrixXd K = lu.kernel();
	std::cout << "Kernel: " << std::endl;
	std::cout << " Dim: " << lu.dimensionOfKernel() << std::endl;
	std::cout << K << "\n\n";

	for (auto i = 0; i < lu.dimensionOfKernel(); i++) {
		K.col(i) = K.col(i) + S;
	}

	std::cout << "Solution space S: " << std::endl;
	std::cout << K << "\n\n";
#endif

#ifdef DEBUG
	Eigen::MatrixXd C = (A * A.transpose())/static_cast<double>(A.cols() - 1);
	std::cout << "Sample covariance matrix C:" << std::endl;
	std::cout << C << "\n\n";
#endif

	return S;
}

inline Eigen::MatrixXd llt_reduce(Eigen::MatrixXd A, Eigen::VectorXd B) {
	Eigen::LLT<Eigen::MatrixXd> llt;
	llt.compute(A);

	Eigen::MatrixXd S = llt.solve(B);

	std::cout << "The solution is: \n" << S << std::endl;
	
	return S;
}