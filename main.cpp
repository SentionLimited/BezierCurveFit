#include "curve_fit.h"
#include <armadillo>

using arma::Col;
using arma::Mat;

using TD = double;

int main(int argc, char *argv[])
{
	const int deg = 5;
	const bool bounded = true;

	// Get input data
	Mat<TD> data;
	data.load("test_data.dat", arma::raw_ascii);

	auto cur = CurveFit<TD>();
	auto bez = cur.bezier_fit(data, deg, bounded);

	bez.save("build/output/bezier_fit.dat", arma::raw_ascii);

	return 0;
}
