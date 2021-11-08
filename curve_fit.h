#ifndef CURVE_FIT_H
#define CURVE_FIT_H

#include <vector>
#include <assert.h>

#include <armadillo>

using arma::Col;
using arma::Mat;

template <typename T>
class CurveFit
{

private:
    Col<T> bezerror(const int deg, const Mat<T> &data, const Mat<T> &C, const Mat<T> &P);
    Mat<T> bezjac(const int deg, const Mat<T> &data, const Mat<T> &B, const Col<T> &w, const Col<T> &f, const bool &bounded);
    Mat<T> bern(const int deg, const Col<T> &t);
    T nchoosek(const int _n, const int _k);

public:
    CurveFit(){};
    Mat<T> bezier_fit(const Mat<T> &data, const int &deg, const bool &bounded);
};

#endif // CURVE_FIT_H
