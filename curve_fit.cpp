#include "curve_fit.h"

using arma::Col;
using arma::Mat;

template <typename T>
T CurveFit<T>::nchoosek(const int _n, const int _k)
{

    T n = _n;
    T k = _k;
    const T zero = 0;
    const T one = 1;
    const T two = 2;

    if (k > n)
        return zero;
    if (k * two > n)
        k = n - k;
    if (k == zero)
        return one;

    T result = n;
    for (T i = two; i <= k; ++i)
    {
        result *= (n - i + one);
        result /= i;
    }

    return result;
}

template <typename T>
Mat<T> CurveFit<T>::bern(const int deg, const Col<T> &t)
{

    const int num_nodes = t.n_elem;

    assert(num_nodes > 0 && "There must be at least 1 output node from Bezier fit");

    Mat<T> B(num_nodes, deg + 1);

    for (int i = 0; i < (deg + 1); ++i)
    {
        B.col(i) = nchoosek(deg, i) * (arma::pow(t, i) % arma::pow(1. - t, deg - i));
    }

    return B;
}

template <typename T>
Col<T> CurveFit<T>::bezerror(const int deg, const Mat<T> &data, const Mat<T> &C, const Mat<T> &P)
{

    // Error from data
    Col<T> fd = arma::vectorise(C - data);

    // Error from control points
    Mat<T> dP = arma::diff(P);
    Col<T> fp = arma::vectorise(arma::diff(dP));

    return arma::join_cols(fd, fp);
}

template <typename T>
Mat<T> CurveFit<T>::bezjac(const int deg, const Mat<T> &data, const Mat<T> &B, const Col<T> &w, const Col<T> &f, const bool &bounded)
{

    int num_errs = f.n_elem;
    int num_vars = 2 * (deg + 1);

    T incr = 1e-5;

    Mat<T> J(num_errs, num_vars, arma::fill::zeros);

    for (int i = 0; i < num_vars; ++i)
    {

        if ((!bounded) ||
            (i != 0 &&
             i != deg &&
             i != (deg + 1) &&
             i != (num_vars - 1)))
        {

            Col<T> w_pert = w;
            w_pert(i) += incr;

            Mat<T> P_pert = arma::reshape(w_pert, deg + 1, 2);
            Mat<T> C_pert = B * P_pert;

            Col<T> df = bezerror(deg, data, C_pert, P_pert) - f;

            J.col(i) = df / incr;
        }
    }

    return J;
}

template <typename T>
Mat<T> CurveFit<T>::bezier_fit(
    const Mat<T> &data,
    const int &deg,
    const bool &bounded)
{

    int num_points = data.n_rows;

    // Initialize array for control points
    // Mat<T> P(deg + 1, 2, arma::fill::randu);
    Mat<T> P(deg + 1, 2, arma::fill::ones);

    // Set the first and last control points
    P.head_rows(1) = data.head_rows(1);
    P.tail_rows(1) = data.tail_rows(1);

    // Initialize t vector
    Col<T> t(num_points, arma::fill::zeros);
    Col<T> local_dists = arma::sqrt(arma::sum(arma::square(data.tail_rows(num_points - 1) - data.head_rows(num_points - 1)), 1));
    Col<T> cum_local_dists = arma::cumsum(local_dists);
    t.tail(num_points - 1) = (cum_local_dists / cum_local_dists(num_points - 2));

    // Initialize Bernstein array
    Mat<T> B = bern(deg, t);

    // Initialize Bezier curve coordinates
    Mat<T> C(num_points, 2);
    C = B * P;

    // Initialize error vector
    Col<T> f = bezerror(deg, data, C, P);

    // Vectorize P
    Col<T> w = arma::vectorise(P);

    // Construct Jacobian
    Mat<T> J = bezjac(deg, data, B, w, f, bounded);

    // Update w
    w -= arma::pinv(J) * f;

    // Calculate new error
    P = arma::reshape(w, deg + 1, 2);

    // Draw the fitting curve in the resolution of the original image
    Col<T> t_final = arma::linspace(0, 1, num_points);
    Mat<T> B_final = bern(deg, t_final);

    return B_final * P;
}

////
//// Explicit instantiation for the template class
////
template class CurveFit<double>;