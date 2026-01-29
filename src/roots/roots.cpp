#include "roots.hpp"
#include <functional>
#include <cmath>
#include <iostream>

const int max_iter = 100000;
const double tol_f    = 1e-10; 
const double tol_x    = 1e-10;


bool bisection(std::function<double(double)> f,
               double a, double b,
               double *root) {
    int n = 0;
    const double tolerance = 1e-10;
    double c;
    if (f(a) * f(b) >= 0) {
	    return 0;
    }


    while (( n < max_iter)) {


	c = (a+b)/2;

	if (std::abs(f(c)) <= tolerance) {
		*root = c;
		return 1;
	} 

	if (f(a) * f(c) < 0) {
		b = c;
	} else {
		a = c;
	}


	n++;
    }

   *root = (a+b)/2;

    return 0;
}

bool regula_falsi(std::function<double(double)> f,
                  double a, double b,
                  double *root) {






	    if (a > b) std::swap(a, b);


	    double fa = f(a);
	    double fb = f(b);

	    if (fa * fb > 0.0) return false;

	    double c = a;
            // This converages very slowly for this curve
	    for (int n = 0; n < (max_iter * 10); n++) {
		// c = a - f(a)*(b-a)/(f(b)-f(a))

		c = a - fa * (b - a) / (f(b) - f(a));
		double fc = f(c);

		// stopping conditions
		if (std::fabs(fc) <= tol_f) {
		    *root = c;
		    return true;
		}

		// keep the bracket
		if (fa * fc < 0) {
		    b  = c;
		    fb = fc;
		} else {
		    a  = c;
		    fa = fc;
		}
	    }

	    *root = c;      // best estimate after max iterations
	    return false;   // didn't hit tolerance within max_iter

}

bool newton_raphson(std::function<double(double)> f,
                    std::function<double(double)> g,
                    double a, double b, double c,
                    double *root) {


    const double eps_df   = 1e-14;

    double x = c;

    for (int n = 0; n < max_iter; ++n) {
        const double fx = f(x);
        if (std::abs(fx) <= tol_f) {
            *root = x;
            return true;
        }

        const double dfx = g(x);
        if (std::abs(dfx) <= eps_df) {
            return false; // derivative too small
        }

        double x_next = x - fx / dfx;

        if (x_next < a || x_next > b) {
            return false;
        }

        if (std::abs(x_next - x) <= tol_x) {
            *root = x_next;
            return true;
        }

        x = x_next;
    }

    *root = x;      // best estimate after max iterations
    return false;





}

bool secant(std::function<double(double)> f,
            double a, double b, double c,
            double *root) {


    const double eps_denom = 1e-14;  // guard against division by ~0

    double x_prev = a;   // x{n-1}
    double x_curr = c;   // x_n

    double f_prev = f(x_prev);
    double f_curr = f(x_curr);
    for (int n = 0; n < max_iter; ++n) {

        if (std::abs(f_curr) <= tol_f) {
            *root = x_curr;
            return true;
        }

        const double denom = (f_curr - f_prev);
        if (std::abs(denom) <= eps_denom) {
            return false; // secant slope too small
        }

        // Secant update
        double x_next = x_curr - f_curr * (x_curr - x_prev) / denom;

        if (x_next < a || x_next > b) {
		// Fall back to bisection
		x_next = 0.5 * (a+b);
        }

        if (std::abs(x_next - x_curr) <= tol_x) {
            *root = x_next;
            return true;
        }

        x_prev = x_curr;
        f_prev = f_curr;

        x_curr = x_next;
        f_curr = f(x_curr);
    }

    *root = x_curr;   // best estimate

	
    return false;
}

