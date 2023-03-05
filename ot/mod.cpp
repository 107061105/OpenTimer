#ifndef FLOAT_MOD
#define FLOAT_MOD

#include <vector>
#include <valarray>
#include <complex>
#include <algorithm>
#include <math.h>
#include <iostream>

#define CONVOLVE_MODE 0
#define TIME_STEP 0.01
// careful with TIME_STEP, double arithmetic can cause imprecision


// using namespace std;
// const double PI = 3.141592653589793238460;
typedef std::complex<double> Complex;


enum min_max_mode {
    Min,
    Max
};


std::vector<double> convolve(std::vector<double>& a, std::vector<double>& b) {
    int n = a.size(), m = b.size();
    std::vector<double> res(n+m-1);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            res[i+j] += a[i] * b[j];
        }
    }
    return res;
}

// inplace FFT
// sign 1 => forward FFT
// sign -1 => inverse FFT
void FFT(std::vector<Complex>& x, int sign) {
    int n = x.size();
    if (n <= 1) return;
    std::vector<Complex> even(n/2), odd(n/2);
    for (int i = 0; i < n/2; i++) {
        even[i] = x[2*i];
        odd[i] = x[2*i+1];
    }
    FFT(even, sign);
    FFT(odd, sign);
    Complex w(1, 0);
    Complex wn(cos(2*M_PI/n), sign*sin(2*M_PI/n));
    for (int i = 0; i < n/2; i++) {
        x[i] = even[i] + w*odd[i];
        x[i+n/2] = even[i] - w*odd[i];
        w *= wn;
    }
}

std::vector<double> fft_convolve(std::vector<double>& a, std::vector<double>& b) {
    int nn = a.size() + b.size() - 1;
    int n = 1;
    while (n < nn) n <<= 1; // find FFT min length
    std::vector<Complex> x(n), y(n);
    for (int i = 0; i < a.size(); i++) x[i] = Complex(a[i], 0);
    for (int i = 0; i < b.size(); i++) y[i] = Complex(b[i], 0);
    FFT(x, 1);
    FFT(y, 1);
    for (int i = 0; i < n; i++) x[i] *= y[i]; // multiplication
    FFT(x, -1); // inverse FFT
    std::vector<double> res(nn);
    for (int i = 0; i < nn; i++) res[i] = x[i].real() / n; // take the real part and normalize
    for (int i = nn; i < n; i++) if (x[i].real() != 0) std::cout << "error! fft convolution has non-zero value outside res_dist (not in n+m-1 values)";
    return res;
}


class float_mod
{
private:
    std::vector<double> dist;
    double st;
    double get_et() {return st + TIME_STEP * (dist.size() - 1);}

    friend float_mod float_mod_mm(float_mod&, float_mod&, bool);
    // friend class ot::Pin;
public:
	float_mod(const std::vector<double>& Dist = std::vector<double>(1, 1), double St = 0.0);
	float_mod(const float_mod& x) {dist = x.dist; st = x.st;}
	float_mod(double St = 0.0);
	float_mod() = default;
	float_mod& operator=(const double&);
	float_mod& operator=(const float_mod&);
    float_mod operator+(float_mod&);
    float_mod operator-(float_mod&);
    float_mod operator+(float);
    float_mod operator-(float);
    operator float () const {return st;} // TODO XD
	void dist_norm();
	void shrink();
	void negate();
    void print_info();
};

float_mod::float_mod(const std::vector<double>& Dist, double St) {
	dist = Dist;
	st = St;
}

float_mod::float_mod(double St) {
    std::vector<double> temp(1, 1);
	dist = temp;
	st = St;
}

float_mod::float_mod() {
    std::vector<double> temp(1, 1);
	dist = temp;
	st = 0.0;
}

float_mod& float_mod::operator=(const double& St) {
    std::vector<double> temp(1, 1);
	dist = temp;
	st = St;
    return *this;
}

float_mod& float_mod::operator=(const float_mod& x) {
	dist = x.dist;
	st = x.st;
    return *this;
}

float_mod float_mod::operator+(float_mod& x) {
	double res_st = st + x.st;
	std::cout << "convolution st: " << res_st << std::endl;

	std::vector<double> res;
	if (CONVOLVE_MODE == 0) {
		res = convolve(dist, x.dist);
	} else if (CONVOLVE_MODE == 1) {
		res = fft_convolve(dist, x.dist);
	}
	
	std::cout << "convolution result: ";
	for (int i = 0; i < res.size(); i++) {
		std::cout << res[i] << " ";
	}
	std::cout << std::endl;

	float_mod temp(res, res_st);
	return temp;
}

float_mod float_mod::operator-(float_mod& x) {
	x.negate();
    float_mod res = *this + x;
    x.negate();
    return res;
}

float_mod float_mod::operator+(float x) {
	st += x;
    return *this;
}

float_mod float_mod::operator-(float x) {
	st -= x;
    return *this;
}

void float_mod::dist_norm() {
    double sum = 0;
    for (auto& n : dist) sum += n;
	for (auto& n : dist) n /= sum;
}

void float_mod::negate() {
	// set start time to negative end time
	st = -get_et();
    // reverse the distribution
	std::reverse(dist.begin(), dist.end());
}

void float_mod::shrink() {
	// get new start time and end time, skipping zeros (or low values)
    double threshold = 0;
    int new_st_index = 0;
    while (dist[new_st_index] <= threshold) {
        new_st_index++;
        st += TIME_STEP;
    }
    dist.erase(dist.begin(), dist.begin() + new_st_index);
    while (dist.back() <= threshold) {
        dist.pop_back();
    }
}

void float_mod::print_info() {
	std::cout << "st: " << st << std::endl;
	std::cout << "dist: ";
	for (int i = 0; i < dist.size(); i++) {
		std::cout << dist[i] << " ";
	}
	std::cout << std::endl;
}

float_mod float_mod_mm(float_mod& x, float_mod&y, bool mode) {
    double st = std::min(x.st, y.st);
    double et = std::max(x.get_et(), y.get_et());
    int length = std::round((et - st) / TIME_STEP) + 1;
    std::cout << "min/max info (st, et, length): " << st << ", " << et << ", " << length << std::endl;

    int x_offset = std::round((x.st-st)/TIME_STEP);
    int y_offset = std::round((y.st-st)/TIME_STEP);
    std::vector<double> res_dist(length, 0);
    std::vector<int> index_order;
    double px_sum = 0; // x's probability sum
    double py_sum = 0;
    if (mode == Min) {
        for (int i = 0; i < length; i++) index_order.push_back(i);
    } else if (mode == Max) {
        for (int i = length - 1; i >= 0; i--) index_order.push_back(i);
    }
    for (const auto& i : index_order) {
        double xp = 0; // x's probability at this time
        double yp = 0;
        if (i >= x_offset && i - x_offset < x.dist.size()) xp = x.dist[i - x_offset];
        if (i >= y_offset && i - y_offset < x.dist.size()) yp = y.dist[i - y_offset];
        res_dist[i] = py_sum * xp + px_sum * yp + xp * yp;
        px_sum = px_sum + xp;
        py_sum = py_sum + yp;
    }

    float_mod res(res_dist, st);
    res.shrink();
    return res;
}

#endif

// int main(int argc, char const *argv[])
// {
//     double a_st = 5;
//     double b_st = 5;
// 	std::vector<double> a = {1, 2, 3, 4};
//     std::vector<double> b = {1, 1, 10};
//     float_mod x(a, a_st);
//     float_mod y(b, b_st);
//     float_mod z = x + y;
//     z = x - y;
//     z.dist_norm();
//     z.print_info();
//     x.dist_norm();
//     y.dist_norm();
//     z = float_mod_mm(x, y, Min);
//     z.print_info();
//     z = float_mod_mm(x, y, Max);
//     z.print_info();
//     return 0;
// }