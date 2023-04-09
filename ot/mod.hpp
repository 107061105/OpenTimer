#ifndef FLOAT_MOD
#define FLOAT_MOD

#include <vector>
#include <complex>
#include <algorithm>
#include <math.h>
#include <iostream>

#define CONVOLVE_MODE 0
#define TIME_STEP 0.001
#define SHRINK_THRESHOLD 1e-6;
// careful with TIME_STEP, double arithmetic can cause imprecision


typedef std::complex<double> Complex;


enum min_max_mode {
    Min,
    Max
};


std::vector<double> convolve(std::vector<double>&, std::vector<double>&);

// inplace FFT
// sign 1 => forward FFT
// sign -1 => inverse FFT
void FFT(std::vector<Complex>&, int);

std::vector<double> fft_convolve(std::vector<double>&, std::vector<double>&);


class float_mod
{
private:
    std::vector<double> dist;
    double st;
    double get_et() {return st + TIME_STEP * (dist.size() - 1);}

    double time_unit = 1;

    friend float_mod float_mod_mm(float_mod&, float_mod&, bool);
public:
	float_mod(const std::vector<double>& Dist, double St = 0.0);
	float_mod(const float_mod& x) {dist = x.dist; st = x.st;}
	float_mod(double);
	float_mod();
	float_mod& operator=(const double&);
	float_mod& operator=(const float_mod&);
    float_mod operator+(const float_mod&) const;
    float_mod operator-(const float_mod&) const;
    float_mod operator+(float) const;
    float_mod operator-(float) const;
    // operator float () const {return st;} // TODO XD
	void dist_norm();
	void shrink();
	void negate();
    void print_info() const;
    const float mean() const;
    const float variance() const;
    const float accum(float) const;
    const float upper() const {return accum(0.9987);}
    const float lower() const {return accum(0.0013);}
    float_mod sacle_time(float);
};


float_mod float_mod_mm(float_mod& x, float_mod&y, bool mode);

std::ostream& operator<<(std::ostream& os, const float_mod& obj);

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
//     std::cout << x.mean() << std::endl;
//     std::cout << x.variance() << std::endl;
//     std::cout << x.accum(0.03) << std::endl;
//     std::cout << x.accum(0.5) << std::endl;
//     std::cout << x.accum(0.997) << std::endl;
//     std::cout << x.accum(1.5) << std::endl;
//     return 0;
// }