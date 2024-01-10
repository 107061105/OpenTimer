#ifndef OT_LIBERTY_STATISICAL_UTILITY_HPP_
#define OT_LIBERTY_STATISICAL_UTILITY_HPP_

// #include <ot/headerdef.hpp>
#include <iostream>
#include <fstream>
#include <vector>
#include <complex>

typedef std::complex<float> Complex;

// The data structure of the statisical analysis,
// used on statisical static analysis
namespace Statisical {

#define TIME_STEP 0.1
#define SAMPLE_NUM 100000
#define SHRINK_THRESHOLD 1e-6

enum Distribution_type {
    Constant,
    Gaussian,
	SkewNormal
};

enum Split {
    MIN,
    MAX
};

std::vector<float> generate_SN_samples(int, float, float, float);

std::vector<float> generate_Gaussian_samples(int, float, float);

std::vector<float> calculateProbabilityDensity(const std::vector<float>&, int*);

void saveSamplesToFile(const std::vector<float>&, const std::string&);

void saveProbabilityDensityToFile(const std::vector<float>&, const std::string&, int);

void fft(std::vector<Complex>&, bool);

void fft(std::vector<Complex>&);

void ifft(std::vector<Complex>&);

void fft_convolve(const std::vector<float>&, const std::vector<float>&, std::vector<float>&);

}

#endif
