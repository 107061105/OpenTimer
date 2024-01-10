#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <complex>
#include <algorithm>
#include <iomanip>
#include <random>
#include <boost/math/distributions/skew_normal.hpp>
#include "stats_util.hpp"

const float PI = 3.141592653589793238462643383279502884;

namespace Statisical {

/**
 * @brief Generate samples of Skew Normal distribution
 * 
 * @param mean  location parameter 
 * @param stdev scale parameter, which is calculated by
 *              stdev * sqrt(1 - (2.0 / PI) * (skew / sqrt(1 + skew^2)))
 * @param skew  skewness
 * @return std::vector<float>, vector of samples
 */
std::vector<float> generate_SN_samples(int num_samples, float mean, float stdev, float skew) 
{
    std::vector<float> samples;
    samples.reserve(num_samples);

    // Convert mean, stdev, and skew to location, scale, and shape parameters
    float location = mean;
    float scale = stdev * std::sqrt(1 - (2.0 / M_PI) * (skew / std::sqrt(1 + skew * skew)));
    float shape = skew;

    // Setup generators
    std::random_device rd;
    std::default_random_engine noise_generator;

    // Sample from a uniform distribution i.e., [0,1)
    std::uniform_real_distribution<float> uniform_dist(0, 1.0);

    auto skew_norm_dist = boost::math::skew_normal_distribution<float>(location, scale, shape);

    // Use the probability from the uniform distribution with the percent point
    // function of the skew_normal
    int i = 0;
    while (i < num_samples) {
        noise_generator.seed(rd());
        auto probability = uniform_dist(noise_generator);
        float skew_normal_sample_point = boost::math::quantile(skew_norm_dist, probability);
        if (skew_normal_sample_point != INFINITY && skew_normal_sample_point != -INFINITY) {
            samples.push_back(skew_normal_sample_point);
            i++;
        }
    }

    return samples;
}

/**
 * @brief Generate samples of Gaussian distribution
 * 
 * @param mean 
 * @param stdev standard deviation
 * @return std::vector<float>, vector of samples
 */
std::vector<float> generate_Gaussian_samples(int num_samples, float mean, float stdev) 
{
    std::vector<float> samples;
    samples.reserve(num_samples);

    // Setup generators
    std::random_device rd;
    std::default_random_engine normal_generator(rd());

    // Sample from a normal distribution
    std::normal_distribution<float> normal_dist(mean, stdev);

    int i = 0;
    while (i < num_samples) {
        float normal_sample_point = normal_dist(normal_generator);
        if (normal_sample_point != INFINITY && normal_sample_point != -INFINITY) {
            samples.push_back(normal_sample_point);
            i++;
        }
    }

    return samples;
}

/**
 * @brief Calculate and return the probability density
 *
 * @param data vector that recording samples
 * @param start record the starting point
 * @return std::vector<float>, the probability of data
 */
std::vector<float> calculateProbabilityDensity(const std::vector<float>& data, int* start) 
{
    int min = floor(*min_element(data.begin(), data.end()) / TIME_STEP);
    int max = ceil(*max_element(data.begin(), data.end()) / TIME_STEP);
    int numBins = max - min;
    
    *start = min;

    // std::cout << "max value, min value: " << *max_element(data.begin(), data.end()) << ", " << *min_element(data.begin(), data.end()) << std::endl;
    // std::cout << "max, min, num: " <<  max << ", " << min << ", " << numBins << std::endl;

    std::vector<int> histogram(numBins, 0);

    for (float value : data) {
        int binIndex = static_cast<int>((value / TIME_STEP) - min);
        if (binIndex >= numBins) {
            std::cerr << "value, max, min, idx: " << value << ", " << max << ", " << min << ", " << binIndex << std::endl;
        }
        histogram[binIndex]++;
    }

    std::vector<float> probabilityDensity(numBins, 0.0);

    for (int i = 0; i < numBins; i++) {
        // Calculate probability density and store in the vector
        float probability = static_cast<float>(histogram[i]) / (data.size() * TIME_STEP);
        probabilityDensity[i] = probability;
    }

    return probabilityDensity;
}

/**
 * @brief Save samples to a text file
 * 
 * @param samples vector of samples
 * @param filename output file name
 */
void saveSamplesToFile(const std::vector<float>& samples, const std::string& filename) 
{
    std::ofstream file(filename.c_str()); // Use c_str() to convert to const char*

    if (file.is_open()) {
        for (const auto &sample : samples) {
            file << sample << "\n";
        }
        file.close();
        // std::cout << "Saved samples to " << filename << std::endl;
    } else {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }
}

/**
 * @brief Save a pdf to file
 * 
 * @param dist vector that recording a pdf
 * @param filename output file name
 * @param st starting time of the pdf
 */
void saveProbabilityDensityToFile(const std::vector<float>& dist, const std::string& filename, int st) 
{
    std::ofstream outputFile(filename);

    if (!outputFile.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    // cout << *max_element(data.begin(), data.end()) << endl;
    for (int i = 0; i < static_cast<int>(dist.size()); i++) {
        float binStart = (st + i) * TIME_STEP;
        float binEnd = binStart + TIME_STEP;
        outputFile << std::fixed << std::setprecision(6) << "[" << binStart << ", " << binEnd << "): " << dist[i] << std::endl;
    }

    outputFile.close();
}

/**
 * @brief Recursive FFT function
 * 
 * @param dist vector to be FFT
 * @param invert if true, perfrom inverse FFT
 */
void fft(std::vector<Complex>& dist, bool invert) 
{
    int n = dist.size();
    if (n <= 1) {
        return;
    }

    std::vector<Complex> a0(n / 2), a1(n / 2);
    for (int i = 0, j = 0; i < n; i += 2, ++j) {
        a0[j] = dist[i];
        a1[j] = dist[i + 1];
    }

    fft(a0, invert);
    fft(a1, invert);

    float angle = 2 * PI / n * (invert ? -1 : 1);
    Complex w(1), wn(cos(angle), sin(angle));

    for (int i = 0; i < n / 2; ++i) {
        Complex t = w * a1[i];
        dist[i] = a0[i] + t;
        dist[i + n / 2] = a0[i] - t;
        w *= wn;
    }
}

/**
 * @brief FFT wrapper function
 * 
 * @param dist vector to be FFT
 */
void fft(std::vector<Complex>& dist) 
{
    fft(dist, false);
}

/**
 * @brief Inverse FFT wrapper function
 * 
 * @param dist vector to be iFFT
 */
void ifft(std::vector<Complex>& dist) 
{
    fft(dist, true);
    for (size_t i = 0; i < dist.size(); ++i) {
        dist[i] /= dist.size();
    }
}

/**
 * @brief Convolve two distribution using FFT convolution
 * 
 * @param dist1
 * @param dist2
 * @param result the convolved result
 */
void fft_convolve(const std::vector<float>& dist1, const std::vector<float>& dist2, std::vector<float>& result) 
{
    size_t n = 1;
    while (n < std::max(dist1.size(), dist2.size())) {
        n *= 2;
    }
    n *= 2;

    std::vector<Complex> dist1_complex(n, 0.0), dist2_complex(n, 0.0);

    // Copy the input sequences into complex vectors for FFT
    for (size_t i = 0; i < dist1.size(); ++i) {
        dist1_complex[i] = dist1[i];
    }

    for (size_t i = 0; i < dist2.size(); ++i) {
        dist2_complex[i] = dist2[i];
    }

    // Perform FFT
    fft(dist1_complex);
    fft(dist2_complex);

    // Point-wise multiplication
    for (size_t i = 0; i < n; ++i) {
        dist1_complex[i] *= dist2_complex[i];
    }

    // Inverse FFT
    ifft(dist1_complex);

    // Copy the real part of the result back to the float vector
    size_t result_size = dist1.size() + dist2.size() - 1;
    result.resize(result_size);
    for (size_t i = 0; i < result_size; ++i) {
        result[i] = dist1_complex[i].real() * TIME_STEP;
    }
}

}