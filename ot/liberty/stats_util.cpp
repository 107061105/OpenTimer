#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <complex>
#include <algorithm>
#include <iomanip>
#include <random>
#include <ot/static/logger.hpp>
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
 * @brief Generate samples of micmic Skew Normal distribution
 * 
 * @param mean  location parameter
 * @return std::vector<float>, vector of samples
 */
std::vector<float> generate_MicMic_SN_samples(int num_samples, ot::Tran rf, float mean) 
{
    std::vector<float> samples;
    samples.reserve(num_samples);
    float stdev = 0.0f, skew = 0.0f;
    if (VDD == 0.5) 
    {
        if (rf == ot::Tran::FALL) {
            stdev = 0.1817071 * mean;
            skew  = 1.02;
        }
        else {
            stdev = 0.2469286 * mean;
            skew  = 1.42; 
        }
    } 
    else if (VDD == 0.4) 
    {
        if (rf == ot::Tran::FALL) {
            stdev = 0.5110573 * mean;
            skew  = 2.32;
        } else {
            stdev = 0.6690626 * mean;
            skew  = 2.73;
        }
    }
    else 
    {
        std::cerr << "Undefined VDD!!!\n";
    }
    // OT_LOGD("Generate MicMic SN samples, mean/stdev/skew = ", mean, "/", stdev, "/", skew);

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
 * @brief Get the cumulative distribution at point x for a Gaussian distribution
 *
 * @param x
 * @param mean      Mean parameter
 * @param stdev     Standard deviation parameter
 * @return float
 */
float get_Gaussian_cdf(float x, float mean, float stdev) 
{
    boost::math::normal_distribution<float> dist(mean, stdev);
    return boost::math::cdf(dist, x);
}

/**
 * @brief Get the probability density at point x for a Gaussian distribution
 *
 * @param x
 * @param mean      Mean parameter
 * @param stdev     Standard deviation parameter
 * @return float
 */
float get_Gaussian_pdf(float x, float mean, float stdev) 
{
    boost::math::normal_distribution<float> dist(mean, stdev);
    return boost::math::pdf(dist, x);
}

/**
 * @brief Generate pdf of Gaussian distribution
 * 
 * @param mean      Mean parameter 
 * @param start     Record the starting point
 * @return std::vector<float>, vector of samples
 */
std::vector<float> generate_Gaussian_pdf(ot::Tran rf, float mean, int* start) 
{
    float stdev = 0.0f;
    float th = 1e-6;    // threshold
    std::vector<float> pdf;

    if (VDD == 0.5) 
    {
        if (rf == ot::Tran::FALL) {
            stdev = 0.1817071 * mean;
        }
        else {
            stdev = 0.2469286 * mean;
        }
    } 
    else if (VDD == 0.4) 
    {
        if (rf == ot::Tran::FALL) {
            stdev = 0.5110573 * mean;
        } else {
            stdev = 0.6690626 * mean;
        }
    }
    else 
    {
        std::cerr << "Undefined VDD!!!\n";
    }
    // std::cout << "Generate Gaussian pdf, mean/stdev = ";
    // std::cout << mean << "/" << stdev << std::endl;

    // Iterate until CDF is smaller than the threshold
    int init = floor(mean / TIME_STEP);

    while (get_Gaussian_cdf(init * TIME_STEP, mean, stdev) > th) {
        // Subtract the time step from x
        init--;
    }
    // record the starting point
    *start = init;
    // std::cout << "PDF Starting point: " << init << ", cdf: " << get_Gaussian_cdf(init * TIME_STEP, mean, stdev) << std::endl;

    // get the probability density function
    float x = init * TIME_STEP;

    while (get_Gaussian_cdf(x, mean, stdev) <= (1.0f - th)) 
    {
        float val = get_Gaussian_pdf(x, mean, stdev);
        pdf.push_back(val);
        x += TIME_STEP;
    }

    return pdf;
}

/**
 * @brief Get the cumulative distribution by at point x
 *
 * @param x
 * @param location  location parameter
 * @param scale     scale parameter
 * @param shape     shape parameter
 * @return float
 */
float get_SN_cdf(float x, float location, float scale, float shape) 
{
    boost::math::skew_normal_distribution<float> dist(location, scale, shape);
    return boost::math::cdf(dist, x);
}

/**
 * @brief Get the probability density by at point x
 *
 * @param x
 * @param location  location parameter
 * @param scale     scale parameter
 * @param shape     shape parameter
 * @return float
 */
float get_SN_pdf(float x, float location, float scale, float shape) 
{
    boost::math::skew_normal_distribution<float> dist(location, scale, shape);
    return boost::math::pdf(dist, x);
}

/**
 * @brief Generate pdf of micmic Skew Normal distribution
 * 
 * @param mean  location parameter 
 * @param start record the starting point
 * @return std::vector<float>, vector of samples
 */
std::vector<float> generate_MicMic_SN_pdf(ot::Tran rf, float mean, int* start) 
{
    float stdev = 0.0f, skew = 0.0f;
    float th = 1e-6;    // threshold
    std::vector<float> pdf;

    if (VDD == 0.5) 
    {
        if (rf == ot::Tran::FALL) {
            stdev = 0.1817071 * mean;
            skew  = 1.02;
        }
        else {
            stdev = 0.2469286 * mean;
            skew  = 1.42; 
        }
    } 
    else if (VDD == 0.4) 
    {
        if (rf == ot::Tran::FALL) {
            stdev = 0.5110573 * mean;
            skew  = 2.32;
        } else {
            stdev = 0.6690626 * mean;
            skew  = 2.73;
        }
    }
    else 
    {
        std::cerr << "Undefined VDD!!!\n";
    }
    // std::cout << "Generate MicMic SN pdf, mean/stdev/skew = ";
    // std::cout << mean << "/" << stdev << "/" << skew << std::endl;

    // Convert mean, stdev, and skew to location, scale, and shape parameters
    float location = mean;
    float scale = stdev * std::sqrt(1 - (2.0 / M_PI) * (skew / std::sqrt(1 + skew * skew)));
    float shape = skew;

    // Iterate until CDF is smaller than the threshold
    int init = floor(mean / TIME_STEP);

    while (get_SN_cdf(init * TIME_STEP, location, scale, shape) > th) {
        // Subtract the time step from x
        init--;
    }
    // record the starting point
    *start = init;
    // std::cout << "PDF Starting point: " << init << ", cdf: " << get_SN_cdf(init * TIME_STEP, location, scale, shape) << std::endl;

    // get the probability density function
    float x = init * TIME_STEP;

    while (get_SN_cdf(x, location, scale, shape) <= (1.0f - th)) 
    {
        float val = get_SN_pdf(x, location, scale, shape);
        pdf.push_back(val);
        x += TIME_STEP;
    }

    return pdf;
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


void convolution(const std::vector<float>& dist1, const std::vector<float>& dist2, std::vector<float>& result) {
    int len1 = dist1.size();
    int len2 = dist2.size();
    int resultSize = len1 + len2 - 1;

    result.resize(resultSize, 0.0);

    for (int i = 0; i < len1; ++i) {
        for (int j = 0; j < len2; ++j) {
            result[i + j] += dist1[i] * dist2[j];
        }
    }

    for (auto &x: result) {
        x *= TIME_STEP;
    }
}


}