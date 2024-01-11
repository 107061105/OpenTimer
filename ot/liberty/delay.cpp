#include <iostream>
#include <cmath>
#include <vector>
#include <numeric>
#include <algorithm>
#include <boost/math/distributions/skew_normal.hpp>
#include "delay.hpp"

namespace Statisical {

/**
 * @brief Initialize with a value, distribution has no variation
 * 
 * @param value constant value
 */
Distribution::Distribution(float value = 0.0f)
{
    _type  = Distribution_type::Constant;
    _pdf   = std::vector<float>(1, 1.0f);
    _start = value / TIME_STEP;
}

/**
 * @brief Initialize with a pdf and its start time
 * 
 * @param pdf probability density
 * @param st  starting time
 */
Distribution::Distribution(const std::vector<float> &pdf, int st = 0)
{
    _type  = Distribution_type::SkewNormal;
    _pdf   = pdf;
    _start = st;
}

/**
 * This distribution type is to model the skewed distriubtions of academic 
 * libraries, since there is only nldm data we can get, therefore, we use 
 * the coeffcient of variation (CV) to get the approximated stdev and skew.
 * 
 * @brief Initialize with given micmic SN distribution mean
 * 
 * @param mean mean
 * @param sta standard deviation
 * @param skew skewness, which is 0 if Gaussian distribution
 */
Distribution::Distribution(Distribution_type type, ot::Tran rf, float mean)
{
    std::vector<float> samples;
    float std = 0.0f, skew = 0.0f;

    assert(type == Distribution_type::MicmicSN);
    _type = type;
    samples = generate_MicMic_SN_samples(SAMPLE_NUM, rf, mean);
    _pdf = calculateProbabilityDensity(samples, &_start);
}


/**
 * @brief Initialize with given Gaussian/SN distribution moments
 * 
 * @param mean mean
 * @param sta standard deviation
 * @param skew skewness, which is 0 if Gaussian distribution
 */
Distribution::Distribution(Distribution_type type, float mean, float std, float skew = 0.0)
{
    std::vector<float> samples;

    if (type == Distribution_type::Gaussian) {
        assert(skew == 0.0); // skewness is zero if it is Gaussian distribution
        _type = type;
        samples = generate_Gaussian_samples(SAMPLE_NUM, mean, std);
    } 
    else if (type == Distribution_type::SkewNormal) {
        _type = type;
        samples = generate_SN_samples(SAMPLE_NUM, mean, std, skew);
    } else {
        std::cerr << "Invalid distribution type!!!\n";
    }

    _pdf = calculateProbabilityDensity(samples, &_start);
}

/**
 * @brief Overload operator + for Distribution, 
 *
 * @param value the constant value
 * @return Distribution
 */
Distribution Distribution::operator+(float value) 
{
    return Distribution(_pdf, _start + static_cast<int>(value / TIME_STEP));
}

/**
 * @brief Overload operator + for Distribution
 *
 * @param rhs the distribution
 * @return Distribution
 */
Distribution Distribution::operator+(Distribution &rhs) 
{
    if (rhs.get_type() == Distribution_type::Constant) 
    {
        assert(rhs.get_bin_num() == 1);
        return Distribution(_pdf, _start + rhs.get_start_point());
    } 
    else 
    {
        std::vector<float> result;
        int start_time = 0;

        start_time = _start + rhs.get_start_point();
        fft_convolve(_pdf, rhs.get_pdf(), result);
        Distribution temp(result, start_time);
        temp.shrink();

        return temp;
    }
}

/**
 * @brief Overload operator - for Distribution, 
 *
 * @param value the constant value
 * @return Distribution
 */
Distribution Distribution::operator-(float value) 
{
    return Distribution(_pdf, _start - static_cast<int>(value / TIME_STEP));
}

/**
 * @brief Overload operator - for Distribution
 *
 * @param rhs the distribution
 * @return Distribution
 */
Distribution Distribution::operator-(Distribution &rhs) 
{
    Distribution temp = rhs;
    temp.negate();
    return *this + temp;
}

/**
 * @brief Normalize the pdf
 */
void Distribution::norm()
{
    float sum = 0.0;

    for (float &value : _pdf) {
        value = std::max(value, 0.f);
        sum += value * TIME_STEP;
    }

    for (float &value : _pdf) {
        value /= sum;
    }
}

/**
 * @brief Negate the pdf
 */
void Distribution::negate()
{
    // set start time to negative end time
    _start = -1 * get_end_point();
    // reverse the distribution
    std::reverse(_pdf.begin(), _pdf.end());
}

/**
 * @brief Shrink the pdf, ignore the tail under threshold
 */
void Distribution::shrink()
{
    // get new start time and end time, skipping zeros (or low values)
    while (!_pdf.empty() && _pdf.front() < SHRINK_THRESHOLD)
    {
        _pdf.erase(_pdf.begin());
        _start++;
    }
    while (!_pdf.empty() && _pdf.back() < SHRINK_THRESHOLD)
    {
        _pdf.pop_back();
    }
    norm();
}

/**
 * @brief Get pdf[i]
 *
 * @param i 
 * @return float
 */
float Distribution::get_ith_pdf(int i) const
{
    if (i < 0)
        std::cerr << "Function \"get_ith_pdf\": index i is < 0 !!\n";
    else if (i >= get_bin_num())
        std::cerr << "Function \"get_ith_pdf\": index i is out of range !!\n";
    return _pdf[i];
}

/**
 * @brief Get the ±3-sigma of the pdf
 *
 * @param type Max is 3 sigma (99.865%) while Min is -3 sigma (0.135%)
 * @return float
 */
float Distribution::get_3_sigma(ot::Split type)
{
    float total = 0.0f;
    float sum = 0.0f;

    total = std::accumulate(_pdf.begin(), _pdf.end(),
                                decltype(_pdf)::value_type(0));

    if (type == ot::Split::MAX) {
        for (int i = 0; i < get_bin_num(); i++)
        {
            sum += _pdf[i];
            if (sum / total >= 0.99865) return TIME_STEP * (i + _start);
        }
    } else if (type == ot::Split::MIN) {
        for (int i = get_bin_num() - 1; i >= 0; i--)
        {
            sum += _pdf[i];
            if (sum / total >= 0.99865) return TIME_STEP * (i + _start);
        }
    }
    std::cerr << "error, in Statisical_delay::accum(p), overflow" << std::endl;
    return 0;
}

/**
 * @brief Print the status of the pdf
 */
void Distribution::print_status()
{
    std::cout << "***************************************************\n";
    std::cout << "Start time: " << get_start_time() << ", End time: " << get_end_time() << std:: endl;
    std::cout << "number of bins: " << get_bin_num() << std::endl;
    std::cout << "3 sigma: " << get_3_sigma(ot::Split::MAX) << std::endl;
    std::cout << "***************************************************\n";
}

/**
 * @brief Perform max operation of two distribution
 *
 * @param dist1
 * @param dist2
 * @return Distribution of max
 */
Distribution max(const Distribution &dist1, const Distribution &dist2) 
{
    int st = std::min(dist1.get_start_point(), dist2.get_start_point());
    int et = std::max(dist1.get_end_point(), dist2.get_end_point());
    int offset1 = dist1.get_start_point() - st;
    int offset2 = dist2.get_start_point() - st;
    int length  = et - st + 1;

    std::vector<float> res_dist(length, 0);

    float p1_sum = 0;
    float p2_sum = 0;
    
    for (int i = 0; i < length; i++)
    {
        float p1 = 0;
        float p2 = 0;

        if (i >= offset1 && i - offset1 < dist1.get_bin_num()) {
            p1 = dist1.get_ith_pdf(i - offset1);
        }
        if (i >= offset2 && i - offset2 < dist2.get_bin_num()) {
            p2 = dist2.get_ith_pdf(i - offset2);
        }

        res_dist[i] = (p2_sum * p1 + p1_sum * p2) * TIME_STEP;
        p1_sum += p1;
        p2_sum += p2;
    }

    Distribution temp(res_dist, st);

    temp.shrink();

    return temp;
}

/**
 * @brief Perform min operation of two distribution
 *
 * @param dist1
 * @param dist2
 * @return Distribution of min
 */
Distribution min(const Distribution &dist1, const Distribution &dist2) 
{
    int st = std::min(dist1.get_start_point(), dist2.get_start_point());
    int et = std::max(dist1.get_end_point(), dist2.get_end_point());
    int offset1 = dist1.get_start_point() - st;
    int offset2 = dist2.get_start_point() - st;
    int length  = et - st + 1;

    std::vector<float> res_dist(length, 0);

    float p1_sum = 0;
    float p2_sum = 0;
    
    for (int i = length - 1; i >= 0; i--)
    {
        float p1 = 0;
        float p2 = 0;

        if (i >= offset1 && i - offset1 < dist1.get_bin_num()) {
            p1 = dist1.get_ith_pdf(i - offset1);
        }
        if (i >= offset2 && i - offset2 < dist2.get_bin_num()) {
            p2 = dist2.get_ith_pdf(i - offset2);
        }

        res_dist[i] = (p2_sum * p1 + p1_sum * p2) * TIME_STEP;
        p1_sum += p1;
        p2_sum += p2;
    }

    Distribution temp(res_dist, st);

    temp.shrink();

    return temp;
}

}
