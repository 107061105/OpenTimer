#include <iostream>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <ot/static/logger.hpp>
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
    _value = value;
}

/**
 * @brief Initialize with a pdf and its start time
 * 
 * @param pdf probability density
 * @param st  starting time
 */
Distribution::Distribution(std::vector<float> &pdf, float value, int st)
{
    _type  = Distribution_type::SkewNormal;
    _value = value;
    _pdf   = std::move(pdf);
    _start = st;
}

/**
 * This distribution type is to model the skewed distriubtions of academic 
 * libraries, since there is only nldm data we can get, therefore, we use 
 * the coeffcient of variation (CV) to get the approximated stdev and skew.
 * 
 * @brief Initialize with given micmic SN distribution mean
 * 
 * @param type need to be Distribution_type::MicmicSN
 * @param rf   rise/fall
 * @param mean mean
 */
Distribution::Distribution(Distribution_type type, ot::Tran rf, float mean)
{
    std::vector<float> samples, pdf;
    int start_time = 0;

    assert(type == Distribution_type::MicmicSN);
    pdf     = generate_MicMic_SN_pdf(rf, mean, &start_time);

    _type  = type;
    _value = mean;
    _pdf   = std::move(pdf);
    _start = start_time;
}

/**
 * @brief Initialize with given Gaussian/SN distribution moments
 * 
 * @param mean mean
 * @param sta standard deviation
 * @param skew skewness, which is 0 if Gaussian distribution
 */
// Distribution::Distribution(Distribution_type type, float mean, float std, float skew = 0.0)
// {
//     std::vector<float> samples, pdf;
//     int start_time = 0;
// 
//     if (type == Distribution_type::Gaussian) {
//         assert(skew == 0.0); // skewness is zero if it is Gaussian distribution
//         samples = generate_Gaussian_samples(SAMPLE_NUM, mean, std);
//     } 
//     else if (type == Distribution_type::SkewNormal) {
//         samples = generate_SN_samples(SAMPLE_NUM, mean, std, skew);
//     } else {
//         std::cerr << "Invalid distribution type!!!\n";
//     }
//     pdf = calculateProbabilityDensity(samples, &start_time);
//     samples.clear();
// 
//     _type  = type;
//     _pdf   = std::move(pdf);
//     _start = start_time;
// }

/**
 * @brief Overload operator + for Distribution, 
 *
 * @param value the constant value
 * @return Distribution
 */
// Distribution Distribution::operator+(float value) 
// {
//     OT_LOGE("Use constant operator +");
//     int offset = static_cast<int>(value / TIME_STEP);
//     std::vector<float> copy = _pdf;
//     return Distribution(_pdf, get_value() + value, get_start_point() + offset);
// }

/**
 * @brief Overload operator + for Distribution
 *
 * @param rhs the distribution
 * @return Distribution
 */
Distribution Distribution::operator+(const Distribution &rhs) 
{
    if (is_constant()) {
        if (rhs.is_constant()) {
            return Distribution(get_constant() + rhs.get_constant());
        }
        else {
            Distribution temp = rhs;
            return temp + *this;
        }
    }
    else if (rhs.is_constant()) {
        assert(!is_constant());
        int offset = static_cast<int>(rhs.get_constant() / TIME_STEP);
        std::vector<float> copy = _pdf;
        return Distribution(copy, get_constant() + rhs.get_constant(), get_start_point() + offset);
    } 
    else {
        std::vector<float> result;
        int start_time = 0;

        // get the sum of two distribution by convolution
        start_time = get_start_point() + rhs.get_start_point();
        fft_convolve(_pdf, rhs.get_pdf(), result);

        // create resulting distribution
        Distribution temp(result, get_constant() + rhs.get_constant(), start_time);
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
// Distribution Distribution::operator-(float value) 
// {
//     OT_LOGE("Use constant operator -");
//     int offset = static_cast<int>(value / TIME_STEP);
//     std::vector<float> copy = _pdf;
//     return Distribution(_pdf, get_value() - value, get_start_point() - offset);
// }

/**
 * @brief Overload operator - for Distribution
 *
 * @param rhs the distribution
 * @return Distribution
 */
Distribution Distribution::operator-(const Distribution &rhs) 
{
    Distribution temp = rhs;
    temp.negate();
    return *this + temp;
}

/**
 * @brief Overload operator <, only for constant Distribution
 *
 * @param rhs the constant distribution
 * @return bool
 */
// bool Distribution::operator<(const Distribution &rhs) 
// {
//     if (get_type() != rhs.get_type()) {
//         std::cerr << "Invalid comparison for operator < !!\n";
//         std::cerr << "lhs is " << to_string(get_type()) << ", and ";
//         std::cerr << "rhs is " << to_string(rhs.get_type()) << std::endl;
//     }
//     OT_LOGD("constant distribution comparator <");
//     return get_value() < rhs.get_value();
// }

/**
 * @brief Overload operator >, only for constant Distribution
 *
 * @param rhs the constant distribution
 * @return bool
 */
// bool Distribution::operator>(const Distribution &rhs) 
// {
//     if (get_type() != rhs.get_type()) {
//         std::cerr << "Invalid comparison for operator > !!\n";
//         std::cerr << "lhs is " << to_string(get_type()) << ", and ";
//         std::cerr << "rhs is " << to_string(rhs.get_type()) << std::endl;
//     }
//     OT_LOGD("constant distribution comparator >");
//     return get_value() > rhs.get_value();
// }

/**
 * @brief Normalize the pdf
 */
void Distribution::norm()
{
    if (!is_constant()) 
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
}

/**
 * @brief Negate the pdf
 */
void Distribution::negate()
{
    _value = -1 * _value;
    // set start time to negative end time
    if (_start) {
        _start = -1 * get_end_point();
    }
    // reverse the distribution
    std::reverse(_pdf.begin(), _pdf.end());
}

/**
 * @brief Shrink the pdf, ignore the tail under threshold
 */
void Distribution::shrink()
{
    if (!is_constant()) {
        // get new start time and end time, skipping zeros (or low values)
        while (_pdf.front() < SHRINK_THRESHOLD)
        {
            _pdf.erase(_pdf.begin());
            (*_start)++;
        }
        while (_pdf.back() < SHRINK_THRESHOLD)
        {
            _pdf.pop_back();
        }
        norm();
    }
}

/**
 * @brief Get pdf[i]
 *
 * @param i 
 * @return float
 */
float Distribution::get_ith_pdf(int i) const
{
    assert(!is_constant()); // get_ith_pdf
    if (i < 0)
        std::cerr << "Function \"get_ith_pdf\": index i is < 0 !!\n";
    else if (i >= get_bin_num())
        std::cerr << "Function \"get_ith_pdf\": index i is out of range !!\n";
    return _pdf[i];
}

/**
 * @brief Get the ±3 sigma of the pdf
 *
 * @param type Max is 3 sigma (99.865%) while Min is -3 sigma (0.135%)
 * @return float
 */
float Distribution::get_3_sigma(ot::Split el) const
{
    float total = 0.0f;
    float sum = 0.0f;

    assert(!is_constant()); // get_3_sigma
    total = std::accumulate(_pdf.begin(), _pdf.end(),
                            decltype(_pdf)::value_type(0));

    if (el == ot::Split::MAX) {
        for (int i = 0; i < get_bin_num(); i++)
        {
            sum += _pdf[i];
            if (sum / total >= 0.99865) return TIME_STEP * (i + *_start);
        }
    } else if (el == ot::Split::MIN) {
        for (int i = get_bin_num() - 1; i >= 0; i--)
        {
            sum += _pdf[i];
            if (sum / total >= 0.99865) return TIME_STEP * (i + *_start);
        }
    }
    std::cout << sum << " " << total << _pdf.size() << std::endl;
    std::cerr << "error, in Distribution::get_3_sigma, overflow" << std::endl;
    return 0;
}

/**
 * @brief Get the estimasted value of the pdf, return _value if the 
 *        distribution is constant, ±3 sigma if others
 *
 * @param type Max is 3 sigma (99.865%) while Min is -3 sigma (0.135%)
 * @return float
 */
float Distribution::get_value(ot::Split el) const 
{
    if (is_constant()) {
        return _value;
    } else {
        return get_3_sigma(el);
    }
}

/**
 * @brief Print the status of the pdf
 */
void Distribution::print_status()
{
    OT_LOGD("***************************************************");
    OT_LOGD("Distribution type: ", to_string(get_type()));
    OT_LOGD("Constant value: ", get_constant());
    if (!is_constant()) {
        OT_LOGD("Start time: ", get_start_time(), ", End time: ", get_end_time());
        OT_LOGD("number of bins: ", get_bin_num());
        OT_LOGD("3 sigma: ", get_3_sigma(ot::Split::MAX));
    }
    OT_LOGD("***************************************************");
}

/**to_string(get_type())
 * @brief Perform max operation of two distribution
 *
 * @param dist1
 * @param dist2
 * @return Distribution of max
 */
Distribution max(const Distribution &dist1, const Distribution &dist2) 
{
    assert(dist1.get_type() == dist2.get_type());   // max operation need same type
    float value = (dist1.get_constant() > dist2.get_constant()) ? dist1.get_constant() : dist2.get_constant();
    

    if (dist1.is_constant()) {
        // OT_LOGD("Max operation of two constant dists");
        return Distribution(value);
    } 
    else {
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

        Distribution temp(res_dist, value, st);

        temp.shrink();

        return temp;
    }
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
    assert(dist1.get_type() == dist2.get_type());   // max operation need same type
    float value = (dist1.get_constant() < dist2.get_constant()) ? dist1.get_constant() : dist2.get_constant();

    if (dist1.is_constant()) {
        // OT_LOGD("Min operation of two constant dists");
        return Distribution(value);
    } 
    else {
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

        Distribution temp(res_dist, value, st);

        temp.shrink();

        return temp;
    }
}

}
