#ifndef OT_LIBERTY_DELAY_HPP_
#define OT_LIBERTY_DELAY_HPP_

// #include <ot/headerdef.hpp>
#include <iostream>
#include <vector>
#include "stats_util.hpp"

typedef std::complex<float> Complex;

// The data structure of the statisical analysis,
// used on statisical static analysis
namespace Statisical 
{

class Distribution {
public:
    // Constructor
    Distribution(float);
    Distribution(const std::vector<float>&, int);
    Distribution(Distribution_type, float, float, float);
    // Destructor
    ~Distribution() { _pdf.clear(); }

    // Summation operation of the distributions
    Distribution operator+(float);
    Distribution operator+(Distribution &);
    // Subtraction operation of the distributions
    Distribution operator-(float);
    Distribution operator-(Distribution &);

    // Print the status
    void print_status();

    // Distribution operation functions
    void norm();
    void negate();
    void shrink();

    // Query functions
    inline Distribution_type get_type () const { return _type; }
    inline float get_start_time       () const { return _start * TIME_STEP; }
    inline float get_end_time         () const { return (_start + static_cast<int>(_pdf.size())) * TIME_STEP; }
    inline int get_start_point        () const { return _start; }
    inline int get_end_point          () const { return _start + static_cast<int>(_pdf.size()); }
    inline int get_bin_num            () const { return static_cast<int>(_pdf.size()); }
    const std::vector<float>& get_pdf () const { return _pdf; }
    float get_ith_pdf(int i) const;
    float get_3_sigma(Split);

private:
    // Type of distribution
    Distribution_type _type;
    // The probability density function of the distribution
    std::vector<float> _pdf;
    // The starting/end time of the distribution
    // Start time = _start * TIME_STEP
    // End time = (_start + bin number) * TIME_STEP
    int _start;
};

Distribution max(const Distribution&, const Distribution&);
Distribution min(const Distribution&, const Distribution&);

}

#endif

