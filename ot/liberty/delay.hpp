#ifndef OT_LIBERTY_DELAY_HPP_
#define OT_LIBERTY_DELAY_HPP_

#include <ot/liberty/stats_util.hpp>
#include <ot/headerdef.hpp>
#include <iostream>
#include <optional>
#include <vector>
#include <cassert>

typedef std::complex<float> Complex;

// The data structure of the statisical analysis,
// used on statisical static analysis
namespace Statisical 
{

class Distribution {
public:
    // Constructor
    Distribution(float);
    Distribution(std::vector<float>&, float, int);
    Distribution(Distribution_type, ot::Tran, float);
    // Distribution(Distribution_type, float, float, float);
    // Destructor
    ~Distribution() { _pdf.clear(); _start.reset(); }

    // Summation operation of the distributions
    Distribution operator+(float);
    Distribution operator+(const Distribution &);
    // Subtraction operation of the distributions
    Distribution operator-(float);
    Distribution operator-(const Distribution &);
    // Comparation operation of the CONSTANT distributions
    bool operator<(const Distribution &);
    bool operator>(const Distribution &);

    // Print the status
    void print_status();

    // Distribution operation functions
    void norm();
    void negate();
    void shrink();

    // Query functions
    inline Distribution_type get_type () const;
    inline float get_start_time       () const;
    inline float get_end_time         () const;
    inline float get_value            () const;
    inline int   get_start_point      () const;
    inline int   get_end_point        () const;
    inline int   get_bin_num          () const;
    inline bool  is_constant          () const;
    inline const std::vector<float>& get_pdf () const;
    float get_ith_pdf(int i) const;
    float get_3_sigma(ot::Split);

private:
    // Type of distribution
    Distribution_type _type;
    // Mean value of the distribution
    float _value;
    // The probability density function of the distribution
    std::vector<float> _pdf;
    // The starting/end time of the distribution
    // Start time = _start * TIME_STEP
    // End time = (_start + bin number) * TIME_STEP
    std::optional<int> _start;
};

// Function: get_type
inline Distribution_type Distribution::get_type() const { 
    return _type; 
}

// Function: get_start_time
inline float Distribution::get_start_time() const {
    assert(_start);
    return get_start_point() * TIME_STEP; 
}

// Function: get_end_time
inline float Distribution::get_end_time() const { 
    assert(_start);
    return get_end_point() * TIME_STEP; 
}

// Function: get_value
inline float Distribution::get_value() const { 
    return _value;
}

// Function: 
inline int Distribution::get_start_point() const { 
    assert(_start);
    return *_start; 
}

// Function: 
inline int Distribution::get_end_point() const {
    assert(_start);
    return *_start + static_cast<int>(_pdf.size()); 
}

// Function: 
inline int Distribution::get_bin_num() const { 
    return static_cast<int>(_pdf.size()); 
}

// Function: 
inline bool Distribution::is_constant() const { 
    if (_type == Distribution_type::Constant) {
        assert(!_start && _pdf.empty());
        return true;
    }
    return false; 
}
    
inline const std::vector<float>& Distribution::get_pdf() const { 
    return _pdf; 
}

Distribution max(const Distribution&, const Distribution&);
Distribution min(const Distribution&, const Distribution&);

}

#endif

