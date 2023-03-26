#include <ot/liberty/delay.hpp>
#include <cmath>

namespace ot {
    
// Constructor
Statisical_delay::Statisical_delay(float nom, float ms, float std, float skew) :
    _nominal {nom},
    _mean_shift {ms},
    _std {std},
    _skew {skew} {}

// Return nominal value of delay distribution
float Statisical_delay::nominal() const {
    return _nominal;
}

// Return mean shift of delay distribution
float Statisical_delay::meanshift() const {
    return _mean_shift;
}

// Return mean of delay distribution
float Statisical_delay::mean() const {
    return _nominal + _mean_shift;
}

// Return standard deviation of delay distribution
float Statisical_delay::stdev() const {
    return _std;
}

// Return skewness of delay distribution
float Statisical_delay::skew() const {
    return _skew;
}

// Operator
Statisical_delay Statisical_delay::operator* ( float s ) const {
    return Statisical_delay(_nominal * s, _mean_shift * s, _std * s, _skew * s); 
}

// Operator
Statisical_delay Statisical_delay::operator+ ( Statisical_delay dist ) const {
    float nominal = _nominal + dist.nominal();
    float meanshift = _mean_shift + dist.meanshift();
    float stdev = sqrt(_std * _std + dist.stdev() * dist.stdev());
    float skew = _skew + dist.skew();

    return Statisical_delay(nominal, meanshift, stdev, skew);
}

// Operator
Statisical_delay Statisical_delay::operator- ( Statisical_delay dist ) const {
    float nominal = _nominal - dist.nominal();
    float meanshift = _mean_shift - dist.meanshift();
    float stdev = sqrt(_std * _std + dist.stdev() * dist.stdev());
    float skew = _skew + dist.skew();

    return Statisical_delay(nominal, meanshift, stdev, skew);
}

// Operator
Statisical_delay Statisical_delay::operator+ ( float f ) const {
    return Statisical_delay(_nominal + f, _mean_shift, _std, _skew);
}

// Operator
Statisical_delay Statisical_delay::operator- ( float f ) const {
    return Statisical_delay(_nominal - f, _mean_shift, _std, _skew);
}

} // end of namespace ot. -----------------------------------------------------------------------

