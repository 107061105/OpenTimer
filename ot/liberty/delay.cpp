#include <ot/liberty/delay.hpp>
#include <ot/static/logger.hpp>
#include <ot/headerdef.hpp>
#include <cmath>

namespace ot {
    
// Constructor
Statisical_delay::Statisical_delay(float nom, float ms, float std, float skew) :
    _nominal {nom},
    _mean_shift {ms},
    _std {std},
    _skew {skew} {}

// for debugging
void Statisical_delay::display() {
    OT_LOGD("Statisical delay\n"
            "nom:  ", nominal(), "\n",
            "ms:   ", meanshift(), "\n",
            "std:  ", stdev(), "\n",
            "skew: ", skew(), "\n\n");
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
    float skew = _skew - dist.skew();

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

