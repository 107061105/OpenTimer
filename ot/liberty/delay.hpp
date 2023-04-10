#ifndef OT_LIBERTY_DELAY_HPP_
#define OT_LIBERTY_DELAY_HPP_

#include <ot/headerdef.hpp>

namespace ot {

// Class Statisical_delay
class Statisical_delay {

    public:

        Statisical_delay(float, float, float, float);

        inline float nominal() const;
        inline float meanshift() const;
        inline float mean() const;
        inline float stdev() const;
        inline float skew() const;
        inline float max() const;
        inline float min() const;

        // for debugging
        void display();

        Statisical_delay operator* ( float s ) const ;
        Statisical_delay operator+ ( float v ) const ;
        Statisical_delay operator- ( float v ) const ;
        Statisical_delay operator+ ( Statisical_delay dist ) const ;
        Statisical_delay operator- ( Statisical_delay dist ) const ;

    private:
 
        float _nominal;
        float _mean_shift;
        float _std;
        float _skew;

};

// ------------------------------------------------------------------------------------------------


// Return nominal value of delay distribution
inline float Statisical_delay::nominal() const {
  return _nominal;
}

// Return mean shift of delay distribution
inline float Statisical_delay::meanshift() const {
  return _mean_shift;
}

// Return mean of delay distribution
inline float Statisical_delay::mean() const {
  return _nominal + _mean_shift;
}

// Return standard deviation of delay distribution
inline float Statisical_delay::stdev() const {
  return _std;
}

// Return skewness of delay distribution
inline float Statisical_delay::skew() const {
  return _skew;
}

inline float Statisical_delay::max() const {
  return _nominal + _mean_shift + 3 * _std;
}

inline float Statisical_delay::min() const {
  return _nominal + _mean_shift - 3 * _std;
}

// ------------------------------------------------------------------------------------------------

};  // end of namespace ot. -----------------------------------------------------------------------

#endif
