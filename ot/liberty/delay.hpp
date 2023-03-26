#ifndef OT_LIBERTY_DELAY_HPP_
#define OT_LIBERTY_DELAY_HPP_

namespace ot {

// Class Statisical_delay
class Statisical_delay {

    public:

        Statisical_delay(float, float, float, float);

        float nominal() const;
        float meanshift() const;
        float mean() const;
        float stdev() const;
        float skew() const;

        Statisical_delay operator* ( float s ) const ;
        Statisical_delay operator+ ( Statisical_delay dist ) const ;
        Statisical_delay operator- ( Statisical_delay dist ) const ;
        Statisical_delay operator+ ( float v ) const ;
        Statisical_delay operator- ( float v ) const ;

        operator float () {
          return _nominal;
        }

    private:
 
        float _nominal;
        float _mean_shift;
        float _std;
        float _skew;

};


};  // end of namespace ot. -----------------------------------------------------------------------

#endif
