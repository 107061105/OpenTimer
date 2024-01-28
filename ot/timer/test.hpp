#ifndef OT_TIMER_TEST_HPP_
#define OT_TIMER_TEST_HPP_

#include <ot/liberty/celllib.hpp>
#include <ot/timer/statisical.hpp>

namespace ot {

// Forward declaration
class Pin;
class Arc;

// ------------------------------------------------------------------------------------------------

// Class: Test
class Test {
  
  friend class Timer;
  friend class Endpoint;

  friend struct Path;

  public:
    
    Test(Arc&);

    std::optional<Dist > rat(Split, Tran) const;
    std::optional<Dist > slack(Split, Tran) const;
    std::optional<Dist > raw_slack(Split, Tran) const;
    std::optional<float> constraint(Split, Tran) const;
    std::optional<float> cppr_credit(Split, Tran) const;

    const Pin& constrained_pin() const;
    const Pin& related_pin() const;
    const Arc& arc() const;

  private:

    Arc& _arc;
    
    std::optional<std::list<Test>::iterator> _satellite;
    std::optional<std::list<Test*>::iterator> _pin_satellite;
    
    TimingData<std::optional<Dist >, MAX_SPLIT, MAX_TRAN> _rat;
    TimingData<std::optional<Dist >, MAX_SPLIT, MAX_TRAN> _related_at;
    TimingData<std::optional<float>, MAX_SPLIT, MAX_TRAN> _cppr_credit;
    TimingData<std::optional<float>, MAX_SPLIT, MAX_TRAN> _constraint;

    void _reset();
    // PathGuide
    void _fprop_rat(float, bool=false);
    
    Pin& _constrained_pin();
    Pin& _related_pin();
};



};  // end of namespace ot. -----------------------------------------------------------------------


#endif


