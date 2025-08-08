#ifndef LIE_DETECTOR_HH
#define LIE_DETECTOR_HH

#include "ode_curve_observations.hh"

class Lie_detector
{
  
  public:

    Lie_Detector(ode_curve_observations &obs_) :
      obs(obs_)
    {}
    ~Lie_Detector() {}

    ode_curve_observations &obs;
};


#endif
