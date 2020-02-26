#ifndef POSTUPSTROKETIMEADAPTIVITYCONTROLLERMONO_HPP_
#define POSTUPSTROKETIMEADAPTIVITYCONTROLLERMONO_HPP_

#include "AbstractTimeAdaptivityController.hpp"

#include "MonodomainProblem.hpp" ///\todo: template over problem?

class PostUpstrokeTimeAdaptivityControllerCLMono : public AbstractTimeAdaptivityController
{
private:
    double mAllDepolarisedTime;
    bool mHaveSwitchedTimestep;
    double mSmallDt;
    double mLargeDt;
    MonodomainProblem<3>& mrProblem; // For tissue, for cells
    double mPacingFrequency; /**< Time period at which to resume small steps. Default is 0 (never reset) */

public:
    /*
     *  Constructor.
     */
    PostUpstrokeTimeAdaptivityControllerCLMono(double smallDt, double largeDt, MonodomainProblem<3>& rProblem, double pacingFrequency=0.0)
        : AbstractTimeAdaptivityController(smallDt, largeDt),
          mAllDepolarisedTime(-1.0), // Nonsense value to indicate not all depolarised
          mHaveSwitchedTimestep(false),
          mSmallDt(smallDt),
          mLargeDt(largeDt),
          mrProblem(rProblem),
          mPacingFrequency(pacingFrequency)
	  {
	  }



	
  
    /*
     *  The method that returns the time step
     */
    double ComputeTimeStep(double currentTime, Vec currentSolution);
};

#endif //POSTUPSTROKETIMEADAPTIVITYCONTROLLER_HPP_
