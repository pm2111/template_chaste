#ifndef POSTUPSTROKETIMEADAPTIVITYCONTROLLER_HPP_
#define POSTUPSTROKETIMEADAPTIVITYCONTROLLER_HPP_

#include "AbstractTimeAdaptivityController.hpp"

#include "BidomainProblem.hpp" ///\todo: template over problem?

class PostUpstrokeTimeAdaptivityController : public AbstractTimeAdaptivityController
{
private:
    double mAllDepolarisedTime;
    bool mHaveSwitchedTimestep;
    double mSmallDt;
    double mLargeDt;
    BidomainProblem<3>& mrProblem; // For tissue, for cells

public:
    /*
     *  Constructor.
     */
    PostUpstrokeTimeAdaptivityController(double smallDt, double largeDt, BidomainProblem<3>& rProblem)
        : AbstractTimeAdaptivityController(smallDt, largeDt),
          mAllDepolarisedTime(-1.0), // Nonsense value to indicate not all depolarised
          mHaveSwitchedTimestep(false),
          mSmallDt(smallDt),
          mLargeDt(largeDt),
          mrProblem(rProblem)
    {
    }
    /*
     *  The method that returns the time step
     */
    double ComputeTimeStep(double currentTime, Vec currentSolution);
};

#endif //POSTUPSTROKETIMEADAPTIVITYCONTROLLER_HPP_
