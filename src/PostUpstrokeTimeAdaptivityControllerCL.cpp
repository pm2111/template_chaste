#include "PostUpstrokeTimeAdaptivityControllerCL.hpp"

double PostUpstrokeTimeAdaptivityControllerCL::ComputeTimeStep(double currentTime, Vec currentSolution)
{
    if (mAllDepolarisedTime<0.0)
    {
        // Check if all cells are depolarised. This tells us we can switch to a longer time step.
        double min_voltage;
        VecStrideMin( currentSolution, 0, PETSC_NULL, &min_voltage );
        if (min_voltage < 0.0)
        {
            return mSmallDt; // Still got an upstroke
        }
        else
        {
            mAllDepolarisedTime = currentTime; // Won't check this any more, continue
        }
    }
    // Everything is depolarised, at the next opportunity switch to the longer timestep
    if (!mHaveSwitchedTimestep)
    {
        // Only switch over 5 ms after the last zero-crossing AND on a multiple of the print time (hardcoded!)
        if ( currentTime > mAllDepolarisedTime+5.0 && Divides(1.0, currentTime) )
        {
            // Set dt for all cells to large time step
            const std::vector<AbstractCardiacCellInterface*>& r_cells = mrProblem.GetTissue()->rGetCellsDistributed();
            for (unsigned i=0; i<r_cells.size(); i++)
            {
                r_cells[i]->SetTimestep(mLargeDt);
            }
            mHaveSwitchedTimestep = true; // Won't check this any more, continue
        }
        else
        {
            return mSmallDt;
        }
    }
    // Reset when time divides mPacingFrequency (up to half a mSmallDt)
    if (CompareDoubles::IsNearZero(fmod(currentTime, mPacingFrequency), mSmallDt/2))
    {
        // Set flags to initial states
        mAllDepolarisedTime = -1.0;
        mHaveSwitchedTimestep = false;
        // Set dt for all cells to small time step
        const std::vector<AbstractCardiacCellInterface*>& r_cells = mrProblem.GetTissue()->rGetCellsDistributed();
        for (unsigned i=0; i<r_cells.size(); i++)
        {
            r_cells[i]->SetTimestep(mSmallDt);
        }
        return mSmallDt;
    }
    return mLargeDt;
}
