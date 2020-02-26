#include <cxxtest/TestSuite.h>
#include "CellProperties.hpp"
#include "SteadyStateRunner.hpp"
#include "AbstractCvodeCell.hpp"
#include "RegularStimulus.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "Shannon2004Cvode.hpp"
#include "FakePetscSetup.hpp"

class TestSingleCellSimulationTutorial : public CxxTest::TestSuite
{
    public:
        void TestShannonSimulation()
        {

#ifdef CHASTE_CVODE
            boost::shared_ptr<RegularStimulus> p_stimulus(new RegularStimulus(-40.0,1.0,1000.0,0.0));//params are magnitude, duration, period, and start time of stimulus.
            boost::shared_ptr<AbstractIvpOdeSolver> p_solver ;
            boost::shared_ptr<AbstractCvodeCell> p_model(new CellShannon2004FromCellMLCvode(p_solver,p_stimulus));
            
            boost::shared_ptr<RegularStimulus> p_regular_stim = p_model->UseCellMLDefaultStimulus();
            
            double max_timestep = p_regular_stim->GetDuration();
            
            p_model->SetMaxSteps(1e6);
            
            SteadyStateRunner steady_runner(p_model);
            bool result;
            result = steady_runner.RunToSteadyState();
            TS_ASSERT_EQUALS(result,false); //check that the model has NOT reached steady state 
            
            p_model->SetMaxTimestep(max_timestep); 
            
            double sampling_timestep = max_timestep;
            double start_time = 0.0;
            double end_time = 10000.0; //this is in ms i assume
            OdeSolution solution = p_model->Compute(start_time,end_time,sampling_timestep);
            
            solution.WriteToFile("TestCvodeCells","Shannon2004Cvode","ms");
        
        
#else
            std::cout << "Cvode is not enabled.\n";
            
#endif

        }
        
};
