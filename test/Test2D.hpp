#include <cxxtest/TestSuite.h>
#include "BidomainProblem.hpp"
#include "RegularStimulus.hpp"
#include "PetscSetupAndFinalize.hpp"
#include <Hdf5DataWriter.hpp>
#include "PostUpstrokeTimeAdaptivityControllerCL2D.hpp"
#include "PostProcessingWriter.hpp"
#include "ToRORd_fkatp_endoCvodeOpt.hpp"
#include "ORd2011endo_fkatpCvodeOpt.hpp"
#include "TetrahedralMesh.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"
#include "ConductivitiesModifierEndo2D.hpp"
#include <algorithm>    // std::find
#include <cmath>        // std::abs
#include "cell_based/src/common/SimulationTime.hpp"
#include <vector>
class PointStimulus2dCellFactory : public AbstractCardiacCellFactory<2>
{
private:


public:
    PointStimulus2dCellFactory()
        : AbstractCardiacCellFactory<2>()
     
    {
    }

	   AbstractCvodeCell* CreateCardiacCellForTissueNode(Node<2>* pNode)
    {
        double x = pNode->rGetLocation()[0];
        double y = pNode->rGetLocation()[1];
 	   AbstractCvodeCell* p_cell;
	boost::shared_ptr<RegularStimulus> p_stim(new RegularStimulus(-200000.0, 1.0, 1000,0.5));
        if (x<0.02+1e-6 && y<0.02+1e-6) // ie if x<=0.02 and y<=0.02 (and we are assuming here x,y>=0).
        { 

	   p_cell = new CellToRORd_fkatp_endoFromCellMLCvodeOpt(mpSolver, mpZeroStimulus);
	   p_cell->SetStimulusFunction(p_stim);
           return p_cell;
        }
        else
        {
 	   AbstractCvodeCell* p_cell;
	   p_cell = new CellToRORd_fkatp_endoFromCellMLCvodeOpt(mpSolver, mpZeroStimulus);
           return p_cell;	
     
        }
    }


};



class TestRunningBidomainSimulationsTutorial : public CxxTest::TestSuite
{
public:
    void TestSimpleSimulation()
    {
      //  HeartEventHandler::BeginEvent(HeartEventHandler::EVERYTHING);

        HeartConfig::Instance()->SetSimulationDuration(3000.0); //ms
	std::string filepath = "mesh/test/data/2D_0_to_1mm_800_elements";
   //     HeartConfig::Instance()->SetMeshFileName(filepath);
       HeartConfig::Instance()->SetMeshFileName(filepath);

       std:: cout << "i get to specify mesh name" << std::endl;
	std::string output_folder = "BidomainTutorial";
        HeartConfig::Instance()->SetOutputDirectory(output_folder);
        HeartConfig::Instance()->SetOutputFilenamePrefix("results");

        HeartConfig::Instance()->SetVisualizeWithMeshalyzer(false);
        HeartConfig::Instance()->SetVisualizeWithCmgui(false);
        HeartConfig::Instance()->SetVisualizeWithVtk(false);
        //HeartConfig::Instance()->SetVisualizeWithParallelVtk(true);
        // HeartEventHandler::EndEvent(HeartEventHandler::INITIALISE);
        // Problem immediately starts the EVERYTHING timer again
       // HeartEventHandler::EndEvent(HeartEventHandler::EVERYTHING);
    // Add this to the IN_MESH time
//       HeartConfig::Instance()->SetMeshFileName(filepath, cp::media_type::NoFibreOrientation);

	 // HeartEventHandler::BeginEvent(HeartEventHandler::READ_MESH);
     //   DisributedTetrahedralMesh<2,2> mesh;
       // TrianglesMeshReader<2,2> mesh_reader(filepath);
         //mesh.ConstructFromMeshReader(mesh_reader);
       TrianglesMeshReader<2,2> reader(filepath);
      DistributedTetrahedralMesh<2,2> mesh;
       mesh.ConstructFromMeshReader(reader);
//	HeartEventHandler::EndEvent(HeartEventHandler::READ_MESH);
       std:: cout << "i get to cell fac" << std::endl;

        PointStimulus2dCellFactory cell_factory;
       std:: cout << "i get past cell fac" << std::endl;

        BidomainProblem<2> bidomain_problem( &cell_factory );
       std:: cout << "i create bid prob" << std::endl;

        //rtEventHandler::EndEvent(HeartEventHandler::INITIALISE);
        // Problem immediately starts the EVERYTHING timer again
        
	
        // bidomain_problem.Solve();
	double scale_cond=1.0;
	c_vector<double,2> extracellular_conductivities=Create_c_vector(3.64*1.5*scale_cond, 2.03*scale_cond);
        c_vector<double,2> intracellular_conductivities=Create_c_vector(1.5, 0.45);

	c_vector<double,2> monodomain_conductivities;
	
            for (unsigned dim=0; dim<2; dim++)
            {
                monodomain_conductivities[dim] = intracellular_conductivities[dim]*extracellular_conductivities[dim]
                                               / (intracellular_conductivities[dim] + extracellular_conductivities[dim]);
            }

	HeartConfig::Instance()->SetIntracellularConductivities(intracellular_conductivities);

        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1400); // 1/cm
        HeartConfig::Instance()->SetCapacitance(1.0); // uF/cm^2

        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.01, 1);


       std::vector<double> upstroke_time_map;
       upstroke_time_map.push_back(0.0); // 0 mV threshold
      HeartConfig::Instance()->SetUpstrokeTimeMaps(upstroke_time_map);

	std:: cout << "i get to time adap controller" << std::endl;
        PostUpstrokeTimeAdaptivityControllerCL2D time_controller(0.01, 0.25, bidomain_problem, 1000);
        bidomain_problem.SetUseTimeAdaptivityController(true, &time_controller);


        bidomain_problem.SetWriteInfo();
	
        bidomain_problem.SetUseHdf5DataWriterCache(true);

        bidomain_problem.Initialise();
	std::vector<unsigned> perm_vec = mesh.rGetNodePermutation();
//Modify CONDUCTIVITY PARAMETERS
       //read the coords of the fibrotic nodes
	
        
	 std::vector<unsigned int>  fibrosis_ids;
   
                std::ifstream infile;
                const std::string filename = "/hppfs/work/pn34qa/di39wun2/meshes/2D/fibrotic_ids.txt";

                infile.open(filename.c_str());
                if (!infile) {
                    std::cout << "Unable to open file datafile.txt";
                    exit(1);   // call system to stop
                }

		//assign the order of the region in the regions array


            while(!infile.eof())
            {

                     int a;
                     infile >> a ;
		     //std::cout << a<< " "  << b  << "  "<<  c << "  " << d << std::endl;
                     fibrosis_ids.push_back(a);



                     //std::cout << i <<std::endl;
            }   
		 infile.close();

	double fibrosis=0;

	if (fibrosis==1)
		{
			std::cout <<" i enter the condctivities modifier" << std::endl;





			Cond2 conduction_modifier(mesh,fibrosis_ids,perm_vec);

			std::cout <<" i finish with the condctivities modifier" << std::endl;

			BidomainTissue<2>* p_bidomain_tissue = bidomain_problem.GetBidomainTissue();
			std::cout <<" i get bidomain tissue" << std::endl;
			p_bidomain_tissue->SetConductivityModifier( &conduction_modifier );
			std::cout <<" i set conductivities via modifier" << std::endl;
		}
        bidomain_problem.Solve();

        ReplicatableVector res_repl(bidomain_problem.GetSolution());
        for (unsigned i=0; i<res_repl.GetSize(); i++)
        {
        //    std::cout << res_repl[i] << "\n";
        }

        HeartEventHandler::Headings();
        HeartEventHandler::Report();
        FileFinder data_folder(output_folder, RelativeTo::ChasteTestOutput);
       /* Do postprocessing (in a block to isolate HDF5 writers */
 /*
142	         * APD AND ACTIVATION MAPS
	         */
	
	        std::vector<std::pair<double,double> > apd_maps;
	        apd_maps.push_back(std::pair<double, double>(90,0)); //90 %, -30 mV
	        HeartConfig::Instance()->SetApdMaps(apd_maps);
		
	{
	            PostProcessingWriter<2,2> ppw(mesh, data_folder, "results");
	            ppw.WritePostProcessingFiles();
       }
    }

};
