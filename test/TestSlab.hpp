#include <cxxtest/TestSuite.h>
#include <boost/assign.hpp>
#include <boost/shared_ptr.hpp>
#include "PetscSetupAndFinalize.hpp"
#include "BidomainProblem.hpp"
#include "HeartGeometryInformation.hpp"
#include "CardiacSimulationArchiver.hpp"
#include "HeartEventHandler.hpp"
#include "SingleTraceOutputModifier.hpp"
#include "SimpleStimulus.hpp"
#include "ORdGksVarierARVC05Small.hpp"
#include "PostUpstrokeTimeAdaptivityController.hpp"
#include "SlabCellFactory.hpp"
#include <string>
//#include "ConductivitiesModifierAnnulusIndexed3DExperimental.hpp"
class TestSolveTorso : public CxxTest::TestSuite
{

public:
    void TestSolve() throw(Exception)
    {




	double sim_duration;
        bool duration = CommandLineArguments::Instance()->OptionExists("--sim_duration");
        if (duration == true) {
            char* val = CommandLineArguments::Instance()->GetValueCorrespondingToOption("--sim_duration");
            sim_duration = atof(val);
        }
        std :: cout << "sim duration is " << sim_duration<< "ms \n";
        double cond_scaling;
        bool scaling = CommandLineArguments::Instance()->OptionExists("--cond_scaling");
        char* val;
 
	 if (scaling == true) {
//	    char* val;
            val = CommandLineArguments::Instance()->GetValueCorrespondingToOption("--cond_scaling");
            cond_scaling = atof(val);
        }
        std :: cout << "conductivity scaling factor is " << cond_scaling<< "ms \n";
	std::string cond_string(val);	
//Commandline arguments
/*   	std::vector <unsigned> nodes_present;
	unsigned dummy;
        bool root0 = CommandLineArguments::Instance()->OptionExists("--lv_basal_anterior_root");
        if (root0 == true) {
            char* val = CommandLineArguments::Instance()->GetValueCorrespondingToOption("--lv_basal_anterior_root");
            dummy = std::atoi(val);
	    nodes_present.push_back(dummy);
        }
	else
	{
	  std::cout << "lv_basal_anterior_root not specified in commandline args" <<std::endl;
	  exit(1);
	}

        bool root1 = CommandLineArguments::Instance()->OptionExists("--lv_apical_posterior_root");
        if (root1 == true) {
            char* val = CommandLineArguments::Instance()->GetValueCorrespondingToOption("--lv_apical_posterior_root");
            dummy = std::atoi(val);
	    nodes_present.push_back(dummy);
        }
	else
	{
	  std::cout << "lv_apical_posterior_root not specified in commandline args" <<std::endl;
	  exit(1);
	}
        bool root2 = CommandLineArguments::Instance()->OptionExists("--lv_mid_posterior_root");
        if (root2 == true) {
            char* val = CommandLineArguments::Instance()->GetValueCorrespondingToOption("--lv_mid_posterior_root");
            dummy = std::atoi(val);
	    nodes_present.push_back(dummy);
        }
	else
	{
	  std::cout << "lv_mid_posterior_root not specified in commandline args" <<std::endl;
	  exit(1);
	}
        bool root3 = CommandLineArguments::Instance()->OptionExists("--lv_septal_root");
        if (root3 == true) {
            char* val = CommandLineArguments::Instance()->GetValueCorrespondingToOption("--lv_septal_root");
            dummy = std::atoi(val);
	    nodes_present.push_back(dummy);
        }
	else
	{
	  std::cout << "lv_septal_root not specified in commandline args" <<std::endl;
	  exit(1);
	}
        bool root4 = CommandLineArguments::Instance()->OptionExists("--rv_mid_anterolateral_root");
        if (root4 == true) {
            char* val = CommandLineArguments::Instance()->GetValueCorrespondingToOption("--rv_mid_anterolateral_root");
            dummy = std::atoi(val);
	    nodes_present.push_back(dummy);
        }
	else
	{
	  std::cout << "rv_mid_anterolateral_root not specified in commandline args" <<std::endl;
	  exit(1);
	}
        bool root5 = CommandLineArguments::Instance()->OptionExists("--rv_basal_posterolateral_root");
        if (root5 == true) {
            char* val = CommandLineArguments::Instance()->GetValueCorrespondingToOption("--rv_basal_posterolateral_root");
            dummy = std::atoi(val);
	    nodes_present.push_back(dummy);
        }
	else
	{
	  std::cout << "rv_basal_posterolateral_root not specified in commandline args" <<std::endl;
	  exit(1);
	}

        bool root6 = CommandLineArguments::Instance()->OptionExists("--rv_septal_root");
        if (root6 == true) {
            char* val = CommandLineArguments::Instance()->GetValueCorrespondingToOption("--rv_septal_root");
            dummy = std::atoi(val);
	    nodes_present.push_back(dummy);
        }
	else
	{
	  std::cout << "rv_septal_root not specified in commandline args" <<std::endl;
	  exit(1)i;
	}

	int size_roots = nodes_present.size();
	for (int i=0; i<size_roots;i++)
	{
		std:: cout << nodes_present[i] <<std::endl;
	}*/
        
        /* Mesh to use */
        std::string filepath = "/gpfs/work/pn34qa/di39wun/meshes/mesh0.4mm/slab_0.4mm" ; 
        HeartConfig::Instance()->SetMeshFileName(filepath, cp::media_type::Orthotropic);

        // Add this to the IN_MESH time
        //HeartEventHandler::BeginEvent(HeartEventHandler::READ_MESH);
        DistributedTetrahedralMesh<3,3> mesh;
        TrianglesMeshReader<3,3> mesh_reader(filepath);
        mesh.ConstructFromMeshReader(mesh_reader);
       // HeartEventHandler::EndEvent(HeartEventHandler::READ_MESH);
       // HeartEventHandler::BeginEvent(HeartEventHandler::INITIALISE);

	 /* Construct HeartGeometryInformation */
        const std::string epi_face_file = filepath+".epi";
        const std::string rv_face_file = filepath+".rv";
        const std::string lv_face_file = filepath+".lv";
        //HeartGeometryInformation<3> heart_geom_info(mesh, epi_face_file, lv_face_file, rv_face_file, true);
	unsigned num_nodes = mesh.GetNumNodes();

        /* Simulation options */
        HeartConfig::Instance()->SetSimulationDuration(sim_duration); // ms
        std::string output_dir = "TestSlabSize0.4NormalCV"+cond_string;
        HeartConfig::Instance()->SetOutputDirectory(output_dir);
        HeartConfig::Instance()->SetOutputFilenamePrefix("results");
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.008, 0.008, 1.0);
        HeartConfig::Instance()->SetKSPPreconditioner("blockdiagonal");
	HeartConfig::Instance()->SetKSPSolver("symmlq");
        /* Conductivities */
        // Intra conductivities result in close to 67 30 17 cm/s conduction velocities as in pig
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(1.5*cond_scaling, 0.45*cond_scaling, 0.225*cond_scaling));
        // OL conductivity is 3.64x larger than IL as per Clerc.
        // OT conductivities (added in quadrature) 2.69x smaller than OL conductivity -> 2*0.91 and 0.91
        // Single OT conductivity (axisymmetric) -> 2.03 for both
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(3.64*1.5*cond_scaling, 2.03*cond_scaling, 2.03*cond_scaling));

        /* Defining the heart regions */
        std::set<unsigned> tissue_ids;
        tissue_ids.insert(0); // Same as default value defined in HeartConfig

        /* Defining the torso regions */
        std::set<unsigned> bath_ids;
        bath_ids.insert(1);
        //bath_ids.insert(2);
        //bath_ids.insert(3);
        //HeartConfig::Instance()->SetTissueAndBathIdentifiers(tissue_ids,bath_ids);

        /* Defining the torso conductivities in the corresponding regions (mS/cm) */
        std::map<unsigned, double> multiple_bath_conductivities;
        multiple_bath_conductivities[1] = 2.16;  // Bulk torso
        //multiple_bath_conductivities[2] = 0.2;   // Bones
        //multiple_bath_conductivities[3] = 0.389; // Lung
       // HeartConfig::Instance()->SetBathMultipleConductivities(multiple_bath_conductivities);


        //Modify CONDUCTIVITY PARAMETERS
       //read the coords of the 20 regions 
        std::vector<std::vector<double> >  regions(num_nodes,std::vector<double>(4,1));
   
  
    /* Visualisers */
        HeartConfig::Instance()->SetVisualizeWithMeshalyzer(false);
        HeartConfig::Instance()->SetVisualizeWithCmgui(false);
        HeartConfig::Instance()->SetVisualizeWithVtk(false);
        HeartConfig::Instance()->SetVisualizeWithParallelVtk(false);
        HeartConfig::Instance()->SetUseStateVariableInterpolation(false);

	
        /* Body surface traces */
        /*std::vector<unsigned> electrodes;
        // Load all surface nodes
        std::ifstream surface_nodes("projects/pm2111/test/data/nejib_20131010/surface_nodes.txt");
        std::string node;
        while ( surface_nodes >> node )
        {
            electrodes.push_back(std::atoi(node.c_str()));
        }
        if( !surface_nodes.eof() )
        {
            EXCEPTION("Couldn't read surface node file.");
        }
        HeartConfig::Instance()->SetRequestedNodalTimeTraces(electrodes);*/
	
        /* Output variables */
        /*std::vector<std::string> output_variables;
        output_variables.push_back("membrane_slow_delayed_rectifier_potassium_current_conductance");
        HeartConfig::Instance()->SetOutputVariables( output_variables );*/

        /* APD map */
        /*std::vector<std::pair<double,double> > apd_map;
        apd_map.push_back(std::pair<double, double>(90.0, 0.0)); // APD90, 0 mV threshold
        HeartConfig::Instance()->SetApdMaps(apd_map);*/

        /* Activation time map */
        std::vector<double> upstroke_time_map;
        upstroke_time_map.push_back(0.0); // 0 mV threshold
        HeartConfig::Instance()->SetUpstrokeTimeMaps(upstroke_time_map);

        // Gks varier (apicobasal, transmural, interventricular)

  
 std::vector<unsigned int> rootnodes;    // a vector to hold apexbase data
        
        std::vector<unsigned> lv_root_nodes;
        lv_root_nodes.push_back(8673u); // (5.577 9.896 6.501) LV anterior base
        /*lv_root_nodes.push_back(242233u);  // (2.230 7.500 3.506) LV posterior 1
        lv_root_nodes.push_back(258676u); // (2.180 6.441 4.740) LV posterior 2
        lv_root_nodes.push_back(237849u);  // (6.333 5.949 3.740) LV mid septum (slightly anterior)*/
        std::vector<unsigned> rv_root_nodes;
        rv_root_nodes.push_back(12661u); // (7.656 3.143 5.224) RV
        /*rv_root_nodes.push_back(399977u); // (4.938 0.704 6.102) RV 2
        rv_root_nodes.push_back(355733u); // (5.299 4.200 4.020) RV septum*/


  
        DijkstraStimuliCellFactory cell_factory( mesh.rGetNodePermutation()); // 150 cm/s

 
        // Problem
        //HeartEventHandler::EndEvent(HeartEventHandler::INITIALISE);
        // Problem immediately starts the EVERYTHING timer again
       // HeartEventHandler::EndEvent(HeartEventHandler::EVERYTHING);
        BidomainProblem<3> bidomain_problem( &cell_factory );
        //bidomain_problem.PrintOutput(false); // Turn OFF output

        // Adaptive time controller
        PostUpstrokeTimeAdaptivityController time_controller(0.008, 1.0, bidomain_problem);
        bidomain_problem.SetUseTimeAdaptivityController(true, &time_controller);

            std::vector<unsigned> perm_vec = mesh.rGetNodePermutation();
	//ConductivitiesModifierAnnulus3D conduction_modifier(mesh,regions,perm_vec);
        /* Modifier output */
        /*boost::shared_ptr<SingleTraceOutputModifier> trace_modifier_LA(new SingleTraceOutputModifier("trace_LA.txt", mesh.rGetNodePermutation()[9470282u]));
        bidomain_problem.AddOutputModifier(trace_modifier_LA);
        boost::shared_ptr<SingleTraceOutputModifier> trace_modifier_RA(new SingleTraceOutputModifier("trace_RA.txt", mesh.rGetNodePermutation()[9470311u]));
        bidomain_problem.AddOutputModifier(trace_modifier_RA);
        boost::shared_ptr<SingleTraceOutputModifier> trace_modifier_LL(new SingleTraceOutputModifier("trace_LL.txt", mesh.rGetNodePermutation()[9463801u]));
        bidomain_problem.AddOutputModifier(trace_modifier_LL);
        boost::shared_ptr<SingleTraceOutputModifier> trace_modifier_RL(new SingleTraceOutputModifier("trace_RL.txt", mesh.rGetNodePermutation()[9463839u]));
        bidomain_problem.AddOutputModifier(trace_modifier_RL);
        boost::shared_ptr<SingleTraceOutputModifier> trace_modifier_V1(new SingleTraceOutputModifier("trace_V1.txt", mesh.rGetNodePermutation()[9467954u]));
        bidomain_problem.AddOutputModifier(trace_modifier_V1);
        boost::shared_ptr<SingleTraceOutputModifier> trace_modifier_V2(new SingleTraceOutputModifier("trace_V2.txt", mesh.rGetNodePermutation()[9467789u]));
        bidomain_problem.AddOutputModifier(trace_modifier_V2);
        boost::shared_ptr<SingleTraceOutputModifier> trace_modifier_V3(new SingleTraceOutputModifier("trace_V3.txt", mesh.rGetNodePermutation()[9467331u]));
        bidomain_problem.AddOutputModifier(trace_modifier_V3);
        boost::shared_ptr<SingleTraceOutputModifier> trace_modifier_V4(new SingleTraceOutputModifier("trace_V4.txt", mesh.rGetNodePermutation()[9466470u]));
        bidomain_problem.AddOutputModifier(trace_modifier_V4);
       
        boost::shared_ptr<SingleTraceOutputModifier> trace_modifier_V5(new SingleTraceOutputModifier("trace_V5.txt", mesh.rGetNodePermutation()[9466680u]));
        bidomain_problem.AddOutputModifier(trace_modifier_V5);
        boost::shared_ptr<SingleTraceOutputModifier> trace_modifier_V6(new SingleTraceOutputModifier("trace_V6.txt", mesh.rGetNodePermutation()[9466972u]));
        bidomain_problem.AddOutputModifier(trace_modifier_V6);
	*/
        bidomain_problem.SetMesh( &mesh );
	std::cout << "I load the mesh" << std::endl;
        bidomain_problem.SetWriteInfo();
	std::cout << "I s et write info" << std::endl;
        //bidomain_problem.SetUseHdf5DataWriterCache(true);

        bidomain_problem.Initialise();
//	std::cout << "I initialise the problem" << std::endl;
      // Get the tissue from the problem and apply the modifier
//	 BidomainTissue<3>* p_bidomain_tissue = bidomain_problem.GetBidomainTissue();
//	std::cout <<" i get bidomain tissue" << std::endl;
//	p_bidomain_tissue->SetConductivityModifier( &conduction_modifier );
//	std::cout <<" i set conductivities via modifier" << std::endl;

        /* Ground the nodes of an element containing the (standard) right leg electrode */
       // std::vector<unsigned> grounded_nodes;
/*for patient 02 i m unsure which triangular element is closest to surface electrode so i ground 2 of them*/
        /*grounded_nodes.push_back(mesh.rGetNodePermutation()[9463931u]); //in chaste numbering
        grounded_nodes.push_back(mesh.rGetNodePermutation()[9463825u]);
        grounded_nodes.push_back(mesh.rGetNodePermutation()[9463838u]);
        grounded_nodes.push_back(mesh.rGetNodePermutation()[9463838u]);
        grounded_nodes.push_back(mesh.rGetNodePermutation()[9463732u]);
        grounded_nodes.push_back(mesh.rGetNodePermutation()[9463748u]);*/
	std::cout << "I select grounded nodes" << std::endl;
	
        //bidomain_problem.SetFixedExtracellularPotentialNodes(grounded_nodes);

 
        // Solve call -> the solver -> the assembler -> modified conductivities

	
        bidomain_problem.Solve();
	std::cout << "I solve the problem" << std::endl;
        /* Archive */
//        CardiacSimulationArchiver<BidomainProblem<3> >::Save(bidomain_problem, output_dir+"/saved_simulation");

    //    HeartEventHandler::Headings();
  //      HeartEventHandler::Report();
//        HeartEventHandler::Headings();

        /* Save permutation just in case */
        OutputFileHandler output_file_handler(output_dir, false); // collective
        if ( PetscTools::AmMaster() )
        {
            out_stream perm_out_stream = output_file_handler.OpenOutputFile("permutation.txt", std::ios::out);

            for (unsigned i=0; i<perm_vec.size(); ++i)
            {
                (*perm_out_stream) << perm_vec[i] << "\n";
            }
        }

//	 BidomainProblem<3>* p_bidomain_problem = CardiacSimulationArchiver<BidomainProblem<3> >::Load(output_dir+"/saved_simulation");
	
//	 HeartConfig::Instance()->SetSimulationDuration(10); //ms
  //       HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.02, 0.02, 1.0);
    //     p_bidomain_problem->Solve();

      //   delete p_bidomain_problem;

    }
};
