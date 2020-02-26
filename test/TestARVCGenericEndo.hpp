#include <cxxtest/TestSuite.h>

#include "TetrahedralMesh.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "BidomainProblem.hpp"
#include "HeartGeometryInformation.hpp"
#include "CardiacSimulationArchiver.hpp"
#include "HeartEventHandler.hpp"
#include "SingleTraceOutputModifier.hpp"
#include "RegularStimulus.hpp"
#include "ORdGksVarierARVC018.hpp"
#include "PostUpstrokeTimeAdaptivityControllerCL.hpp"
#include "DijkstraARVC018.hpp"
#include "PseudoEcgCalculator.hpp"
#include "OutputFileHandler.hpp"
#include "FileComparison.hpp"
#include "HeartConfig.hpp"
#include <Hdf5DataWriter.hpp>
#include "ConductivitiesModifierEndo.hpp"
#include <algorithm> //for search of root node vec
#include <string>     // std::string, std::stoi

class TestSolveTorso : public CxxTest::TestSuite
{

public:
    void TestSolve() throw(Exception)
    {

	//Commandline arguments
   	std::vector <unsigned> nodes_present;
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
	  exit(1);
	}

       bool reg = CommandLineArguments::Instance()->OptionExists("--regions_file");
	std::string regions_file;
        if (reg == true) {
	
            regions_file = CommandLineArguments::Instance()->GetStringCorrespondingToOption("--regions_file");

        }
	else
	{
	  std::cout << "regions_file not specified in commandline args" <<std::endl;
	  exit(1);
	}


	std::string output_dir;
        output_dir = CommandLineArguments::Instance()->GetStringCorrespondingToOption("--save_dir");
        std::string apexbase_path;
        apexbase_path = CommandLineArguments::Instance()->GetStringCorrespondingToOption("--apexbase");

        std::string mindists_path;
        mindists_path = CommandLineArguments::Instance()->GetStringCorrespondingToOption("--mindists");




        std::string edge_nodes_path;
        edge_nodes_path = CommandLineArguments::Instance()->GetStringCorrespondingToOption("--edge_nodes");
	    int fibrosis;
            char* val = CommandLineArguments::Instance()->GetValueCorrespondingToOption("--fibrosis");
            fibrosis = std::atoi(val);
           double scale_cond;
            char* val2 = CommandLineArguments::Instance()->GetValueCorrespondingToOption("--scale_cond");
            scale_cond = std::atoi(val2);
           double time_step=0.008;
        //    char* val3 = CommandLineArguments::Instance()->GetValueCorrespondingToOption("--pde_timestep");
          //  time_step = std::atoi(val3);



	

        HeartEventHandler::BeginEvent(HeartEventHandler::EVERYTHING);

        /* Mesh to use */
        std::string filepath = edge_nodes_path+"HEART" ; 
        HeartConfig::Instance()->SetMeshFileName(filepath, cp::media_type::Orthotropic);


        // Add this to the IN_MESH time
        HeartEventHandler::BeginEvent(HeartEventHandler::READ_MESH);
        DistributedTetrahedralMesh<3,3> mesh;
        TrianglesMeshReader<3,3> mesh_reader(filepath);
        mesh.ConstructFromMeshReader(mesh_reader);
        HeartEventHandler::EndEvent(HeartEventHandler::READ_MESH);
        HeartEventHandler::BeginEvent(HeartEventHandler::INITIALISE);
	unsigned num_nodes = mesh.GetNumNodes();
	 /* Construct HeartGeometryInformation */
        const std::string epi_face_file = filepath+".epi";
        const std::string rv_face_file = filepath+".rv";
        const std::string lv_face_file = filepath+".lv";
        HeartGeometryInformation<3> heart_geom_info(mesh, epi_face_file, lv_face_file, rv_face_file, true);


        /* Simulation options */
        HeartConfig::Instance()->SetSimulationDuration(1400); // ms
	double period = 800.0;//stimulation period in ms, 
        HeartConfig::Instance()->SetOutputDirectory(output_dir);
        HeartConfig::Instance()->SetOutputFilenamePrefix("results");
	//double time_step=time_step;
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(time_step, time_step, 1.0);
        HeartConfig::Instance()->SetKSPPreconditioner("blockdiagonal");
	HeartConfig::Instance()->SetKSPSolver("symmlq");
        /* Conductivities */
        // Intra conductivities result in close to 67 30 17 cm/s conduction velocities as in pig

	//find the equivalent conductivities for the bidomain model to reduce to the monodomain

	c_vector<double,3> extracellular_conductivities=Create_c_vector(3.64*1.5*scale_cond, 2.03*scale_cond, 2.03*scale_cond);
        c_vector<double,3> intracellular_conductivities=Create_c_vector(1.5, 0.45, 0.225);

	c_vector<double,3> monodomain_conductivities;
	
            for (unsigned dim=0; dim<3; dim++)
            {
                monodomain_conductivities[dim] = intracellular_conductivities[dim]*extracellular_conductivities[dim]
                                               / (intracellular_conductivities[dim] + extracellular_conductivities[dim]);
            }

	HeartConfig::Instance()->SetIntracellularConductivities(intracellular_conductivities);
        // OL conductivity is 3.64x larger than IL as per Clerc.
        // OT conductivities (added in quadrature) 2.69x smaller than OL conductivity -> 2*0.91 and 0.91
        // Single OT conductivity (axisymmetric) -> 2.03 for both
        HeartConfig::Instance()->SetExtracellularConductivities(extracellular_conductivities);

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
	
        std::vector<std::vector<double> >  regions(num_nodes,std::vector<double>(2,1));
   
                std::ifstream infile;
                const std::string filename = regions_file;

                infile.open(filename.c_str());
                if (!infile) {
                    std::cout << "Unable to open file datafile.txt";
                    exit(1);   // call system to stop
                }

		//assign the order of the region in the regions array


            for(int i=0; i<num_nodes;i++)
            {

                     double a,b;
                     infile >> a >> b;
		     //std::cout << a<< " "  << b  << "  "<<  c << "  " << d << std::endl;
                     regions[i][0]=a;
                     regions[i][1]=b;


                     //std::cout << i <<std::endl;
            }   
		 infile.close();



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

        ORdGksVarierDTI024 gks_varier(&heart_geom_info, true, false, false,apexbase_path);  //transmural heterogeneities are defined in the Dijkstra file

	/* Root node specification*/
	
    	std::vector<unsigned int> rootnodes;    // a vector to hold apexbase data
	// open file    
	std::ifstream inputFile(edge_nodes_path+"roots/start.txt");

	if (inputFile) {       
	int value;
	
	// read the elements in the file into a vector  
	while ( inputFile >> value ) {
            unsigned int value2 = value;
	rootnodes.push_back(value2);
	}
	}
  	std::vector<unsigned int> rootnodes_permuted;    // a vector to hold apexbase data
	// open file    
	std::ifstream inputFile1(edge_nodes_path+"roots/permuted.txt");
	inputFile.close();
	if (inputFile1) {       
	int value;
	
	// read the elements in the file into a vector  
	while ( inputFile1 >> value ) {
            unsigned int value2 = value;
	rootnodes_permuted.push_back(value2);
	}
	}
	inputFile1.close();
       // LV and RV activation root nodes
        std::vector<unsigned> lv_root_nodes;
	for (int i =0; i<4;i++)
{
	if ( nodes_present[i]==1 )
	{
      	  lv_root_nodes.push_back(rootnodes[i]); // (5.577 9.896 6.501) LV anterior base
	}
	  if( nodes_present[i] ==2)
	{
		 lv_root_nodes.push_back(rootnodes_permuted[i]); // (5.577 9.896 6.501) LV anterior base
	}  
	
}

        std::vector<unsigned> rv_root_nodes;
	for (int i =4; i<7;i++)
{
	if ( nodes_present[i]==1 )
	{
      	  rv_root_nodes.push_back(rootnodes[i]); // (5.577 9.896 6.501) LV anterior base
	}
	  if( nodes_present[i] ==2)
	{
		 rv_root_nodes.push_back(rootnodes_permuted[i]); // (5.577 9.896 6.501) LV anterior base
	}  

}
	



       // std::vector<unsigned> lv_root_nodes;
        //lv_root_nodes.push_back(286462u); // (5.577 9.896 6.501) LV anterior base
       /* lv_root_nodes.push_back(237401u);  // (2.230 7.500 3.506) LV posterior 1
        lv_root_nodes.push_back(255997u); // (2.180 6.441 4.740) LV posterior 2
        lv_root_nodes.push_back(239979u);  // (6.333 5.949 3.740) LV mid septum (slightly anterior)
        std::vector<unsigned> rv_root_nodes;
       // rv_root_nodes.push_back(362760u); // (7.656 3.143 5.224) RV
        //rv_root_nodes.push_back(383391u); // (4.938 0.704 6.102) RV 2
        //rv_root_nodes.push_back(340698u); // (5.299 4.200 4.020) RV septum*/


	std::vector<unsigned> perm_vec = mesh.rGetNodePermutation();
	std::map<unsigned,unsigned> MapNewOldIndex;
	for (unsigned int i=0; i<perm_vec.size();i++)
	{			
		unsigned new_index = perm_vec[i];
		MapNewOldIndex[new_index] =i;
	}

        DijkstraStimuliCellFactory cell_factory(&gks_varier,
                                                mesh.rGetNodePermutation(),
                                                lv_root_nodes,
                                                rv_root_nodes,
                                                0.150,period,num_nodes,regions,MapNewOldIndex,edge_nodes_path,mindists_path); // 150 cm/s
        std::cout <<" i finish with cell fac" << std::endl;

        cell_factory.SetHeartGeometryInformation( &heart_geom_info );
        // Problem
        HeartEventHandler::EndEvent(HeartEventHandler::INITIALISE);
        // Problem immediately starts the EVERYTHING timer again
        HeartEventHandler::EndEvent(HeartEventHandler::EVERYTHING);
        BidomainProblem<3> bidomain_problem( &cell_factory );
        //bidomain_problem.PrintOutput(false); // Turn OFF output

        // Adaptive time controller
        PostUpstrokeTimeAdaptivityControllerCL time_controller(time_step, 0.25, bidomain_problem, period);
       	bidomain_problem.SetUseTimeAdaptivityController(true, &time_controller);

	//Endocardial Layer Disruption (modify RV endo conductivities by factor of 4 higher)
	
        std::cout <<" i enter the condctivities modifier" << std::endl;

	const std::vector<unsigned> endo =  heart_geom_info.rGetNodesOnRVSurface();


	
	
        bidomain_problem.SetMesh( &mesh );
	
        bidomain_problem.SetWriteInfo();
	
        bidomain_problem.SetUseHdf5DataWriterCache(true);
	bidomain_problem.SetHdf5DataWriterTargetChunkSizeAndAlignment(1048576);
        bidomain_problem.Initialise();


	// Get the tissue from the problem and apply the modifier
	ConductivitiesModifierEndo conduction_modifier(mesh,endo,perm_vec);

        std::cout <<" i finish with the condctivities modifier" << std::endl;

	BidomainTissue<3>* p_bidomain_tissue = bidomain_problem.GetBidomainTissue();
	std::cout <<" i get bidomain tissue" << std::endl;
	p_bidomain_tissue->SetConductivityModifier( &conduction_modifier );
	std::cout <<" i set conductivities via modifier" << std::endl;
	
        /* Ground the nodes of an element containing the (standard) right leg electrode */
        
/*for patient 02 i m unsure which triangular element is closest to surface electrode so i ground 2 of them*/
        /*grounded_nodes.push_back(mesh.rGetNodePermutation()[9463931u]); //in chaste numbering
        grounded_nodes.push_back(mesh.rGetNodePermutation()[9463825u]);
        grounded_nodes.push_back(mesh.rGetNodePermutation()[9463838u]);
        grounded_nodes.push_back(mesh.rGetNodePermutation()[9463838u]);
        grounded_nodes.push_back(mesh.rGetNodePermutation()[9463732u]);
        grounded_nodes.push_back(mesh.rGetNodePermutation()[9463748u]);*/
	//std::cout << "I select grounded nodes" << std::endl;
	
        //bidomain_problem.SetFixedExtracellularPotentialNodes(grounded_nodes);

 
        // Solve call -> the solver -> the assembler -> modified conductivities
	
	//delete p_bidomain_tissue; //free up memory
        std::cout <<"i start solving the problem" << std::endl;

	bidomain_problem.Solve();

	
//	mesh.SetMeshHasChangedSinceLoading(); //need thhis so that mesh is saved in permuted form and not copied!!
	//std::cout << "archiving mesh" <<std::endl;
	  /* Archive */
  /*      CardiacSimulationArchiver<BidomainProblem<3> >::Save(bidomain_problem, output_dir+"/saved_simulation");
    	std::cout << "loading saved sim" <<std::endl;
	 BidomainProblem<3>* p_bidomain_problem = CardiacSimulationArchiver<BidomainProblem<3> >::Load(output_dir+"/saved_simulation");
	 HeartConfig::Instance()->SetSimulationDuration(1000); //ms
         HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.02, 0.02, 1.0);
         HeartConfig::Instance()->SetOutputFilenamePrefix("results2");

	 BidomainTissue<3>* p_bidomain_tissue1 = p_bidomain_problem->GetBidomainTissue();
	 std::cout <<" i get bidomain tissue" << std::endl;
	 p_bidomain_tissue1->SetConductivityModifier( &conduction_modifier );
//	 std::cout <<" i set conductivities via modifier" << std::endl;
 //   	 conduction_modifier.~ConductivitiesModifierAnnulus3Dnew();

	 regions.clear();

         p_bidomain_problem->Solve();

         delete p_bidomain_problem;

 
	cell_factory.~DijkstraStimuliCellFactory();
	//conduction_modifier.~ConductivitiesModifierAnnulus3D();
*/	
	/*Pseudo ECG*/

/*	ChastePoint<3> point1(-2.3,-14.4,0.2); //V1
        PseudoEcgCalculator<3,3,1> ecg_calculator1 (mesh, point1, FileFinder(output_dir,RelativeTo::ChasteTestOutput), "results");
	ecg_calculator1.WritePseudoEcg();
	ecg_calculator1.~PseudoEcgCalculator();
	ChastePoint<3> point2(1.5,-14.5,-0.3); //V2
        PseudoEcgCalculator<3,3,1> ecg_calculator2 (mesh, point2, FileFinder(output_dir,RelativeTo::ChasteTestOutput), "results");
	ecg_calculator2.WritePseudoEcg();
	ecg_calculator2.~PseudoEcgCalculator();
	ChastePoint<3> point3(4.8,-15.1,-2.1); //V3
        PseudoEcgCalculator<3,3,1> ecg_calculator3 (mesh, point3, FileFinder(output_dir,RelativeTo::ChasteTestOutput), "results");
	ecg_calculator3.WritePseudoEcg();
	ecg_calculator3.~PseudoEcgCalculator();
	ChastePoint<3> point4(7.9,-15.1,-4.1); //V4
        PseudoEcgCalculator<3,3,1> ecg_calculator4 (mesh, point4, FileFinder(output_dir,RelativeTo::ChasteTestOutput), "results");
	ecg_calculator4.WritePseudoEcg();
	ecg_calculator4.~PseudoEcgCalculator();
	ChastePoint<3> point5(14.8,-10.3,-3.4); //V5
        PseudoEcgCalculator<3,3,1> ecg_calculator5 (mesh, point5, FileFinder(output_dir,RelativeTo::ChasteTestOutput), "results");
	ecg_calculator5.WritePseudoEcg();
	ecg_calculator5.~PseudoEcgCalculator();
	ChastePoint<3> point6(17.0,-2.7,-2.1); //V6
        PseudoEcgCalculator<3,3,1> ecg_calculator6 (mesh, point6, FileFinder(output_dir,RelativeTo::ChasteTestOutput), "results");
	ecg_calculator6.WritePseudoEcg();
	ecg_calculator6.~PseudoEcgCalculator();
	ChastePoint<3> point7(14.6,-6.0,15.8); //LA
        PseudoEcgCalculator<3,3,1> ecg_calculator7 (mesh, point7, FileFinder(output_dir,RelativeTo::ChasteTestOutput), "results");
	ecg_calculator7.WritePseudoEcg();
	ecg_calculator7.~PseudoEcgCalculator();
	ChastePoint<3> point8(-14.8,-5.6,-15.7); //RA
        PseudoEcgCalculator<3,3,1> ecg_calculator8 (mesh, point8, FileFinder(output_dir,RelativeTo::ChasteTestOutput), "results");
	ecg_calculator8.WritePseudoEcg();
	ecg_calculator8.~PseudoEcgCalculator();
	ChastePoint<3> point9(12.4,-8.2,-23.2); //LL
        PseudoEcgCalculator<3,3,1> ecg_calculator9 (mesh, point9, FileFinder(output_dir,RelativeTo::ChasteTestOutput), "results");
	ecg_calculator9.WritePseudoEcg();
	ecg_calculator9.~PseudoEcgCalculator();
	ChastePoint<3> point10(-11.7,-9.4,-22.2); //RL
        PseudoEcgCalculator<3,3,1> ecg_calculator10 (mesh, point10, FileFinder(output_dir,RelativeTo::ChasteTestOutput), "results");
	ecg_calculator10.WritePseudoEcg();
	ecg_calculator10.~PseudoEcgCalculator();
        * Archive */
       /* CardiacSimulationArchiver<BidomainProblem<3> >::Save(bidomain_problem, output_dir+"/saved_simulation");

        HeartEventHandler::Headings();
        HeartEventHandler::Report();
        HeartEventHandler::Headings();

        * Save permutation just in case */
       OutputFileHandler output_file_handler(output_dir, false); // collective
        if ( PetscTools::AmMaster() )
        {
            out_stream perm_out_stream = output_file_handler.OpenOutputFile("permutation.txt", std::ios::out);
            for (unsigned i=0; i<perm_vec.size(); i++)
            {
                (*perm_out_stream) << perm_vec[i] << "\n";
            }
        }



    }
};
