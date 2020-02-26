#include <cxxtest/TestSuite.h>
#include "PetscSetupAndFinalize.hpp"
#include "MonodomainProblem.hpp"
#include "HeartGeometryInformation.hpp"
#include "CardiacSimulationArchiver.hpp"
#include "HeartEventHandler.hpp"
#include "SingleTraceOutputModifier.hpp"
#include "RegularStimulus.hpp"
#include "ORdGksVarier.hpp"
#include "PostUpstrokeTimeAdaptivityControllerCL_mono.hpp"
#include "TypeAsigment_stimulus.hpp"
#include "PseudoEcgCalculator.hpp"

class TestSolveTorso : public CxxTest::TestSuite
{

public:
    void TestSolve() throw(Exception)
    {
        HeartEventHandler::BeginEvent(HeartEventHandler::EVERYTHING);

        /* Mesh to use */
        std::string filepath = "projects/anamin/test/data/HCM_Meshes/HEART";
        HeartConfig::Instance()->SetMeshFileName(filepath, cp::media_type::Orthotropic);

        // Add this to the IN_MESH time
        HeartEventHandler::BeginEvent(HeartEventHandler::READ_MESH);
        DistributedTetrahedralMesh<3,3> mesh;
        TrianglesMeshReader<3,3> mesh_reader(filepath);
        mesh.ConstructFromMeshReader(mesh_reader);
        HeartEventHandler::EndEvent(HeartEventHandler::READ_MESH);
        HeartEventHandler::BeginEvent(HeartEventHandler::INITIALISE);

        /* Construct HeartGeometryInformation */
        const std::string epi_face_file = filepath+".epi";
        const std::string rv_face_file = filepath+".rv";
        const std::string lv_face_file = filepath+".lv";
        HeartGeometryInformation<3> heart_geom_info(mesh, epi_face_file, lv_face_file, rv_face_file, true);

        /* Simulation options */
        HeartConfig::Instance()->SetSimulationDuration(1600.0); // ms
        std::string output_dir = "TestECGAlfonso";
        HeartConfig::Instance()->SetOutputDirectory(output_dir);
        HeartConfig::Instance()->SetOutputFilenamePrefix("results");
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.02, 0.02, 1.0);
        //HeartConfig::Instance()->SetKSPPreconditioner("blockdiagonal");
	//HeartConfig::Instance()->SetKSPPreconditioner("ldufactorisation");
	HeartConfig::Instance()->SetKSPSolver("symmlq");

        /* Conductivities */
        // Intra conductivities result in close to 67 30 17 cm/s conduction velocities as in pig
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(1.5, 0.45, 0.225));

        /* Visualisers */
        HeartConfig::Instance()->SetVisualizeWithMeshalyzer(false);
        HeartConfig::Instance()->SetVisualizeWithCmgui(false);
        HeartConfig::Instance()->SetVisualizeWithVtk(false);
        HeartConfig::Instance()->SetVisualizeWithParallelVtk(true);
        HeartConfig::Instance()->SetUseStateVariableInterpolation(false);
    

        /* Output variables */
        /*std::vector<std::string> output_variables;
        output_variables.push_back("membrane_slow_delayed_rectifier_potassium_current_conductance");
        HeartConfig::Instance()->SetOutputVariables( output_variables );*/

        /* APD map */
        std::vector<std::pair<double,double> > apd_map;
        apd_map.push_back(std::pair<double, double>(90.0, 0.0)); // APD90, 0 mV threshold
        HeartConfig::Instance()->SetApdMaps(apd_map);

        /* Activation time map */
        std::vector<double> upstroke_time_map;
        upstroke_time_map.push_back(0.0); // 0 mV threshold
        HeartConfig::Instance()->SetUpstrokeTimeMaps(upstroke_time_map);

        // Gks varier (apicobasal, transmural, interventricular)
        ORdGksVarier gks_varier(&heart_geom_info, true, false, false);
        // LV and RV activation root nodes
	std::vector<unsigned> lv_root_nodes;
        lv_root_nodes.push_back(299671u); // (5.577 9.896 6.501) LV anterior base
        lv_root_nodes.push_back(240981u);  // (2.230 7.500 3.506) LV posterior 1
        lv_root_nodes.push_back(267696u); // (2.180 6.441 4.740) LV posterior 2
        lv_root_nodes.push_back(250820u);  // (6.333 5.949 3.740) LV mid septum (slightly anterior)
        std::vector<unsigned> rv_root_nodes;

        rv_root_nodes.push_back(380419u); // (7.656 3.143 5.224) RV
        rv_root_nodes.push_back(404783u); // (4.938 0.704 6.102) RV 2
        rv_root_nodes.push_back(356982u); // (5.299 4.200 4.020) RV septum

	TypeAsigment_stimulus cell_factory(&gks_varier,
                                                mesh.rGetNodePermutation(),
                                                lv_root_nodes,
                                                rv_root_nodes,
                                                0.179); // 150 cm/s

        cell_factory.SetHeartGeometryInformation( &heart_geom_info );
        // Problem
        HeartEventHandler::EndEvent(HeartEventHandler::INITIALISE);
        // Problem immediately starts the EVERYTHING timer again
        HeartEventHandler::EndEvent(HeartEventHandler::EVERYTHING);
        MonodomainProblem<3> monodomain_problem( &cell_factory );
        //bidomain_problem.PrintOutput(false); // Turn OFF output

        // Adaptive time controller
        PostUpstrokeTimeAdaptivityControllerCL_mono time_controller(0.02, 1.0, monodomain_problem, 800);
        monodomain_problem.SetUseTimeAdaptivityController(true, &time_controller);


        monodomain_problem.SetMesh( &mesh );
        monodomain_problem.SetWriteInfo();
        monodomain_problem.SetUseHdf5DataWriterCache(true);
        monodomain_problem.Initialise();


        // Solve call -> the solver -> the assembler -> modified conductivities
        monodomain_problem.Solve();

        /* Archive */
        // CardiacSimulationArchiver<BidomainWithBathProblem<3> >::Save(bidomain_problem, output_dir+"/saved_simulation");

        HeartEventHandler::Headings();
        HeartEventHandler::Report();
        HeartEventHandler::Headings();

        /* Save permutation just in case */
        OutputFileHandler output_file_handler(output_dir, false); // collective
        if ( PetscTools::AmMaster() )
        {
            out_stream perm_out_stream = output_file_handler.OpenOutputFile("permutation.txt", std::ios::out);
            std::vector<unsigned> perm_vec = mesh.rGetNodePermutation();
            for (unsigned i=0; i<perm_vec.size(); ++i)
            {
                (*perm_out_stream) << perm_vec[i] << "\n";
            }
        }
    }
};
