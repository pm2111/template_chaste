#include <cxxtest/TestSuite.h>
#include "PetscSetupAndFinalize.hpp"
#include "PetscTools.hpp"
#include "TetrahedralMesh.hpp"
#include "StreeterFibreGenerator.hpp"
#include "HeartEventHandler.hpp"
#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <functional>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_expression.hpp>
#include <boost/numeric/ublas/traits/c_array.hpp>
#include <boost/numeric/ublas/assignment.hpp>
#include <boost/assign.hpp>
#include <boost/assign/std/vector.hpp> // for 'operator+=()'
#include <boost/assert.hpp> 
#include <HeartGeometryInformation.hpp> 
using namespace boost::assign; // bring 'operator+=()' into scope

class TestDoStreeterFibres : public CxxTest::TestSuite
{
public:
	void TestMakeTheFibres() throw(Exception)
    {

        std::string output_dir;
        output_dir = CommandLineArguments::Instance()->GetStringCorrespondingToOption("--save_dir");
        std::string apexbase_path;
        apexbase_path = CommandLineArguments::Instance()->GetStringCorrespondingToOption("--apexbase");
        char* dummy;
        dummy = CommandLineArguments::Instance()->GetValueCorrespondingToOption("--dists");
	unsigned dist=std::atoi(dummy);

        //const std::string path="/home/scratch/hdd/meshes/mikael2/tmpmesh_triangles_NoBath";
	const std::string path=output_dir+"HEART";
        TrianglesMeshReader<3,3> mesh_reader(path);
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        const std::string epi_face_file = path+".epi";
        const std::string rv_face_file = path+".rv";
        const std::string lv_face_file = path+".lv";
	std::vector<double> dist_epi;
        std::vector<double> dist_endo;
        std::vector<double> dist_RV;
        std::vector<double> dist_LV;

	if (dist==1)
{
	
        HeartGeometryInformation<3> heart_geom_info(mesh, epi_face_file, lv_face_file, rv_face_file, true);
	dist_epi=heart_geom_info.rGetDistanceMapEpicardium();
	dist_endo=heart_geom_info.rGetDistanceMapEndocardium();
        dist_LV=heart_geom_info.rGetDistanceMapLeftVentricle();
        dist_RV=heart_geom_info.rGetDistanceMapRightVentricle();
        OutputFileHandler output_file_handler(output_dir, false); // collective

	out_stream perm_out_stream = output_file_handler.OpenOutputFile("dists.txt", std::ios::out);

	for (int i=0; i<dist_epi.size();i++)
	{	
	            (*perm_out_stream) << dist_epi[i] << " " << dist_endo[i] << " " << dist_LV[i] << " "<< dist_RV[i] << "\n";


	}
}

        StreeterFibreGenerator<3> fibre_generator(mesh);
        fibre_generator.SetSurfaceFiles(epi_face_file, rv_face_file, lv_face_file, true);

	std::vector<double> rainfall;    // a vector to hold apexbase data
	// open file    
	std::ifstream inputFile(apexbase_path);

	// test file open   
	if (inputFile) {       
	double value;
	
	// read the elements in the file into a vector  
	while ( inputFile >> value ) {
	rainfall.push_back(value);
	}
	}
	
	// close the file
	
	/*for(int i=0;i<6;i++){
	std::cout<<rainfall[i]<<std::endl;
	}*/
	/*const std::vector<double,3> apexbase;
	apexbase.push_back(rainfall[3]-rainfall[0]);
	apexbase.push_back(rainfall[4]-rainfall[1]);
	apexbase.push_back(rainfall[5]-rainfall[2]);*/
	//const std::c_vector<double,3> apexbase;
	
	double magnitude = std::sqrt(	(rainfall[3]-rainfall[0]) *(rainfall[3]-rainfall[0]) + (rainfall[4]-rainfall[1])*(rainfall[4]-rainfall[1]) + (rainfall[5]-rainfall[2])*(rainfall[5]-rainfall[2]));
	
	boost::numeric::ublas::c_vector<double,3> apexbase;
	

	apexbase <<= (rainfall[3]-rainfall[0])/magnitude,(rainfall[4]-rainfall[1])/magnitude,(rainfall[5]-rainfall[2])/magnitude ;
	
	
	// = { rainfall[3]-rainfall[0], rainfall[4]-rainfall[1], rainfall[5]-rainfall[2]} ;
        fibre_generator.SetApexToBase(apexbase);
//        fibre_generator.GenerateOrthotropicFibreOrientation("TestDoStreeterFibres", "UPF_human_0.4mm_NoBath_apex.ortho", true);
        OutputFileHandler handler("streeter", false);
        fibre_generator.WriteData(handler, "HEART.ortho");

        HeartEventHandler::Headings();
//        HeartEventHandler::Report();
    }
};
