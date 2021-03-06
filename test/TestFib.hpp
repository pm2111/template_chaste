#include <cxxtest/TestSuite.h>
	#include "PetscSetupAndFinalize.hpp"
	#include "TetrahedralMesh.hpp"
	#include "StreeterFibreGenerator.hpp"
	#include "HeartEventHandler.hpp"
	
	class TestDoStreeterFibres : public CxxTest::TestSuite
	{
	public:
	        void TestMakeTheFibres() throw(Exception)
	    {
	        //const std::string path="/home/scratch/hdd/meshes/mikael2/tmpmesh_triangles_NoBath";
	        const std::string path="projects/pm2111/test/data/ARVC_002/HEART";
	        TrianglesMeshReader<3,3> mesh_reader(path);
	        TetrahedralMesh<3,3> mesh;
	        mesh.ConstructFromMeshReader(mesh_reader);
	       
	        const std::string epi_face_file = path+".epi";
	        const std::string rv_face_file = path+".rv";
	        const std::string lv_face_file = path+".lv";
		
	        StreeterFibreGenerator<3> fibre_generator(mesh);
	        
		fibre_generator.SetSurfaceFiles(epi_face_file, rv_face_file, lv_face_file, true);
	        //fibre_generator.SetApexToBase(2);

		c_vector<double, 3> axis;
                        axis[0] = -50.027;
                        axis[1] = -16.9187;
                        axis[2] = 39.3525;
                fibre_generator.SetApexToBase(axis);
	        //fibre_generator.GenerateOrthotropicFibreOrientation("TestDoStreeterFibres", "UPF_human_0.4mm_NoBath_apex.ortho", true);
        
OutputFileHandler handler("streeter_out", true);

 fibre_generator.SetLogInfo(false);
	        fibre_generator.WriteData(handler, "HEART.ortho");

	//        HeartEventHandler::Headings();
	//        HeartEventHandler::Report();

/* 
TrianglesMeshReader<3,3> TORSOmesh_reader( "projects/anamin/test/data/HCM_Meshes/TORSO" );
TetrahedralMesh<3,3> TORSOmesh;
TORSOmesh.ConstructFromMeshReader( TORSOmesh_reader );

//printf("\n\nen TORSOmesh hay    %d    elementos\n\n",TORSOmesh.GetNumElements() );


FILE *fid;
fid = fopen( "../testoutput/streeter_out/TORSO.ortho" , "a" );

for( unsigned e = mesh.GetNumElements() ; e < TORSOmesh.GetNumElements() ; e++ ){
  fprintf( fid , "1 0 0 0 1 0 0 0 1\n" );
}

fclose( fid );
*/


	    }
};
