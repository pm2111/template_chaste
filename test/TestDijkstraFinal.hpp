		#ifndef TESTDIJKSTRA_HPP_
        #define TESTDIJKSTRA_HPP_

        #include <cstdio>
        #include <fstream> // for ifstream
        #include <cxxtest/TestSuite.h>
        #include "PetscSetupAndFinalize.hpp"
        #include "PetscTools.hpp"
        #include <stack> // For stack
        #include "mpi.h"
        /*#include "TetrahedralMesh.hpp"
        #include "HeartGeometryInformation.hpp"*/

        //#include <boost/graph/use_mpi.hpp>
        //#include <boost/mpi/communicator.hpp>
        //#include <boost/graph/distributed/mpi_process_group.hpp>
        #include <boost/graph/graph_traits.hpp>
        #include <boost/graph/dijkstra_shortest_paths.hpp>
        #include <boost/graph/adjacency_list.hpp>
        //#include <boost/graph/distributed/adjacency_list.hpp> // For distributedS

        class TestDijkstraFinal : public CxxTest::TestSuite
        {
        public:
            void Final() throw(Exception)
                {

		
		double *buffer;
                buffer= (double*) malloc (sizeof(double)*num_nodes);
                size_t result;
                FILE * min_dist_file1 = fopen((save_dir+"lv.mindists").c_str(), "w");
                 int numprocs2=48*16;

                FILE* pFile;
                for (int i=0; i<numprocs2;i++)
                {	 
			 std::cout <<i <<std::endl;
                         pFile = fopen ( ("lv.mindists"+std::to_string(i)).c_str() , "rb" );

                         result=fread(buffer,sizeof(double),num_nodes,pFile);
                         fwrite(&result, sizeof(double), num_nodes , min_dist_file1);
			 fclose(pFile);

                }
		}
};

#endif 
