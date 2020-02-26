#ifndef TESTDIJKSTRA_HPP_
        #define TESTDIJKSTRA_HPP_

        #include <cstdio>
        #include <fstream> // for ifstream
        #include <cxxtest/TestSuite.h>
        #include "PetscSetupAndFinalize.hpp"
        #include "PetscTools.hpp"
        #include <stack> // For stack
//      #include "boost.hpp"    
        /*#include "TetrahedralMesh.hpp"
        #include "HeartGeometryInformation.hpp"*/

        #include <boost/graph/use_mpi.hpp>
        #include <boost/mpi/communicator.hpp>
        #include <boost/graph/distributed/mpi_process_group.hpp>
        #include <boost/graph/graph_traits.hpp>
        #include <boost/graph/dijkstra_shortest_paths.hpp>
        #include <boost/graph/adjacency_list.hpp>
        #include <boost/graph/distributed/adjacency_list.hpp> // For distributedS
	
	class TestDijkstra : public CxxTest::TestSuite
	{
	public:
	    void noTestPrintSurfaceGraphs() throw(Exception)
	    {
	        /*EXIT_IF_PARALLEL
                 std::string patient_dir;
                  patient_dir = CommandLineArguments::Instance()->GetStringCorrespondingToOption("--patient_dir");
                  std::string save_dir;
                  save_dir = CommandLineArguments::Instance()->GetStringCorrespondingToOption("--save_dir");
	
	        std::string filepath = "projects/louiecn/test/data/nejib_20131010/FullMeshVolume_bin";
	
	        // Load mesh
	        TetrahedralMesh<3,3> mesh;
	        TrianglesMeshReader<3,3> mesh_reader(filepath);
	        mesh.ConstructFromMeshReader(mesh_reader);
	
	        // Construct heart geometry info for endocardial maps
	        const std::string epi_face_file = filepath+".epi";
	        const std::string rv_face_file = filepath+".rv";
	        const std::string lv_face_file = filepath+".lv";
	        HeartGeometryInformation<3> heart_geom_info(mesh, epi_face_file, lv_face_file, rv_face_file, true);
	
	        for (TetrahedralMesh<3,3>::NodeIterator iter=mesh.GetNodeIteratorBegin();
	             iter != mesh.GetNodeIteratorEnd();
	             ++iter)
	        {
	            unsigned node_index = iter->GetIndex();
	            const double rv_dist = heart_geom_info.rGetDistanceMapRightVentricle()[node_index];
	            const double lv_dist = heart_geom_info.rGetDistanceMapLeftVentricle()[node_index];
	            if(CompareDoubles::IsNearZero(rv_dist, 1e-6))
	            {
	                // This node is on an endo surface
	            }
	        }*/
	    }
	
	    void TestBoostDijkstra() throw(Exception)
	    {
	        /*boost::mpi::communicator comm(PETSC_COMM_WORLD, boost::mpi::comm_attach);
	        typedef boost::graph::distributed::mpi_process_group process_group_type;
	        process_group_type process_group(comm);*/
	
	        typedef boost::adjacency_list<boost::vecS, // previously listS
	                                      boost::vecS,
	                                      //boost::distributedS<process_group_type, boost::vecS>,
	                                      boost::undirectedS,
	                                      boost::no_property,                 // Vertex properties
	                                      boost::property<boost::edge_weight_t, double> // Edge properties
	                                      > graph_t;
	        graph_t g;
	        typedef boost::graph_traits < graph_t >::vertex_descriptor vertex_descriptor;
	        //typedef boost::graph_traits < graph_t >::edge_descriptor edge_descriptor;
	        boost::property_map<graph_t, boost::edge_weight_t>::type edge_distance = get(boost::edge_weight, g); // For alternative
	
	        typedef std::pair<unsigned, unsigned> Edge;
                 std::string patient_dir;
                  patient_dir = CommandLineArguments::Instance()->GetStringCorrespondingToOption("--patient_dir");
                  std::string save_dir;
                  save_dir = CommandLineArguments::Instance()->GetStringCorrespondingToOption("--save_dir");

	        /* Open edge file and read number of edges */
	        std::ifstream edge_file(patient_dir+"rv.edges", std::ios::in|std::ios::binary);
	        char header[36];
	        edge_file.getline(header,36);
	        unsigned num_edges = atoi(header);
	        std::cout<<"Number of edges="<<num_edges<<".\n";
	        std::ifstream weights_file(patient_dir+"rv.weights", std::ios::in|std::ios::binary);
	
	        /* Make arrays */
	        Edge* edge_array = new Edge[num_edges];
	        double* weights = new double[num_edges];
	        std::set<unsigned> nodes;
	
	        typedef boost::graph_traits<graph_t>::edge_descriptor edge_decriptor;
	        /* Read edges into array */
	        for ( unsigned i=0; i<num_edges; ++i )
	        {
	            unsigned edge_packet[2];
	            double weight_packet;
	            edge_file.read((char*)edge_packet, 2*sizeof(unsigned));
	            //if(PetscTools::AmMaster()){std::cout << "Edge has node indices " << edge_packet[0] << " and " << edge_packet[1] << ".\n";}
	            edge_array[i] = Edge(edge_packet[0], edge_packet[1]);
	            nodes.insert(edge_packet[0]);
	            nodes.insert(edge_packet[1]);
	            // Alternatively
	            weights_file.read((char*)&weight_packet, sizeof(double));
	            //if(PetscTools::AmMaster()){std::cout << "Edge has weight " << weight_packet << ".\n";}
	            // Add edge
	            edge_decriptor e = boost::add_edge(edge_packet[0], edge_packet[1], g).first;
	            edge_distance[e] = weight_packet;
	        }
	        edge_file.close();
	        weights_file.close();
	        unsigned num_nodes = nodes.size();
	        std::cout<<"Got "<<num_nodes<<" unique nodes.\n";
	        std::cout<<"Got "<<boost::num_vertices(g)<<" unique nodes.\n";
	
	/*
114	        unsigned num_nodes=5;
115	        //enum nodes { A, B, C, D, E };
116	        char name[] = "ABCDE";
117	
118	        //Edge edge_array[] = { Edge(A, C), Edge(B, B), Edge(B, D), Edge(B, E),
119	        //  Edge(C, B), Edge(C, D), Edge(D, E), Edge(E, A), Edge(E, B)
120	        //};
121	        Edge edge_array[] = { Edge(0, 2), Edge(1, 1), Edge(1, 3), Edge(1, 4),
122	          Edge(2, 1), Edge(2, 3), Edge(3, 4), Edge(4, 0), Edge(4, 1)
123	        };
124	        int weights[] = { 1, 2, 1, 2, 7, 3, 1, 1, 1 };
125	
126	        int num_arcs = sizeof(edge_array) / sizeof(Edge);*/
	
	        //graph_t g(edge_array, edge_array + num_edges, weights, num_nodes);
	
	        // Keeps track of the predecessor of each vertex
	        std::vector<vertex_descriptor> p(num_vertices(g));
	        // Keeps track of the distance to each vertex
	        std::vector<double> d(num_vertices(g));
	
	        FILE * min_dist_file = fopen((save_dir+"rv.mindists").c_str(), "w");
	//        fprintf(min_dist_file, "%u\n", num_nodes);
	
	        //vertex_descriptor s = vertex(A, g);
	        unsigned about_1_pc = num_nodes/100;
	        for(unsigned root=0; root<num_nodes; ++root)
	        {
	            if(root % about_1_pc == 0) std::cout << root/about_1_pc << "%" << std::endl;
	            vertex_descriptor s = vertex(root, g);
	            boost::dijkstra_shortest_paths
	              (g, s,
	               //predecessor_map(
	               //  make_iterator_property_map(p.begin(), get(boost::vertex_index, g))).
	               distance_map(
	                 make_iterator_property_map(d.begin(), get(boost::vertex_index, g)))
	               );
	
	         //   fwrite(&d[0], sizeof(double), d.size() , min_dist_file);
			fwrite(&d[0], sizeof(double), root , min_dist_file);
	        }
	
	        /*std::cout << "distances and parents:" << std::endl;
156	        boost::graph_traits < graph_t >::vertex_iterator vi, vend;
157	        for (boost::tie(vi, vend) = vertices(g); vi != vend; ++vi)
158	        {
159	            std::cout << d[*vi] << std::endl;
160	            //std::cout << "distance(" << name[*vi] << ") = " << d[*vi] << ", ";
161	            //std::cout << "parent(" << name[*vi] << ") = " << name[p[*vi]] << std::endl;
162	        }*/
	
	        delete[] edge_array;
	        delete[] weights;
	    }
	};
	
#endif /*TESTDIJKSTRA_HPP_*/
