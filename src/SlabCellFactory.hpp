#ifndef APEXSTIMULUSCELLFACTORYORdPat05_HPP_
#define APEXSTIMULUSCELLFACTORYORdPat05_HPP_


#include <algorithm>
#include <fstream>
#include <sstream>


#include <ChastePoint.hpp>

#include "ToRORd_fkatp_endo_mid_ACMCvodeOpt.hpp"
#include <algorithm>

class DijkstraStimuliCellFactory : public AbstractCardiacCellFactory<3>
{
private:
    ORdGksVarierDTI024* mpORdGksVarierDTI024;
    const std::vector<unsigned>& mrNodePerm;
    std::vector<unsigned> mLvRootNodeNewIndices; // Indices in distance file of root points
    std::vector<unsigned> mRvRootNodeNewIndices;
    //const std::vector<unsigned> rLvRootNodes;
   // const std::vector<unsigned> rRvRootNodes;
    std::vector<boost::shared_ptr<SimpleStimulus> > mStimuli;
    std::map<unsigned,unsigned> mMappingToLvDists; // Mapping between new index and position in distance file
    std::map<unsigned,unsigned> mMappingToRvDists;
    FILE * mpLvDistancesFile; // File that contains the minimum distances
    FILE * mpRvDistancesFile;
    double mSpeed; // Wave "speed"

    // Function to ease access to off-diagonal triangle binary file
    unsigned RowColToIndexMap(unsigned i, unsigned j)
    {
        assert(i!=j); // Diagonal is zero distance and not saved in file
        if(j>i)
        {
            std::swap(i,j);
        }
        // i<j
        return (i*(i-1))/2 + j; // Sum of 1 to i-1, plus j
    }

    void LoadFiles(std::string nodesFile,
                   std::map<unsigned,unsigned>& rMappingToXvDists,
                   const std::vector<unsigned>& rXvRootNodes,
                   std::vector<unsigned>& rXvRootNodeNewIndices)
    {
	        std::cout << "I m inside loader " << std::endl ;

        // Open Dijkstra files
        // File that gives mapping between original node indices and indices in distances file
        std::ifstream nodes_file(nodesFile.c_str());
               std::cout << "I open lv_edge_nodes file " << std::endl ;

        std::string buffer;
        unsigned i=0;
        while ( nodes_file >> buffer )
        {
            unsigned old_index = std::atoi(buffer.c_str());
//	                std::cout << "old index is  " << old_index << "size of map" <<mrNodePerm.size() << std::endl ;
            // For this old_index, find the new index, and set rMappingToXvDists s.t.
            // rMappingToXvDists[new_index] gives the index for looking up in distance file
//          std::cout << mrNodePerm.size() <<std::endl;  
	   unsigned new_index = mrNodePerm[i];
//	        std::cout <<  "new index  "<< new_index << "old index" << old_index << "size of map" <<rMappingToXvDists.size() << std::endl ;

            //rMappingToXvDists[i] = i;

            rMappingToXvDists[new_index] = i;
  //        std::cout << "get map index  "<< i << std::endl ;

            // Also see if this old_index is one of our root points. If it is, do the same as above.
            for(unsigned r=0; r<rXvRootNodes.size(); ++r)
            {
                if(old_index == rXvRootNodes[r])
                {
                    rXvRootNodeNewIndices.push_back(i);
                    break;
                }
            }
            i++;
        }
        if( !nodes_file.eof() )
        {
            EXCEPTION("Couldn't read lv_edges_nodes file.");
        }
        assert( ! rXvRootNodeNewIndices.empty() ); // Need at least 1 root node.
        nodes_file.close();

    }

    double GetMinimumDistanceFromRootsForNode(unsigned nodeIndex,
                                              std::map<unsigned,unsigned>& mMappingToXvDists,
                                              std::vector<unsigned>& rXvRootNodeNewIndices,
                                              FILE* pXvDistancesFile)
    {
        unsigned new_index;
        try
        {
            new_index = mMappingToXvDists.at(nodeIndex);
        }
        catch (const std::out_of_range& oor)
        {
            return -1.0; // Nonsense
        }
        // Loop over root nodes and look up distance. Keep the smallest.
        double buffer, min_distance=DBL_MAX;
        for(unsigned i=0; i<rXvRootNodeNewIndices.size(); ++i)
        {
            if( new_index==rXvRootNodeNewIndices[i] ) // This is the root node
            {
                return 0.0;
            }
            unsigned tri_index = RowColToIndexMap(new_index, rXvRootNodeNewIndices[i]);
            // Read distance from file
            fseek(pXvDistancesFile, tri_index*sizeof(double), SEEK_SET);
            fread(&buffer, sizeof(double), 1, pXvDistancesFile); // Read a double
            min_distance = buffer<min_distance ? buffer : min_distance;
        }
        return min_distance;
    }

public:
    DijkstraStimuliCellFactory(const std::vector<unsigned>& rNodePerm)
    : AbstractCardiacCellFactory<3>(),
    mrNodePerm(rNodePerm)
	{
	        std::cout << "i get to cell fac" << std::endl;

		
	}
/*
    ~DijkstraStimuliCellFactory()
	{
		        std::vector<int> mMappingToLvDists;
	}
*/

    
        // The following has two functions. First, load the list of LV/RV nodes (1) and use to populate
        // the map so it can be used thusly: posn_in_dist_file = mMappingToXvDists.at(new_node_index)
        // Then does the same for each of the specified root nodes (3) into member var (4).
                
  //      std::cout << "I start opening lv_edge_nodes " << std::endl ;

	//LoadFiles("/data/blanca-rodriguez/shug5389/meshes/chaste_tiny/lv_edge_nodes", mMappingToLvDists, rLvRootNodes, mLvRootNodeNewIndices );
        //LoadFiles("/data/blanca-rodriguez/shug5389/meshes/chaste_tiny/rv_edge_nodes", mMappingToRvDists, rRvRootNodes, mRvRootNodeNewIndices );
    

   
    //    fclose(mpLvDistancesFile);
      //  fclose(mpRvDistancesFile);
    

    AbstractCvodeCell* CreateCardiacCellForTissueNode(Node<3>* pNode)// method creates a cell for a tissue simulation
    {
//	std::cout << "I get to line 134 of Dijkstra Cell Factory " << std::endl ;
        AbstractCvodeCell* p_cell = new CellToRORd_fkatp_endo_mid_ACMFromCellMLCvodeOpt(mpSolver, mpZeroStimulus); //defines a cell model and that the cell will not be stimulated 

        unsigned node_index = pNode->GetIndex(); // New index
/*	//std::cout << "I get to line 140 of cell factory file" << std::endl;
        HeartGeometryInformation<3>* p_heart_geom_info = GetHeartGeometryInformation();
        const double lv_dist = p_heart_geom_info->rGetDistanceMapLeftVentricle()[node_index]; //get distance to the LV
        const double rv_dist = p_heart_geom_info->rGetDistanceMapRightVentricle()[node_index]; //get distance to the RV

        double min_distance = -1.0; // Nonsense

	std::ostringstream strs;
	strs << node_index;
	std::string str = strs.str();
*/
	//std::cout << "I start to compute the min distances for node " << str << std::endl;
          ChastePoint<3> location;
        location = pNode->rGetLocation();
//	std::cout << location[0] << " " << location[1] << location[2] <<std::endl;
          if (location[0]<0.08 && location[1]<0.08 && location[2]<0.08)
        {

                    boost::shared_ptr<SimpleStimulus> p_stim(new SimpleStimulus(-100000.0, 1.0, 1.0));

		    p_cell->SetStimulusFunction(p_stim);//only pick the stimulus time as the time corresponding to the time it take for the signal to travel to it!!
	}
	            //p_cell->SetStimulusFunction(p_stim);//only pick the stimulus time as the time corresponding to the time it take for the signal to travel to it!!

        //std::cout << "I start to compute the min distances for node " << str << std::endl;
        //
        //        /*Cell Assignment*/
        //                //double distance_epi = p_heart_geom_info->rGetDistanceMapEpicardium()[node_index];
        //                        //double distance_endo = std::min(lv_dist, rv_dist);
        //                                //double relative_position = distance_endo / (distance_endo + distance_epi);
        //
        //
        //                                        if (location[0]<0.2 && location[1]<0.2 && location[2]<0.2)
        //                                                {
        //
        //                                                                        p_cell = new CellORd2011epi_fkatpFromCellMLCvodeOpt(mpSolver, mpZeroStimulus);
        //
        //                                                                                }
        //
	/*  if (CompareDoubles::IsNearZero(lv_dist, 1e-3))
        {
            min_distance = GetMinimumDistanceFromRootsForNode(node_index, mMappingToLvDists, mLvRootNodeNewIndices, mpLvDistancesFile); //compute the min distance to the closest root node in LV
        }
        else if (CompareDoubles::IsNearZero(rv_dist, 1e-3))
        {
            min_distance = GetMinimumDistanceFromRootsForNode(node_index, mMappingToRvDists, mRvRootNodeNewIndices, mpRvDistancesFile);//compute the min distance to the closest root node in RV
        }
         if ( min_distance >= 0.0 )   
        {
            // Convert this into a time
           unsigned upstroke_time = min_distance/mSpeed + 0.5; // Result is rounded! In ms.
            // Have we got a stimulus for this time? If not, extend.
           if( upstroke_time >= mStimuli.size() )
            {
                for(unsigned t=mStimuli.size(); t<=upstroke_time; ++t)
                {
                    boost::shared_ptr<SimpleStimulus> p_stim(new SimpleStimulus(-70000.0, 1.0, (double)t));
                   // boost::shared_ptr<SimpleStimulus> p_stim(new SimpleStimulus(-70000.0, 1.0, 0.5)); //stimulate at 0.5ms
                    mStimuli.push_back(p_stim);
               }
            }
            p_cell->SetStimulusFunction(mStimuli[upstroke_time]);//only pick the stimulus time as the time corresponding to the time it take for the signal to travel to it!!
        }*/
//	std::cout << "I have assigned a stimulus for node " << str << std::endl;
        //const double gks = p_cell->GetParameter("membrane_slow_delayed_rectifier_potassium_current_conductance");
        //const double gks_scaling = mpORdGksVarierDTI024->GetScalingFactor1(*pNode);
       // p_cell->SetParameter("membrane_slow_delayed_rectifier_potassium_current_conductance", gks_scaling*gks); //apply a scaling of Iks across the apexbase of the ventricle!

        //p_cell->SetTolerances(1e-4,1e-6);

        return p_cell;
    }
};

#endif//DIJKSTRASTIMULICELLFACTORYORdPat02_HPP_
