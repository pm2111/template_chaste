#ifndef TypeAsigment_stimulus_HPP_
#define TypeAsigment_stimulus_HPP_

#include <algorithm>
#include <fstream>

#include "ORd2011epi_fkatpCvodeOpt.hpp"
#include "ORd2011m_fkatpCvodeOpt.hpp"
#include "ORd2011endo_fkatpCvodeOpt.hpp"

class TypeAsigment_stimulus : public AbstractCardiacCellFactory<3>
{
private:
    ORdGksVarier* mpORdGksVarier;
    const std::vector<unsigned>& mrNodePerm;
    std::vector<unsigned> mLvRootNodeNewIndices; // Indices in distance file of root points
    std::vector<unsigned> mRvRootNodeNewIndices;
    std::vector<boost::shared_ptr<RegularStimulus> > mStimuli;
    std::map<unsigned,unsigned> mMappingToLvDists; // Mapping between new index and position in distance file
    std::map<unsigned,unsigned> mMappingToRvDists;
    FILE * mpLvDistancesFile; // File that contains the minimum distances
    FILE * mpRvDistancesFile;
    double mSpeed; // Wave "speed"

    // Function to ease access to off-diagonal triangle binary file
    unsigned long long RowColToIndexMap(unsigned long long i, unsigned long long j)
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
        // Open Dijkstra files
        // File that gives mapping between original node indices and indices in distances file
        std::ifstream nodes_file(nodesFile.c_str());
        std::string buffer;
        unsigned i=0;
        while ( nodes_file >> buffer )
        {
            unsigned old_index = std::atoi(buffer.c_str());
            // For this old_index, find the new index, and set rMappingToXvDists s.t.
            // rMappingToXvDists[new_index] gives the index for looking up in distance file
            unsigned new_index = mrNodePerm[old_index];
            rMappingToXvDists[new_index] = i;

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
            unsigned long long tri_index = RowColToIndexMap(new_index, rXvRootNodeNewIndices[i]);
            // Read distance from file
            fseek(pXvDistancesFile, tri_index*sizeof(double), SEEK_SET);
            fread(&buffer, sizeof(double), 1, pXvDistancesFile); // Read a double
            min_distance = buffer<min_distance ? buffer : min_distance;
        }
        return min_distance;
    }

public:
    TypeAsigment_stimulus(ORdGksVarier* pGksVarier,
                              const std::vector<unsigned>& rNodePerm,
                              const std::vector<unsigned>& rLvRootNodes,
                              const std::vector<unsigned>& rRvRootNodes,
                              double speed)
    : AbstractCardiacCellFactory<3>(),
    mpORdGksVarier(pGksVarier),
    mrNodePerm(rNodePerm),
    mpLvDistancesFile(fopen("projects/anamin/test/data/HCM_Meshes/lv.mindists","rb")),
    mpRvDistancesFile(fopen("projects/anamin/test/data/HCM_Meshes/rv.mindists","rb")),
    mSpeed(speed)
    {
        // The following has two functions. First, load the list of LV/RV nodes (1) and use to populate
        // the map so it can be used thusly: posn_in_dist_file = mMappingToXvDists.at(new_node_index)
        // Then does the same for each of the specified root nodes (3) into member var (4).
        LoadFiles("projects/anamin/test/data/HCM_Meshes/lv_edge_nodes", mMappingToLvDists, rLvRootNodes, mLvRootNodeNewIndices );
        LoadFiles("projects/anamin/test/data/HCM_Meshes/rv_edge_nodes", mMappingToRvDists, rRvRootNodes, mRvRootNodeNewIndices );
    }

    ~TypeAsigment_stimulus()
    {
        fclose(mpLvDistancesFile);
        fclose(mpRvDistancesFile);
    }

    AbstractCvodeCell* CreateCardiacCellForTissueNode(Node<3>* pNode)
    {
        AbstractCvodeCell* p_cell;

        unsigned node_index = pNode->GetIndex(); // New index

        HeartGeometryInformation<3>* p_heart_geom_info = GetHeartGeometryInformation();
        const double lv_dist = p_heart_geom_info->rGetDistanceMapLeftVentricle()[node_index];
        const double rv_dist = p_heart_geom_info->rGetDistanceMapRightVentricle()[node_index];

        double min_distance = -1.0; // Nonsense
        if (CompareDoubles::IsNearZero(lv_dist, 1e-3))
        {
            min_distance = GetMinimumDistanceFromRootsForNode(node_index, mMappingToLvDists, mLvRootNodeNewIndices, mpLvDistancesFile);
        }
        else if (CompareDoubles::IsNearZero(rv_dist, 1e-3))
        {
            min_distance = GetMinimumDistanceFromRootsForNode(node_index, mMappingToRvDists, mRvRootNodeNewIndices, mpRvDistancesFile);
        }
	// distances
        double distance_epi = p_heart_geom_info->rGetDistanceMapEpicardium()[node_index];
        double distance_endo = std::min(lv_dist, rv_dist);
        double relative_position = distance_endo / (distance_endo + distance_epi);

	if (relative_position > 0.75)
	{
                p_cell = new CellORd2011epi_fkatpFromCellMLCvodeOpt(mpSolver, mpZeroStimulus);
        }
	else if(relative_position>0.45)
	{
		p_cell = new CellORd2011m_fkatpFromCellMLCvodeOpt(mpSolver, mpZeroStimulus);
	}
	else
	{
                p_cell = new CellORd2011endo_fkatpFromCellMLCvodeOpt(mpSolver, mpZeroStimulus);
	}


	if ( min_distance >= 0.0)
        {
            // Convert this into a time
            unsigned upstroke_time = min_distance/mSpeed + 0.5; // Result is rounded! In ms.
            // Have we got a stimulus for this time? If not, extend.
            if( upstroke_time >= mStimuli.size() )
            {
                for(unsigned t=mStimuli.size(); t<=upstroke_time; ++t)
                {
                    boost::shared_ptr<RegularStimulus> p_stim(new RegularStimulus(-120000.0, 1.0, 800, (double)t));
               //boost::shared_ptr<SimpleStimulus_Vcheck2> p_stim(new SimpleStimulus_Vcheck2(-120000.0, 1.0, (double)t, p_cell, -20));
 
		   mStimuli.push_back(p_stim);
                }
            }
            p_cell->SetStimulusFunction(mStimuli[upstroke_time]);
        }

       const double gks = p_cell->GetParameter("membrane_slow_delayed_rectifier_potassium_current_conductance");
       const double gks_scaling = mpORdGksVarier->GetScalingFactor(*pNode);
       p_cell->SetParameter("membrane_slow_delayed_rectifier_potassium_current_conductance", gks_scaling*gks);

        p_cell->SetTolerances(1e-4,1e-6);

        return p_cell;
    }
};

#endif //TypeAsigment_stimulusHPP_
