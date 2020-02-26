#ifndef DIJKSTRASTIMULUSCELLFACTORYORDEPIPAT05APICOBASALGRADTHREECELLTYPES_HPP_
#define DIJKSTRASTIMULUSCELLFACTORYORDEPIPAT05APICOBASALGRADTHREECELLTYPES_HPP_

#include <algorithm>
#include <fstream>
#include <math.h>       /* remainder */


#include "ORd2011epi_fkatpCvodeOpt.hpp"
#include "ORd2011m_fkatpCvodeOpt.hpp"
#include "ORd2011endo_fkatpCvodeOpt.hpp"

class DijkstraStimuliCellFactory : public AbstractCardiacCellFactory<3>
{
private:
    ORdGksVarierDTI024* mpORdGksVarierDTI024;
    const std::vector<unsigned>& mrNodePerm;
    std::vector<unsigned> mLvRootNodeNewIndices; // Indices in distance file of root points
    std::vector<unsigned> mRvRootNodeNewIndices;
    const unsigned mnum_nodes;
    std::vector<unsigned> heterogeneities;
    //const std::vector<unsigned> rLvRootNodes;
    //const std::vector<unsigned> rRvRootNodes;
    std::vector<boost::shared_ptr<RegularStimulus> > mStimuli;
    std::map<unsigned,unsigned> mMappingToLvDists; // Mapping between new index and position in distance file
    std::map<unsigned,unsigned> mMappingToRvDists;
    std::map<unsigned,unsigned> MappingOldIndexNewIndex;
    FILE * mpLvDistancesFile; // File that contains the minimum distances
    FILE * mpRvDistancesFile;
    double mSpeed; // Wave "speed"
    unsigned dummy; //for monitoring of cell production process
    double mperiod;
    const std::vector<std::vector<double> >& mregions;
    const std::map<unsigned,unsigned>& mMapNewOldIndex;

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
 void LoadHeterogeneities(std::string heterogen_file,
                    std::vector<unsigned>& rheterogeneities,std::map<unsigned,unsigned>& rMappingOldIndexNewIndex)
    {
        // Open Dijkstra files
        // File that gives mapping between original node indices and indices in distances file
	//std::cout << "i enter the loader function" << std::endl;
        std::ifstream file(heterogen_file.c_str());
	//std::cout << "i open the fibrosis file" << std::endl;
        std::string buffer1;

        unsigned i=0;
        while ( file >> buffer1 )
        {
            unsigned heterogen_value = std::atoi(buffer1.c_str());

	    std::cout << i << endl;
            // For this old_index, find the new index corresponding to the parallel computation index
            unsigned new_index = mrNodePerm[i];

            rheterogeneities.push_back(heterogen_value);
	    rMappingOldIndexNewIndex[new_index] = i; //this way i have a value (i) corresponding to the non partitioned node ordering

            i++;
        }
        if( !file.eof() )
        {
            EXCEPTION("Couldn't read heterogeneities file.");
        }
        file.close();

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
    DijkstraStimuliCellFactory(ORdGksVarierDTI024* pGksVarier,
                              const std::vector<unsigned>& rNodePerm,
                              const std::vector<unsigned>& rLvRootNodes,
                              const std::vector<unsigned>& rRvRootNodes,
                              double speed,double period,unsigned num_nodes, const std::vector<std::vector<double> >& regions,const std::map<unsigned,unsigned>& MapNewOldIndex)
    : AbstractCardiacCellFactory<3>(),
     mpORdGksVarierDTI024(pGksVarier),
    mrNodePerm(rNodePerm),
    mnum_nodes(num_nodes),
    mperiod(period),
    mregions(regions),
    mMapNewOldIndex(MapNewOldIndex),
    mpLvDistancesFile(fopen("/gpfs/scratch/pr1efz00/pr1efz06/testoutput/lv.mindists","rb")),
    mpRvDistancesFile(fopen("/gpfs/projects/pr1efz00/pr1efz06/testoutput/rv.mindists","rb")),
    mSpeed(speed)
    {
        // The following has two functions. First, load the list of LV/RV nodes (1) and use to populate
        // the map so it can be used thusly: posn_in_dist_file = mMappingToXvDists.at(new_node_index)
        // Then does the same for each of the specified root nodes (3) into member var (4).
	std::cout << "i assign edge nodes" << std::endl;
        LoadFiles("projects/template_chaste/data/ARVC_005/lv_edge_nodes", mMappingToLvDists, rLvRootNodes, mLvRootNodeNewIndices );
        LoadFiles("projects/template_chaste/data/ARVC_005/rv_edge_nodes", mMappingToRvDists, rRvRootNodes, mRvRootNodeNewIndices );
        //LoadHeterogeneities("projects/pm2111/test/data/ARVC_005/fibrosis1.txt", heterogeneities,MappingOldIndexNewIndex);
    }

    ~DijkstraStimuliCellFactory()
    {
        fclose(mpLvDistancesFile);
        fclose(mpRvDistancesFile);
	dummy =0; 
	std::vector<int> mMappingToLvDists;
	std::vector<int> mMappingToRvDists;
	std::vector<int> MappingOldIndexNewIndex;
	std::vector<int> MapNewOldIndex;
	std::vector<int> mMapNewOldIndex;
	std::vector<int> regions;
	std::vector<int> mregions;
	std::vector<int> rNodePerm;
	std::vector<int> mrNodePerm;
    }

    AbstractCvodeCell* CreateCardiacCellForTissueNode(Node<3>* pNode)// method creates a cell for a tissue simulation
    {
	
        AbstractCvodeCell* p_cell;

	/*std::ostringstream strs;
	strs << node_index;
	std::string str = strs.str();
std::cout << "attempt node index "<< str  << std::endl ;
unsigned S2_time = 0;
double S2_r=1.0;
double x_S2=4.8; double y_S2= -9.1; double z_S2= -8.6;
ChastePoint<3> S2_centre (x_S2,y_S2,z_S2);
ChastePoint<3> S2_radius (S2_r,S2_r,S2_r);
ChasteEllipsoid<3> S2_region (S2_centre, S2_radius);
bool cell_is_in_S2Reg = S2_region.DoesContain(this->GetMesh()->GetNode(node_index)->rGetLocation());
double x = this->GetMesh()->GetNode(node_index)->rGetLocation()[0];
double y = this->GetMesh()->GetNode(node_index)->rGetLocation()[1];
double z = this->GetMesh()->GetNode(node_index)->rGetLocation()[2];
if(cell_is_in_S2Reg)
{
	
	
	boost::shared_ptr<SimpleStimulus> p_stim(new SimpleStimulus(-120000.0, 1.0, S2_time));
	p_cell->SetStimulusFunction(p_stim);
	

}
*/

        unsigned node_index = pNode->GetIndex(); // New index


        HeartGeometryInformation<3>* p_heart_geom_info = GetHeartGeometryInformation();
        const double lv_dist = p_heart_geom_info->rGetDistanceMapLeftVentricle()[node_index]; //get distance to the LV
        const double rv_dist = p_heart_geom_info->rGetDistanceMapRightVentricle()[node_index]; //get distance to the RV

        double min_distance = -1.0; // Nonsense





	//std::cout << "I start to compute the min distances for node " << str << std::endl;
  	
	/*Cell Assignment*/
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
	
	/*Cell Stimuli*/
      if (CompareDoubles::IsNearZero(lv_dist, 1e-3))
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
                    boost::shared_ptr<RegularStimulus> p_stim(new RegularStimulus(-120000.0, 1.0, mperiod,(double)t));
                    //boost::shared_ptr<SimpleStimulus> p_stim(new SimpleStimulus(-70000.0, 1.0, 0.5)); //stimulate at 0.5ms
                    mStimuli.push_back(p_stim);
                }
            }
            p_cell->SetStimulusFunction(mStimuli[upstroke_time]);//only pick the stimulus time as the time corresponding to the time it take for the signal to travel to it!!
        }


	//std::cout << "I start to vary Iks for node " << str << std::endl;
        const double gks = p_cell->GetParameter("membrane_slow_delayed_rectifier_potassium_current_conductance");
        const double gks_scaling = mpORdGksVarierDTI024->GetScalingFactor(*pNode);
        p_cell->SetParameter("membrane_slow_delayed_rectifier_potassium_current_conductance", gks_scaling*gks); //apply a scaling of Iks across the apexbase of the ventricle!

        p_cell->SetTolerances(1e-4,1e-6);

	//vary the sodium current conductance
	
	unsigned old_index = mMapNewOldIndex.at(node_index);
	double gna_scaling;

	
		gna_scaling = mregions[old_index][1];
	const double gina = p_cell->GetParameter("membrane_fast_sodium_current_conductance");
	p_cell-> SetParameter("membrane_fast_sodium_current_conductance",gna_scaling*gina);

        return p_cell;
	if (remainder(dummy,100000) ==0)

	{
		std::ostringstream strs;
		strs << dummy;
		std::string str = strs.str();
		std::cout << "Number of cells made : " <<  dummy << std::endl;
	}
	dummy = dummy++;
    
    }
};

#endif //DIJKSTRASTIMULICELLFACTORYORdPat02_HPP_
