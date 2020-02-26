#ifndef TypeAsigment_stimulus_HPP_
#define TypeAsigment_stimulus_HPP_

#include <algorithm>
#include <fstream>

#include "Tomek_model13endo_fkatpCvodeOpt.hpp"
#include "Tomek_model13mid_fkatpCvodeOpt.hpp"
#include "Tomek_model13epi_fkatpCvodeOpt.hpp"

class TypeAsigment_stimulusDTI004 : public AbstractCardiacCellFactory<3>
{
private:
    ORdGksVarierDTI004* mpORdGksVarier;
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
    TypeAsigment_stimulusDTI004(ORdGksVarierDTI004* pGksVarier,
                              const std::vector<unsigned>& rNodePerm,
                              const std::vector<unsigned>& rLvRootNodes,
                              const std::vector<unsigned>& rRvRootNodes,
                              double speed)
    : AbstractCardiacCellFactory<3>(),
    mpORdGksVarier(pGksVarier),
    mrNodePerm(rNodePerm),
    mpLvDistancesFile(fopen("/gpfs/projects/pr1efz00/pr1efz10/DTI004Mesh/lv.mindists","rb")),
    mpRvDistancesFile(fopen("/gpfs/projects/pr1efz00/pr1efz10/DTI004Mesh/rv.mindists","rb")),
    mSpeed(speed)
    {
        // The following has two functions. First, load the list of LV/RV nodes (1) and use to populate
        // the map so it can be used thusly: posn_in_dist_file = mMappingToXvDists.at(new_node_index)
        // Then does the same for each of the specified root nodes (3) into member var (4).
        LoadFiles("/gpfs/projects/pr1efz00/pr1efz10/DTI004Mesh/lv_edge_nodes", mMappingToLvDists, rLvRootNodes, mLvRootNodeNewIndices );
        LoadFiles("/gpfs/projects/pr1efz00/pr1efz10/DTI004Mesh/rv_edge_nodes", mMappingToRvDists, rRvRootNodes, mRvRootNodeNewIndices );
    }

    ~TypeAsigment_stimulusDTI004()
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

	if (relative_position > 0.70)
        {
                p_cell = new CellTomek_model13epi_fkatpFromCellMLCvodeOpt(mpSolver, mpZeroStimulus);
        }
        else if(relative_position>0.45)
        {
                p_cell = new CellTomek_model13mid_fkatpFromCellMLCvodeOpt(mpSolver, mpZeroStimulus);
        }
        else
        {
                p_cell = new CellTomek_model13endo_fkatpFromCellMLCvodeOpt(mpSolver, mpZeroStimulus);
        }
	
/*
//S2 stimulus
unsigned S2_time = 0;
double S2_r=1.0;
double x_S2=6.004576; double y_S2= -4.308411; double z_S2= -2.004203;
ChastePoint<3> S2_centre (x_S2,y_S2,z_S2);
ChastePoint<3> S2_radius (S2_r,S2_r,S2_r);
ChasteEllipsoid<3> S2_region (S2_centre, S2_radius);
bool cell_is_in_S2Reg = S2_region.DoesContain(this->GetMesh()->GetNode(node_index)->rGetLocation());
double x = this->GetMesh()->GetNode(node_index)->rGetLocation()[0];
double y = this->GetMesh()->GetNode(node_index)->rGetLocation()[1];
double z = this->GetMesh()->GetNode(node_index)->rGetLocation()[2];
if(cell_is_in_S2Reg)
{
	if(relative_position==0.0) 
	{
	boost::shared_ptr<SimpleStimulus> p_stim(new SimpleStimulus(-120000.0, 1.0, S2_time));
	p_cell->SetStimulusFunction(p_stim);
	}

}
*/


/*	double isch_r = 3.0;
        double isch_bz_r = 3.5;
        double x_isch = 3.2;double y_isch = 10.3; double z_isch=3.6;
        ChastePoint<3> isch_centre (x_isch,y_isch,z_isch);
        ChastePoint<3> isch_radius (isch_r,isch_r,isch_r);
        ChastePoint<3> isch_bz_radius (isch_bz_r,isch_bz_r,isch_bz_r);
        ChasteEllipsoid<3> isch_region (isch_centre,isch_radius);
        ChasteEllipsoid<3> isch_bz_region (isch_centre,isch_bz_radius);
        bool cell_is_in_ischReg = isch_region.DoesContain(this->GetMesh()->GetNode(node_index)->rGetLocation());
        bool cell_is_in_isch_BZReg = isch_bz_region.DoesContain(this->GetMesh()->GetNode(node_index)->rGetLocation());
	
        double x = this->GetMesh()->GetNode(node_index)->rGetLocation()[0];
        double y = this->GetMesh()->GetNode(node_index)->rGetLocation()[1];
        double z = this->GetMesh()->GetNode(node_index)->rGetLocation()[2];
	// ---------- Ischaemic parameters
	        double k_n = 5.4;
	        double k_i = 9.0;
	        double fkatp_n = 0.0;
	        double fkatp_i = 0.05;
	        double ina_n = 1.0;
        	double ina_i = 0.75;

	if(cell_is_in_isch_BZReg) 
	{

	if(relative_position>0.1)   // Endocardial BZ
        {
	
                if(cell_is_in_ischReg)
                {
                   // Cell is in central ischaemic region
                       p_cell->SetParameter("extracellular_potassium_concentration",k_i);
                       double gna = p_cell->GetParameter("membrane_fast_sodium_current_conductance");
                       double gca = p_cell->GetParameter("membrane_L_type_calcium_current_conductance");
                       p_cell->SetParameter("membrane_fast_sodium_current_conductance",gna*ina_i);
                       p_cell->SetParameter("membrane_L_type_calcium_current_conductance",gca*ina_i);
                       p_cell->SetParameter("membrane_atp_dependent_potassium_current_conductance",fkatp_i);
                 }

		else
	        {
                 // cell is in BZ area
                 //Get BZ properties
                 // we assume INa and ICaL BZ are 1/2 of K+ BZ and that fkatp BZ is a tenth of K+ BZ
                 double distance_from_isch_centre = sqrt(pow((x-x_isch),2)+pow((y-y_isch),2)+pow((z-z_isch),2)); // distance from ischaemic centre
                        // distance from central ischaemic zone
	                        // 0 if at central ischaemic zone
	                        // 0.5 if at normal zone
	                        double distance_from_cz = distance_from_isch_centre - isch_r;
	
	
	                        double bz_k = isch_bz_r-isch_r; // K+ border zone width
	//                        double bz_ina = 0.5*bz_k; // INa and ICaL border zone width
	//                        double distance_from_cz_ina = distance_from_isch_centre - (isch_r; // distance from central ischaemic zone when INa BZ starts
	//                        double bz_fkatp = 0.1*bz_k; // IK(ATP) border zone width
	//                        double distance_from_cz3 = distance_from_cz - bz_fkatp; // distance from central ischaemic zone when fkatp BZ starts
	
	                        // linear ax+b equations, b=ischaemic value
	                        double a_k = (k_n-k_i)/bz_k;
	                        double a_ina = (ina_n-ina_i)/bz_k;
	                        double a_fkatp = (fkatp_n-fkatp_i)/bz_k;
	
	                        // ----- calculate y
	                                    double y_k = a_k*distance_from_cz+k_i; // solution to ax+b=y equation
	                                    double y_ina = a_ina*distance_from_cz+ina_i;
	                                    double y_fkatp = a_fkatp*distance_from_cz+fkatp_i;
	                                    p_cell->SetParameter("extracellular_potassium_concentration",y_k);
	                        double gna = p_cell->GetParameter("membrane_fast_sodium_current_conductance");
	                        double gca = p_cell->GetParameter("membrane_L_type_calcium_current_conductance");
	                        p_cell->SetParameter("membrane_fast_sodium_current_conductance",gna*y_ina);
	                        p_cell->SetParameter("membrane_L_type_calcium_current_conductance",gca*y_ina);
	                        p_cell->SetParameter("membrane_atp_dependent_potassium_current_conductance",y_fkatp);
	
                                }
}
}*/
	if ( min_distance >= 0.0)
        {
            // Convert this into a time
            unsigned upstroke_time = min_distance/mSpeed + 0.5; // Result is rounded! In ms.
            // Have we got a stimulus for this time? If not, extend.
            if( upstroke_time >= mStimuli.size() )
            {
                for(unsigned t=mStimuli.size(); t<=upstroke_time; ++t)
                {
                    boost::shared_ptr<RegularStimulus> p_stim(new RegularStimulus(-120000.0, 1.0, 1000,(double)t));
                    mStimuli.push_back(p_stim);
                }
            }
	    
		p_cell->SetStimulusFunction(mStimuli[upstroke_time]);
	//	p_cell->SetStimulusFunction(mStimuli[upstroke_time+1000]);
        }

       const double gks = p_cell->GetParameter("membrane_slow_delayed_rectifier_potassium_current_conductance");
       const double gks_scaling = mpORdGksVarier->GetScalingFactor(*pNode);
       p_cell->SetParameter("membrane_slow_delayed_rectifier_potassium_current_conductance", gks_scaling*gks);

        p_cell->SetTolerances(1e-4,1e-6);

        return p_cell;
    }
};

#endif //TypeAsigment_stimulusHPP_
