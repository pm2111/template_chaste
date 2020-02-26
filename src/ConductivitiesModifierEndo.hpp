#ifndef CONDUCTIVITIESMODIFIERANNULUSINDEXED3D_HPP_
#define CONDUCTIVITIESMODIFIERANNULUSINDEXED3DNEW_HPP_

#include "AbstractConductivityModifier.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "HeartGeometryInformation.hpp"


class ConductivitiesModifierEndo : public AbstractConductivityModifier<3,3>
{
private:
    const DistributedTetrahedralMesh<3,3>& mrMesh;
    // I'd like the following to be const too, but it would require making rGetDistanceMap*Ventricle() const
    //HeartGeometryInformation<2>& mrHeartGeometryInformation;
   // const std::vector<double>& mrDistMapRV;
    //const std::vector<double>& mrDistMapLV;
    c_matrix<double,3,3> mIntraFastTensor;
    c_matrix<double,3,3> mExtraFastTensor;
    c_matrix<double,3,3> mTensor;
    std::vector<unsigned>& mperm_vec;

    const std::vector<unsigned>& mendo; 
    std::map<unsigned,unsigned> MapNewOldIndex;

public:
    /*
     * Constructor.
     */
    /*ConductivitiesModifierAnnulus( HeartGeometryInformation<2>& rHeartGeomInfo,
                             const DistributedTetrahedralMesh<3,3>& rMesh)
                           : AbstractConductivityModifier<3,3>(),*/
    ConductivitiesModifierEndo(
                             const DistributedTetrahedralMesh<3,3>& rMesh,  const std::vector<unsigned>& endo,std::vector<unsigned>& perm_vec)
                           : AbstractConductivityModifier<3,3>(),
    mrMesh(rMesh),
    //mrHeartGeometryInformation(rHeartGeomInfo),
    //mrDistMapRV( rHeartGeomInfo.rGetDistanceMapRightVentricle() ),
    //mrDistMapLV( rHeartGeomInfo.rGetDistanceMapLeftVentricle() ),
    mIntraFastTensor( zero_matrix<double>(3,3) ),
    mExtraFastTensor( zero_matrix<double>(3,3)),
    mendo(endo),
    mperm_vec(perm_vec),
    mTensor( zero_matrix<double>(3,3) )
    {
        // Set constant tensors now
     //   for ( unsigned dim=0; dim<3; ++dim )
       // {
       //     mIntraFastTensor(dim,dim) = 3.0;
       //     mExtraFastTensor(dim,dim) = 2.03;
       // }
		mIntraFastTensor(0,0) = 1.5;
                mIntraFastTensor(1,1)=0.45;
                mIntraFastTensor(2,2)=0.225;
            mExtraFastTensor(0,0) = 3.64*1.5;
                mExtraFastTensor(1,1)= 2.03;
                mExtraFastTensor(2,2)=2.03;
 		std::cout <<" i finished allocating regions variable inside modifier class" << std::endl;
		for (unsigned int i=0; i<mperm_vec.size();i++)
		{			
			unsigned new_index = mperm_vec[i];
			MapNewOldIndex[new_index] =i;
		}
		        std::cout <<" i ve made a map" << std::endl;
              
}
   ~ConductivitiesModifierEndo()
    {

	std::vector<int> mrMesh;
	std::vector<int> rMesh;
	std::vector<int> mperm_vec;
	std::vector<int> perm_vec;
	std::vector<int> mregions;
	std::vector<int> MapNewOldIndex;


     }
  // mdomainIndex(        bool intracellular = CompareDoubles::WithinAnyTolerance( rOriginalConductivity(0,0), 1.5 ))

    /*
     * Method that returns modified tensor.
     */
    c_matrix<double,3,3>& rCalculateModifiedConductivityTensor(unsigned elementIndex,
                                                               const c_matrix<double,3,3>& rOriginalConductivity,
                                                               unsigned domainIndex);

};

#endif //FASTENDOCARDIUMMODIFIER_HPP_
