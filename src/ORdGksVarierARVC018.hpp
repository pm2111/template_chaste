#ifndef ORdGKSVARIER_HPP_
#define ORdGKSVARIER_HPP_

#include "HeartGeometryInformation.hpp"

class ORdGksVarierDTI024
{
private:
    HeartGeometryInformation<3>* mpHeartGeomInfo;
    bool mApexbase;
    bool mTransmur;
    bool mInterven;
    int mNormaliser;
    std::string mpath;

public:
    /*
     * Constructor. Takes a heart geometry info object and which gradients we're using.
     */
    ORdGksVarierDTI024(HeartGeometryInformation<3>* pHeartGeomInfo, bool ApexBase, bool Transmural, bool Interven,std::string path) :
    mpHeartGeomInfo(pHeartGeomInfo),
    mApexbase(ApexBase),
    mTransmur(Transmural),
    mpath(path),
    mInterven(Interven)
    {
        mNormaliser = (int)mApexbase+(int)mTransmur+(int)mInterven;
    }
    /*
     * Method to return the scaling factor.
     */
    double GetScalingFactor2(const Node<3>& rNode);
};

#endif //TT06GKSVARIER_HPP_
