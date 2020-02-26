#include "ORdGksVarierARVC05.hpp"

/*
 * Method to return the scaling factor.
 */
double ORdGksVarierDTI024::GetScalingFactor(const Node<3>& rNode)
{
    if (mNormaliser==0) return 1.0; // No heterogeneity

    unsigned node_index = rNode.GetIndex();
    HeartRegionType this_region = mpHeartGeomInfo->GetHeartRegion(node_index);
    double gks_scaling = 0.0;
    // Decrease APD at apex and vice versa. Some values are hard-coded for speed.
    if (mApexbase)
    {
        const double x = rNode.GetPoint()[0];
        const double y = rNode.GetPoint()[1];
	const double z = rNode.GetPoint()[2];
	

        // Gks varier (apicobasal, transmural, interventricular)
	std::vector<double> apexbase;    // a vector to hold apexbase data
	// open file    
	std::ifstream inputFile("projects/template_chaste/data/ARVC_005/apexbase");

	// test file open   
	if (inputFile) {       
	double value;
	
	// read the elements in the file into a vector  
	while ( inputFile >> value ) {
	apexbase.push_back(value);
	}
	}
	const double Ax = apexbase[0];
	const double Ay = apexbase[1];
	const double Az = apexbase[2];
	const double Bx = apexbase[3];
	const double By = apexbase[4];
	const double Bz = apexbase[5];
/*const double  Ax= 7.6037445068359357e+00;
const double  Ay= -4.5340534210205075e+00;
const double  Az= -3.4286319732666009e+00;

const double  Bx= 3.1061785340784427e+00;
const double  By= -1.1930589407185264e+00;
const double  Bz= 5.6445072311608602e-01;*/


const double BAx= Bx-Ax;
const double BAy= By-Ay;
const double BAz= Bz-Az;

const double n2BA= BAx*BAx + BAy*BAy + BAz*BAz;
const double Z= ( (x-Ax)*BAx + (y-Ay)*BAy + (z-Az)* BAz )/ n2BA;	

        const double coord = 2.0*Z - 1.0;
        /* AB uses s=0.7160 for a 42 ms range from default */
        gks_scaling += pow(0.2, coord);
    }
    // Increase APD at endocardium and vice versa.
    if (mTransmur)
    {
        const double coord = 1.0-2.0*mpHeartGeomInfo->CalculateRelativeWallPosition(node_index);
        /* TM uses s=0.7160 for a 42 ms range from default */
        gks_scaling += pow(0.15, coord);
    }
    // Increase the scaling factor in the RV (shorter APD) and decrease it in the LV (longer APD).
    if (mInterven)
    {
        /* IV uses s=0.7753 for a 32 ms range from default */
        if ( this_region == mpHeartGeomInfo->RIGHT_VENTRICLE_WALL )
        {
            gks_scaling += 1.0/0.2;
        }
        else
        {
            gks_scaling += 0.2;
        }
    }
    return gks_scaling/mNormaliser;
}
