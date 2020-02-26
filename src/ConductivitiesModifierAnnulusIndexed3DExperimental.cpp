#include "ConductivitiesModifierAnnulusIndexed3DExperimental.hpp"
#include <algorithm>    // std::find
#include <cmath>        // std::abs
#include "cell_based/src/common/SimulationTime.hpp"
/*
bool inSEG07( double x , double y , double z ){

  //SEGMENT 07
  double xm, ym, zm, xr, yr, zr;
  xm = x - (6.67902213614625);
  ym = y - (-5.00772451688348);
  zm = z - (-1.52547143989758);
  xr = xm * (-0.318059449084858) + ym * (0.787748468897617) + zm * (0.527532498143287);
  yr = xm * (-0.804220072437278) + ym * (0.0704877252690338) + zm * (-0.590136895707578);
  zr = xm * (-0.502064001833248) + ym * (-0.611950839803213) + zm * (0.611103843652883);
  xr = xr / (2.41460114345428);
  yr = yr / (1.77879695152772);
  zr = zr / (0.860150719845657);
  return ( xr * xr + yr * yr + zr * zr )<1;
}

bool inSEG11( double x , double y , double z ){

  //SEGMENT 11
  double xm, ym, zm, xr, yr, zr;
  xm = x - (6.91316708907576);
  ym = y - (-8.2209571500698);
  zm = z - (-3.12940263113258);
  xr = xm * (0.538484216000832) + ym * (0.685218196002451) + zm * (0.490418976983064);
  yr = xm * (0.715047580806098) + ym * (-0.679500108111716) + zm * (0.164275866333165);
  zr = xm * (-0.44580456065558) + ym * (-0.262212941982869) + zm * (0.85586369636489);
  xr = xr / (2.67874067627243);
  yr = yr / (1.41421770847719);
  zr = zr / (0.874942132792125);
  return ( xr * xr + yr * yr + zr * zr )<1;
}


bool inSEG16( double x , double y , double z ){

  //SEGMENT 16
  double xm, ym, zm, xr, yr, zr;
  xm = x - (-3.11559147977939);
  ym = y - (-9.32516699946728);
  zm = z - (-2.06152833800312);
  xr = xm * (-0.842246697392823) + ym * (0.452832096284263) + zm * (0.292512552389948);
  yr = xm * (0.0235749791878068) + ym * (-0.511144244377954) + zm * (0.85917156714802);
  zr = xm * (-0.5385765693819) + ym * (-0.730530392259002) + zm * (-0.41983404447321);
  xr = xr / (3.88671876575002);
  yr = yr / (2.7611104834997);
  zr = zr / (0.989126696522687);
  return ( xr * xr + yr * yr + zr * zr )<1;
}

bool inSEG17( double x , double y , double z ){
//SEGMENT 17
  double xm, ym, zm, xr, yr, zr;
  xm = x - (0.860061707872041);
  ym = y - (-10.5016068905292);
  zm = z - (-3.36963685047205);
  xr = xm * (-0.976761271794726) + ym * (0.200872018655396) + zm * (0.0747519233415108);
  yr = xm * (0.00390360153923988) + ym * (-0.332039384663966) + zm * (0.943257445731014);
  zr = xm * (-0.214294609864485) + ym * (-0.921629144045088) + zm * (-0.3235390873291);
  xr = xr / (2.88352555903292);
  yr = yr / (2.07585901444624);
  zr = zr / (0.617161840775838);
  return ( xr * xr + yr * yr + zr * zr )<1;
}



bool inSEG21( double x , double y , double z ){

  //SEGMENT 21
  double xm, ym, zm, xr, yr, zr;
  xm = x - (4.76030861779894);
  ym = y - (-10.32538132068);
  zm = z - (-4.05638269450689);
  xr = xm * (-0.972331014338916) + ym * (-0.216403900739746) + zm * (-0.0879872166810481);
  yr = xm * (0.095335538959528) + ym * (-0.711441378793548) + zm * (0.696248733967705);
  zr = xm * (-0.213268688667492) + ym * (0.668595929007167) + zm * (0.712387499994975);
  xr = xr / (2.59254678500536);
  yr = yr / (1.39115495879548);
  zr = zr / (0.57271600144441);
  return (xr * xr + yr * yr + zr * zr)<1 ;
}
*/


c_matrix<double,3,3>& ConductivitiesModifierAnnulus3D::
                      rCalculateModifiedConductivityTensor(unsigned elementIndex,
                                                           const c_matrix<double,3,3>& rOriginalConductivity,
                                                           unsigned domainIndex)
{


        //std::cout <<" i get to conductivities modifier cpp" << std::endl;
	double conductivity_scaling_factor;
	//std::cout <<" i get to the conductivities changer for node" <<  std::endl;
	//if any of the nodes in the element are fibrotic, change the conductivities of the whole element!
	bool   fibrosis = true;
    	Element<3, 3>* p_element = mrMesh.GetElement(elementIndex);
	//std::cout <<" i get the element index " << elementIndex <<  std::endl;
	int i=0;

//		Node<3> *Nodo;	
    	//while (i<4 )
    	//{


		//try {
		
//			Nodo = mrMesh.GetNode( p_element->GetNodeGlobalIndex(i) );
			unsigned new_index =p_element->GetNodeGlobalIndex(i);
			
//Nodo =  p_element->GetNode(i) ;	// this is a local node        
			//x =(double) Nodo->GetPoint()[0];
			//y = (double) Nodo->GetPoint()[1];

//			double x = Nodo->rGetLocation()[0];
//			double y = Nodo->rGetLocation()[1];
//			double z = Nodo->rGetLocation()[2];

			//std::cout << "coords chaste" << x << " " << y << " " << z << std::endl;
			//std::cout << "coords from regions" << mregions[i][0] << " " << mregions[i][1] << " " << mregions[i][2] <<std::endl;
			//std::cout << x  << "  " << y << std::endl;
		  //  }
		//catch (...) {
			//i++;
			//Nodo = mrMesh.GetNode( p_element->GetNodeGlobalIndex(i) );
		     //}
			


			//unsigned new_index = Nodo->GetIndex();
			//std::cout <<" i get new index " <<  std::endl;		
			//unsigned old_index = mMapNewOldIndex.at(new_index);
			//std::cout << old_index <<  std::endl;		
			//std::cout <<" i get old index " <<  std::endl;		
				//std::cout <<" i get to the conductivities changer for node" << new_index << "old index" << old_index << std::endl;
				/*for (int j=0; j<mregions.size(); j++)
				{
			
					//if ( (std::abs(mfibrosis_x_coord[j] - x) < 0.000001)) //&& ( std::abs(mfibrosis_y_coord[j] - y) < 0.01) ) //something is not working with the y-coordinate
					if ( (std::abs(mregions[j][0] - x) < 0.001) && ( std::abs(mregions[j][1] - y) < 0.001) && ( std::abs(mregions[j][2] - z) < 0.001) ) //something is not working with the y-coordinate				
					{*/
						
						//fibrosis =true;
						conductivity_scaling_factor =mregions[MapNewOldIndex.at(new_index)][3];
						//std::cout << x << " " << y << " "<< z << " " <<conductivity_scaling_factor << std::endl;
						//std::cout << "fibrotic index " << Nodo->GetIndex() << "conductivity scaling " << mregions[j][3] << std::endl;
						//break;
						
					
					//}
					
//}
				/*if (fibrosis==true)
				{
					break;
				}*/

			



    		//x = x + (double) Nodo->GetPoint()[0];
        	//y = y + (double) Nodo->GetPoint()[1];
        	//z = z + (double) Nodo->GetPoint()[2];
		i++;
	
	//}
	
//x = x/3.0;
       // y = y/3.0; //find position of center of the triangle
       // z = z/4.0;

        // Figure out region to determine whether this should be fast conducting
        //double rv_dist = mrDistMapRV.at(global_index);
        //double lv_dist = mrDistMapLV.at(global_index);
	
	



	/*if( x  > 0.0 )
        {
//MARK;
           	fibrosis = true;
		
	//	FILE *fid;
    	//	fid = fopen( "points_in_scar_07_simple.xyz" , "a" );
    	//	fprintf( fid , "%0.7g %0.7g %0.7g\n", x , y , z );
    	//	fclose( fid );
       }*/
   

    if ( fibrosis==true )
    {	
	//change the conductivities tensors according to lge data for that segment
	mIntraFastTensor(0,0) = 1.5*conductivity_scaling_factor;
                mIntraFastTensor(1,1)=0.45*conductivity_scaling_factor;
                mIntraFastTensor(2,2)=0.225*conductivity_scaling_factor;
            mExtraFastTensor(0,0) = 3.64*1.5*conductivity_scaling_factor;
                mExtraFastTensor(1,1)= 2.03*conductivity_scaling_factor;
                mExtraFastTensor(2,2)=2.03*conductivity_scaling_factor;
        return domainIndex ? mExtraFastTensor : mIntraFastTensor; // 0=intra, 1=extra
    }
    else
    {
        // Return the original one
        for ( unsigned dim=0; dim<3; ++dim )
        {
            mTensor(dim,dim) = rOriginalConductivity(dim,dim);
        }
        return mTensor;
    }

}

