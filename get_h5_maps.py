import h5py
import numpy as np
import os
import sys
#sys.path.insert(0,"/home/pr1efz00/pr1efz06/chaste-libs/lib/python2.7/site-packages")
print sys.path
from joblib import Parallel,delayed


def my_fcn(i):

	results = h5py.File('results.h5', 'r')
	
#	atmap = results['/UpstrokeTimeMap_0'][()]
	a,N,b = atmap.shape
	voltages=results['/Data'][i,0:N,0]
        #print(min(voltages))
        #print(max(voltages))
        #print(voltages.shape)

        f=open("voltages"+".bin"+str(i),'w+b')
        binary_format=bytearray(voltages)
        f.write(binary_format)
        f.close()

	
	
	return 0


os.chdir('/data/blanca-rodriguez/shug5389/testoutput/TestHeart09LGEStrong/')
results = h5py.File('results.h5', 'r')
atmap = results['/UpstrokeTimeMap_0'][()]
a,N,b = atmap.shape
Parallel(n_jobs=100)(map(delayed(my_fcn),range(800,1000,1)))#results['/Data'][1,0:N,0]))


#print(voltages.shape)
#apdmap = iresults['/Apd_90_0_Map'][()]
#results.close()

maps = h5py.File('results_maps.h5', 'w')
param = np.empty_like(atmap)
param[:] = atmap

#param1 = np.empty_like(voltages)
#param1[:]=voltages
maps.create_dataset('/UpstrokeTimeMap_0', data=param)
#maps.create_dataset('Voltages',data=param1)
#param = numpy.empty_like(apdmap)
#param[:] = apdmap
#maps.create_dataset('/Apd_90_0_Map', data=param)
maps.close()
