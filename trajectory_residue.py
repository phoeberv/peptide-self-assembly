#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 16:49:34 2023

@author: phoeberose

TRAJECTORY FOR EACH RESIDUE

"""

import mdtraj as md
import numpy as np

FileName = 'QGAVIL_QGAVIL_final.pdb'
t = md.load(FileName)
L = t.unitcell_lengths[0][0]

topo = t.topology

"""
for i in range(1,7):
    for j in range(9,15):
        print("resid "+str(i)+"_"+"resid "+str(j))
"""

tResidue1 = t.atom_slice(topo.select("resid 1"))
tResidue2 = t.atom_slice(topo.select("resid 2"))
tResidue3 = t.atom_slice(topo.select("resid 3"))
tResidue4 = t.atom_slice(topo.select("resid 4"))                      
tResidue5 = t.atom_slice(topo.select("resid 5"))
tResidue6 = t.atom_slice(topo.select("resid 6"))

tResidue9 = t.atom_slice(topo.select("resid 9"))
tResidue10 = t.atom_slice(topo.select("resid 10"))
tResidue11 = t.atom_slice(topo.select("resid 11"))
tResidue12 = t.atom_slice(topo.select("resid 12"))
tResidue13 = t.atom_slice(topo.select("resid 13"))
tResidue14 = t.atom_slice(topo.select("resid 14"))

CM1 =  md.compute_center_of_mass(tResidue1, select=None)
CM2 =  md.compute_center_of_mass(tResidue2, select=None)
CM3 =  md.compute_center_of_mass(tResidue3, select=None)
CM4 =  md.compute_center_of_mass(tResidue4, select=None)
CM5 =  md.compute_center_of_mass(tResidue5, select=None)
CM6 =  md.compute_center_of_mass(tResidue6, select=None)
CM9 =  md.compute_center_of_mass(tResidue9, select=None)

CM10 =  md.compute_center_of_mass(tResidue10, select=None)
CM11 =  md.compute_center_of_mass(tResidue11, select=None)
CM12 =  md.compute_center_of_mass(tResidue12, select=None)
CM13 =  md.compute_center_of_mass(tResidue13, select=None)
CM14 =  md.compute_center_of_mass(tResidue14, select=None)


def distPeriodicBox(cm1,cm2,L):
    """

    Parameters
    ----------
    cm1 : TYPE
        DESCRIPTION.
    cm2 : TYPE
        DESCRIPTION.
    L : TYPE
        DESCRIPTION.

    Returns
    -------
    distStore : TYPE
        DESCRIPTION.

    """
    distStore = np.ones(cm1.shape[0])*np.inf
    (X,Y,Z) = np.meshgrid([-1,0,1],[-1,0,1],[-1,0,1])
    for i in range(3):
        for j in range(3):
            for k in range(3):   
                cm1_shift = 0+cm1
                cm1_shift[:,0]  += L*X[i,j,k]
                cm1_shift[:,1]  += L*Y[i,j,k]
                cm1_shift[:,2]  += L*Z[i,j,k]

                d = np.sqrt(np.sum((cm1_shift-cm2)**2,axis=1))
                distStore = np.min(np.array([d,distStore]),axis=0)
                
    return distStore
    

distanceBetweenRes1_9 = distPeriodicBox(CM1,CM9,L)
distanceBetweenRes1_10 = distPeriodicBox(CM1,CM10,L)
distanceBetweenRes1_11 = distPeriodicBox(CM1,CM11,L)
distanceBetweenRes1_12 = distPeriodicBox(CM1,CM12,L)
distanceBetweenRes1_13 = distPeriodicBox(CM1,CM13,L)
distanceBetweenRes1_14 = distPeriodicBox(CM1,CM14,L)

distanceBetweenRes2_9 = distPeriodicBox(CM2,CM9,L)
distanceBetweenRes2_10 = distPeriodicBox(CM2,CM10,L)
distanceBetweenRes2_11 = distPeriodicBox(CM2,CM11,L)
distanceBetweenRes2_12 = distPeriodicBox(CM2,CM12,L)
distanceBetweenRes2_13 = distPeriodicBox(CM2,CM13,L)
distanceBetweenRes2_14 = distPeriodicBox(CM2,CM14,L)


distanceBetweenRes3_9 = distPeriodicBox(CM3,CM9,L)
distanceBetweenRes3_10 = distPeriodicBox(CM3,CM10,L)
distanceBetweenRes3_11 = distPeriodicBox(CM3,CM11,L)
distanceBetweenRes3_12 = distPeriodicBox(CM3,CM12,L)
distanceBetweenRes3_13 = distPeriodicBox(CM3,CM13,L)
distanceBetweenRes3_14 = distPeriodicBox(CM3,CM14,L)

distanceBetweenRes4_9 = distPeriodicBox(CM4,CM9,L)
distanceBetweenRes4_10 = distPeriodicBox(CM4,CM10,L)
distanceBetweenRes4_11 = distPeriodicBox(CM4,CM11,L)
distanceBetweenRes4_12 = distPeriodicBox(CM4,CM12,L)
distanceBetweenRes4_13 = distPeriodicBox(CM4,CM13,L)
distanceBetweenRes4_14 = distPeriodicBox(CM4,CM14,L)

distanceBetweenRes5_9 = distPeriodicBox(CM5,CM9,L)
distanceBetweenRes5_10 = distPeriodicBox(CM5,CM10,L)
distanceBetweenRes5_11 = distPeriodicBox(CM5,CM11,L)
distanceBetweenRes5_12 = distPeriodicBox(CM5,CM12,L)
distanceBetweenRes5_13 = distPeriodicBox(CM5,CM13,L)
distanceBetweenRes5_14 = distPeriodicBox(CM5,CM14,L)

distanceBetweenRes6_9 = distPeriodicBox(CM6,CM9,L)
distanceBetweenRes6_10 = distPeriodicBox(CM6,CM10,L)
distanceBetweenRes6_11 = distPeriodicBox(CM6,CM11,L)
distanceBetweenRes6_12 = distPeriodicBox(CM6,CM12,L)
distanceBetweenRes6_13 = distPeriodicBox(CM6,CM13,L)
distanceBetweenRes6_14 = distPeriodicBox(CM6,CM14,L)


## SAVING FILENAMES

CMFileName1 = FileName[:-4]+'DistRes_1_9'
CMFileName2 = FileName[:-4]+'DistRes_1_10'
CMFileName3 = FileName[:-4]+'DistRes_1_11'
CMFileName4 = FileName[:-4]+'DistRes_1_12'
CMFileName5 = FileName[:-4]+'DistRes_1_13'
CMFileName6 = FileName[:-4]+'DistRes_1_14'

CMFileName7 = FileName[:-4]+'DistRes_2_9'
CMFileName8 = FileName[:-4]+'DistRes_2_10'
CMFileName9 = FileName[:-4]+'DistRes_2_11'
CMFileName10 = FileName[:-4]+'DistRes_2_12'
CMFileName11 = FileName[:-4]+'DistRes_2_13'
CMFileName12 = FileName[:-4]+'DistRes_2_14'

CMFileName13 = FileName[:-4]+'DistRes_3_9'
CMFileName14 = FileName[:-4]+'DistRes_3_10'
CMFileName15 = FileName[:-4]+'DistRes_3_11'
CMFileName16 = FileName[:-4]+'DistRes_3_12'
CMFileName17 = FileName[:-4]+'DistRes_3_13'
CMFileName18 = FileName[:-4]+'DistRes_3_14'

CMFileName19 = FileName[:-4]+'DistRes_4_9'
CMFileName20 = FileName[:-4]+'DistRes_4_10'
CMFileName21 = FileName[:-4]+'DistRes_4_11'
CMFileName22 = FileName[:-4]+'DistRes_4_12'
CMFileName23 = FileName[:-4]+'DistRes_4_13'
CMFileName24 = FileName[:-4]+'DistRes_4_14'

CMFileName25 = FileName[:-4]+'DistRes_5_9'
CMFileName26 = FileName[:-4]+'DistRes_5_10'
CMFileName27 = FileName[:-4]+'DistRes_5_11'
CMFileName28 = FileName[:-4]+'DistRes_5_12'
CMFileName29 = FileName[:-4]+'DistRes_5_13'
CMFileName30 = FileName[:-4]+'DistRes_5_14'

CMFileName31 = FileName[:-4]+'DistRes_6_9'
CMFileName32 = FileName[:-4]+'DistRes_6_10'
CMFileName33 = FileName[:-4]+'DistRes_6_11'
CMFileName34 = FileName[:-4]+'DistRes_6_12'
CMFileName35 = FileName[:-4]+'DistRes_6_13'
CMFileName36 = FileName[:-4]+'DistRes_6_14'


#SAVE ALL THE FILES

np.save(CMFileName1,distanceBetweenRes1_9)
np.save(CMFileName2,distanceBetweenRes1_10)
np.save(CMFileName3,distanceBetweenRes1_11)
np.save(CMFileName4,distanceBetweenRes1_12)
np.save(CMFileName5,distanceBetweenRes1_13)
np.save(CMFileName6,distanceBetweenRes1_14)

np.save(CMFileName7,distanceBetweenRes2_9)
np.save(CMFileName8,distanceBetweenRes2_10)
np.save(CMFileName9,distanceBetweenRes2_11)
np.save(CMFileName10,distanceBetweenRes2_12)
np.save(CMFileName11,distanceBetweenRes2_13)
np.save(CMFileName12,distanceBetweenRes2_14)

np.save(CMFileName13,distanceBetweenRes3_9)
np.save(CMFileName14,distanceBetweenRes3_10)
np.save(CMFileName15,distanceBetweenRes3_11)
np.save(CMFileName16,distanceBetweenRes3_12)
np.save(CMFileName17,distanceBetweenRes3_13)
np.save(CMFileName18,distanceBetweenRes3_14)

np.save(CMFileName19,distanceBetweenRes4_9)
np.save(CMFileName20,distanceBetweenRes4_10)
np.save(CMFileName21,distanceBetweenRes4_11)
np.save(CMFileName22,distanceBetweenRes4_12)
np.save(CMFileName23,distanceBetweenRes4_13)
np.save(CMFileName24,distanceBetweenRes4_14)

np.save(CMFileName25,distanceBetweenRes5_9)
np.save(CMFileName26,distanceBetweenRes5_10)
np.save(CMFileName27,distanceBetweenRes5_11)
np.save(CMFileName28,distanceBetweenRes5_12)
np.save(CMFileName29,distanceBetweenRes5_13)
np.save(CMFileName30,distanceBetweenRes5_14)

np.save(CMFileName31,distanceBetweenRes6_9)
np.save(CMFileName32,distanceBetweenRes6_10)
np.save(CMFileName33,distanceBetweenRes6_11)
np.save(CMFileName34,distanceBetweenRes6_12)
np.save(CMFileName35,distanceBetweenRes6_13)
np.save(CMFileName36,distanceBetweenRes6_14)


