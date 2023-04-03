#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 16:53:58 2023

@author: phoeberose
"""

import mdtraj as md
import numpy as np

FileName = 'HGAVIL_HGAVIL_final.pdb'
t = md.load(FileName)
L = t.unitcell_lengths[0][0]

topo = t.topology


tPeptide1 = t.atom_slice(topo.select("resid 1 2 3 4 5 6"))
tPeptide2 = t.atom_slice(topo.select("resid 9 10 11 12 13 14"))

CM1 =  md.compute_center_of_mass(tPeptide1, select=None)
CM2 =  md.compute_center_of_mass(tPeptide2, select=None)



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
    

distanceBetweenCM = distPeriodicBox(CM1,CM2,L)


CMFileName = FileName[:-4]+'DistCM'
np.save(CMFileName,distanceBetweenCM)
    
    


"""


md.compute_distances(t, [[40,121], [20,165]], periodic=True)





# To check atom name
#tPeptide1.topology.atom(109)

print (md.compute_phi(t))

#md.compute_distances(traj, atom_pairs, periodic=True, opt=True)

#$am_pairs=np.ndarray, shape=(num_pairs, 2), dtype=int

md.compute_center_of_mass(t, select=None)

#print center of mass
print(md.compute_center_of_mass(t, select=None))

#print 
print(md.compute_contacts(traj, contacts='all', scheme='closest-heavy', ignore_nonprotein=True, periodic=True, soft_min=False, soft_min_beta=20))

"""