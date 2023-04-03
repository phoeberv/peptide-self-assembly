#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  7 17:17:54 2023

@author: phoeberose
"""

import numpy as np
from matplotlib.pyplot import *
import matplotlib.pyplot as plt
from scipy.signal import find_peaks



#LOAD DATA
pep1 = np.load ('HGAVIL_HGAVIL_finalDistCM.npy')
pep1_2 = np.load ('HGAVIL_HGAVIL_finalDistCM_2.npy')
pep1_3 = np.load ('HGAVIL_HGAVIL_finalDistCM_3.npy')

pep2 = np.load ('KGAVIL_KGAVIL_finalDistCM.npy')
pep2_2 = np.load ('KGAVIL_KGAVIL_finalDistCM_2.npy')
pep2_3 = np.load ('KGAVIL_KGAVIL_finalDistCM_3.npy')

pep3 = np.load ('RGAVIL_RGAVIL_finalDistCM.npy')
pep3_2 = np.load ('RGAVIL_RGAVIL_finalDistCM_2.npy')
pep3_3 = np.load ('RGAVIL_RGAVIL_finalDistCM_3.npy')

pep4 = np.load ('GGAVIL_GGAVIL_finalDistCM.npy')
pep4_2 = np.load ('GGAVIL_GGAVIL_finalDistCM_2.npy')
pep4_3 = np.load ('GGAVIL_GGAVIL_finalDistCM_3.npy')

pep8 = np.load ('EGAVIL_EGAVIL_finalDistCM.npy')
pep8_2 = np.load ('EGAVIL_EGAVIL_finalDistCM_2.npy')
pep8_3 = np.load ('EGAVIL_EGAVIL_finalDistCM_3.npy')

pep9 = np.load ('QGAVIL_QGAVIL_finalDistCM.npy')
pep9_2 = np.load ('QGAVIL_QGAVIL_finalDistCM_2.npy')
pep9_3 = np.load ('QGAVIL_QGAVIL_finalDistCM_3.npy')

#print average distance
print('Pep1: ' + str(np.mean(pep1)))
print('Pep1: ' + str(np.mean(pep1_2)))
print('Pep1: ' + str(np.mean(pep1_3)))

print('Pep2: ' + str(np.mean(pep2)))
print('Pep2: ' + str(np.mean(pep2_2)))
print('Pep2: ' + str(np.mean(pep2_3)))

print('Pep3: ' + str(np.mean(pep3)))
print('Pep3: ' + str(np.mean(pep3_2)))
print('Pep3: ' + str(np.mean(pep3_3)))

print('Pep4: ' + str(np.mean(pep4)))
print('Pep4: ' + str(np.mean(pep4_2)))
print('Pep4: ' + str(np.mean(pep4_3)))

print('Pep8: ' + str(np.mean(pep8)))
print('Pep8: ' + str(np.mean(pep8_2)))
print('Pep8: ' + str(np.mean(pep8_3)))


print('Pep9: ' + str(np.mean(pep9)))
print('Pep9: ' + str(np.mean(pep9_2)))
print('Pep9: ' + str(np.mean(pep9_3)))






#make histogram
#combined histogram

X = plt.hist(np.hstack([pep3,pep3_2]),bins=200,range=[0,4],density=True)
plt.ylim(0,2.5)
plt.xlabel('$r_{cm}$ (nm)' )
plt.ylabel('$f~(r_{cm})$')
resolution_value = 600
plt.savefig('Hist_PepP9_y2.5.png', format="png", dpi=resolution_value)
plt.show()

c = 0.5;
f = X[0]
r = X[1]
prob_less_than_c = sum(f[np.nonzero(r<c)[0]])*(r[1]-r[0])

print(prob_less_than_c)

plot(pep1)

"""
from matplotlib.pyplot import plot
X
from matplotlib.pyplot import plot
X[0]
X[A]
X[1]
plot(X[1],X[0])
X[1]
x = (X[1][:-1]+X[1][1:])/2
x
plot(x,X[0])
"""
