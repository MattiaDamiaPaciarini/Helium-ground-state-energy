import numpy as np
import pylab as plt
import math
import random
from matplotlib import pyplot as plt
import array
import scipy
from scipy.constants import c

# Definitions of many functions 

# Definition of a distance funtion

def d(R,i,j):
  d = (R[0,i]-R[0,j])**2 + (R[1,i]-R[1,j])**2 +(R[2,i]-R[2,j])**2
  return d

# Definition of the Wavefunction

def WF(R, Par):
  wf = 1
  for i in range(0,A-1):
    for j in range(i+1,A):
      wf = wf * (math.exp(-Par[0]*d(R,i,j))+Par[1]*math.exp(-Par[2]*d(R,i,j)))
  #print("Eigenfunction value is = ", wf)     
  return wf 


# Definition and Computation of the Local Energy

def V(R):
  V = 0
  for i in range(0,A-1):
    for j in range(i+1,A):
      V = V + 1000*math.exp(-3*d(R,i,j))-165.35*math.exp(-1.05*d(R,i,j))-21.5*math.exp(-0.6*d(R,i,j))-83*math.exp(-0.8*d(R,i,j))-11.5*math.exp(-0.4*d(R,i,j))
  return V

# Derivative of WF, h is the increment

def Diff (R, i, j, h, Par):
  R_P = R.copy()
  R_M = R.copy()
  R_P[i,j] = R[i,j]+h
  R_M[i,j] = R[i,j]-h
  return (WF(R_P, Par)+WF(R_M, Par)-2*WF(R, Par))/h**2 

def K(R, Par):
  K = 0 
  for j in range(0,A):
    for i in range(0,3):
      K = K + Diff(R, i, j, h, Par )*C
  return K



# Some Initializations

NA = 0 # Accepted moves
NP = 300  # Number of moves to find parameters
NM = 10000  # Number of particles' moves
A = 4   # Numbers of nucleons
h = 0.00001  # Differential Increment 
step = 0.55  # Algorithm step for R
hbar = 6.582119*10**(-22)
C = (197.3269804)**2/2/939.56542052 # Constant in frontof the kinetic term, the first term is in Mev*fm
step1 = 1/100
step2 = 1/20
step3 = 1/10
list_gamma = []
list_a  = []
list_beta = []
list_E = []
list_EvalE = []

# Definition of the form factor function and its variational parameters

# gamma = 0.08597 
# a = -0.7191 
# beta = 2.13796 gamma+beta of Guardiola

# Parameters' list

BestPars = [0.08597 , -0.7191 , 2.13796 ]

Bestpars = np.array(BestPars)

Pars = [0.2, -0.2, 1]

Par = np.array(Pars)

Par_new = np.zeros(3)

# List for a particle's coordinates

R = (np.random.rand(3,A)-0.5)*10*step # Coordinates Matrix

# Coordinates list

list_R = []
list_R.append(R)

for k in range(0,NM):
  R_new = R.copy() + step*(np.random.rand(3,A)-0.5)
  if ((WF(R_new, Bestpars)/WF(R, Bestpars))**2>1 or (WF(R_new, Bestpars)/WF(R, Bestpars))**2 > np.random.rand() ):
    R = R_new.copy()
    list_R.append(R)
    NA = NA +1
  else:
    R = R.copy()
    list_R.append(R) 

print("Acceptance ratio =", 100*NA/NM)

for l in range(0,NP):
  Par_new[0] = Par[0].copy() + (np.random.rand()-0.5)*step1
  Par_new[1] = Par[1].copy() + (np.random.rand()-0.5)*step2
  Par_new[2] = Par[2].copy() + (np.random.rand()-0.5)*step3
  E_loc = 0
  D = 0
  for j in range(0,NM):
    E_loc = E_loc + V(list_R[j])*(WF(list_R[j], Par_new)/WF(list_R[j], Bestpars))**2 - K(list_R[j], Par_new)*WF(list_R[j], Par_new)/WF(list_R[j], Bestpars)**2 # Modify how the parameters arrive   
    D = D + (WF(list_R[j], Par_new)/WF(list_R[j], Bestpars))**2
    
  E_mean = E_loc/D

  if (l==0):
    Par = Par_new.copy()  
    list_gamma.append(Par[0])
    list_a.append(Par[1])
    list_beta.append(Par[2])
    list_E.append(E_mean)
    list_EvalE.append(E_mean)
    print("\nHI THERE, LET'S BEGIN!!", "\n", E_mean, "\nWith pars:","\nGamma = ",list_gamma[l],"\na = ",list_a[l],"\nbeta = ",list_beta[l])
  
  else:   
    print("l = ", l) 
    if (list_EvalE[l-1]>E_mean):
      Par = Par_new.copy()
      list_gamma.append(Par[0])
      list_a.append(Par[1])
      list_beta.append(Par[2])
      list_EvalE.append(E_mean)
      list_E.append(E_mean)
      print("Accepted!\nE_mean =", E_mean, "\nWith pars:","\nGamma = ",list_gamma[-1],"\na = ",list_a[-1],"\nbeta = ",list_beta[-1])
    else:
      Par = Par.copy()
      list_EvalE.append(list_EvalE[l-1])
      
      
print("last values:", "\nGamma = ",list_gamma[-1],"\na = ",list_a[-1],"\nbeta = ",list_beta[-1], "\nEnergy = ", list_E[-1])


fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(projection='3d')
cmap = plt.cm.cool
img = ax.scatter(list_gamma, list_a, list_beta, s = 10, c=list_E, vmin = -27, vmax = 0, cmap = cmap)
cbar = fig.colorbar(img, label= "Energy (MeV)", shrink = 0.5, orientation = 'horizontal')
ax.set_xlabel('gamma')
ax.set_ylabel('a')
ax.set_zlabel('beta')
plt.title('Energy Minimum', fontsize = 30)
