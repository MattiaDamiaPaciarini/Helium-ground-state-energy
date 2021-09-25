import numpy as np
import pylab as plt
import math
import random
from matplotlib import pyplot as plt
import array
import scipy
from scipy.constants import c

# Definitions of many functions 

# center of mass

def CM (R):
  x_cm = 0 
  y_cm = 0 
  z_cm = 0
  for k in range(0,A):
    x_cm = x_cm + R[0,k]
    y_cm = y_cm + R[1,k]
    z_cm = z_cm + R[2,k]
  return x_cm/A, y_cm/A, z_cm/A

# Definition of a distance funtion

def d(R,i,j):
  d = (R[0,i]-R[0,j])**2 + (R[1,i]-R[1,j])**2 +(R[2,i]-R[2,j])**2
  return d

# Definition of the Wavefunction

def WF(R, gamma, a, beta):
  wf = 1
  for i in range(0,A-1):
    for j in range(i+1,A):
      wf = wf * (math.exp(-gamma*d(R,i,j))+a*math.exp(-beta*d(R,i,j)))    
  return wf 


# Definition and Computation of the Local Energy

def V(R):
  V = 0
  for i in range(0,A-1):
    for j in range(i+1,A):
      V = V + 1000*math.exp(-3*d(R,i,j))-165.35*math.exp(-1.05*d(R,i,j))-21.5*math.exp(-0.6*d(R,i,j))-83*math.exp(-0.8*d(R,i,j))-11.5*math.exp(-0.4*d(R,i,j))
  return V

# Derivative of WF, h is the increment

def Diff (R, i, j, h):
  R_P = R.copy()
  R_M = R.copy()
  R_P[i,j] = R[i,j]+h
  R_M[i,j] = R[i,j]-h
  return (WF(R_P, gamma, a, beta)+WF(R_M, gamma, a, beta)-2*WF(R, gamma, a, beta))/h**2 

def K(R):
  K = 0 
  for j in range(0,A):
    for i in range(0,3):
      K = K + Diff(R,i,j,h)*C
  return K



# Some Initializations

Len = 1200   # Number of catched points
NM = 12000 # Number of moves
NA = 0  # Number of acceptances
A = 4   # Numbers of nucleons
h = 0.00001  # Differential Increment 
t = 10  # Thermalizer
step = 0.6  # Algorithm step
hbar = 6.582119*10**(-22)
C = (197.3269804)**2/2/939.56542052 # Constant in frontof the kinetic term, the first term is in Mev*fm
E_loc = 0
SE = 0

# Definition of the form factor function and its variational parameters 

# Optimal values

gamma = 0.08597
a =  -0.7191
beta =  2.13796 #gamma+beta of Guardiola

# List for a particle's coordinates

list_X = [] # Respect to an external coordinates system
list_Y = []
list_X_cm = [] # Respect to cm
list_Y_cm = []
list_Z_cm = []
list_X_cm_t = [] # Wrt the CM and thermalized
list_Y_cm_t = []
list_X_t = [] # Respect to an external coordinates system
list_Y_t = []
list_X_cm_t_2 = [] # Second particle's coordinates wrt CM and thermalized
list_Y_cm_t_2 = []
list_X_cm_t_3 = [] # Third particle's coordinates wrt CM and thermalized
list_Y_cm_t_3 = []
list_X_cm_t_4 = [] # Fourth particles' coordinates wrt CM and thermalized
list_Y_cm_t_4 = []

# Generation of the Position Matrix

R = (np.random.rand(3,A)-0.5)*10*step # Coordinates Matrix

# Metropolis Algorithm

for k in range(0,NM):

# Generation of a new configuration

  R_new = R.copy() + step*(np.random.rand(3,A)-0.5)

# Acceptance or not?
  
  if ((WF(R_new, gamma, a, beta)/WF(R, gamma, a, beta))**2>1 or (WF(R_new, gamma, a, beta)/WF(R, gamma, a, beta))**2 > np.random.rand() ):
    R = R_new.copy()
    NA = NA + 1
    if (len(list_X)<Len):
      list_X.append(R[0,0]) # Append new point in the list
      list_Y.append(R[1,0])
      list_X_cm.append(CM(R)[0]-R[0,0]) # Append new position wrt CM
      list_Y_cm.append(CM(R)[1]-R[1,0])
      list_Z_cm.append(CM(R)[2]-R[2,0])
    if (k%t == 0 and len(list_X_cm_t)<Len):
      list_X_cm_t.append(CM(R)[0]-R[0,0]) # Append new position wrt CM and thermalized
      list_Y_cm_t.append(CM(R)[1]-R[1,0])
      list_X_cm_t_2.append(CM(R)[0]-R[0,1])  # Second particle's coordinates wrt CM and thermalized
      list_Y_cm_t_2.append(CM(R)[1]-R[1,1]) 
      list_X_cm_t_3.append(CM(R)[0]-R[0,2]) # Third particle's coordinates wrt CM and thermalized
      list_Y_cm_t_3.append(CM(R)[1]-R[1,2])
      list_X_cm_t_4.append(CM(R)[0]-R[0,3]) # Fouth particle's coordinates wrt CM and thermalized
      list_Y_cm_t_4.append(CM(R)[1]-R[1,3])
      list_X_t.append(R[0,0]) # Respect to an external coordinates system
      list_Y_t.append(R[1,0])
      
  else:
    R = R.copy()

  E_loc = E_loc + V(R) - (K(R)/WF(R, gamma, a, beta))
  SE = SE + E_loc**2
  
print("Mean Energy = ", E_loc/NM, "+/-", math.sqrt((SE/NM)-(E_loc/NM)**2)/NM, "with acceptance ratio", NA/NM*100,"%\n", "With the following paramters:\ngamma =", gamma, ", a =", a, ", beta =", beta)


fig, ax = plt.subplots(figsize=(5, 5), dpi=600)
ax.plot(list_X_t, list_Y_t, color = "black")
ax.axvline(x=0,color='black')
ax.axhline(y=0, color='black')
ax.set_xlabel('x (fm)')
ax.set_ylabel('y (fm)')
plt.title("Particle's position in 2D wrt external source")
plt.show()

fig, ax = plt.subplots(figsize=(5, 5), dpi=600)
ax.scatter(list_X_cm, list_Y_cm, s= 2,c = "black"  )
ax.set_xlim(-3, 3) 
ax.set_ylim(-3, 3) 
ax.vlines(0,-3,3,colors='black')
ax.hlines(0,-3,3,colors='black')
ax.set_xlabel('x (fm)')
ax.set_ylabel('y (fm)')
ax.set_title("Particle's position in 2D wrt cm")
plt.show()

plt.figure(figsize=(5, 5), dpi=600)
plt.scatter(list_X_cm_t, list_Y_cm_t, s= 2,c = "black"  )
plt.xlim(-3, 3)  # Set x-axis limits
plt.ylim(-3, 3)  # Set y-axis limits
plt.vlines(0,-3,3,colors='black')
plt.hlines(0,-3,3,colors='black')
plt.xlabel('x (fm)')
plt.ylabel('y (fm)')
plt.title("Particle's position in 2D wrt cm with thermalization")
plt.show()

fig, ax = plt.subplots(figsize=(5, 5), dpi=600)
ax.plot(list_X_cm_t[0:25], list_Y_cm_t[0:25], color = "blue", label= "First particles' positions")
ax.axvline(x=0,color='black')
ax.axhline(y=0, color='black')
ax.set_xlabel('x (fm)')
ax.set_ylabel('y (fm)')
ax.plot(list_X_cm_t_2[0:25], list_Y_cm_t_2[0:25], c="purple", label = "Second particles' positions")
ax.plot(list_X_cm_t_3[0:25], list_Y_cm_t_3[0:25], c="cyan", label = "Third particles' positions")
ax.plot(list_X_cm_t_4[0:25], list_Y_cm_t_4[0:25], c="violet", label = "Fourth particles' positions")
ax.set_title("Particles' positions in 2D wrt cm")
ax.legend(fontsize = 'x-small')
plt.show()

fig = plt.figure(figsize=(6,6), dpi=100)
ax = fig.add_subplot(projection='3d')
ax.scatter(list_X_cm, list_Y_cm, list_Z_cm, s=2, c="black")
plt.xlabel('x (fm)')
plt.ylabel('y (fm)')
ax.set_zlabel('z (fm)')
plt.title("Particle's position in 3D wrt cm")
plt.show()

