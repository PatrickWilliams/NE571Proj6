#####################################################################
# Author: Zack Taylor
# Class: NE 571
# Assignment: Project 6
#####################################################################

# Libraries
import numpy as np
import scipy.special as sp
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#####################################################################
# Material class for calling material properties
#####################################################################
class Material(object):
	# Takes in a list of values for each fuel element. The list 
	# contains two values. The value of interest for each group
	# Example: D[(fast value),(thermal value)]. The order is 
	# fast then thermal. The element attribute is something I added
	# to help test and make sure GetMaterial was return the right
	# values. This can be removed. 
    def __init__(self,sigTr,siga,vsig,sigS,D,diameter,delR):
        self.sigTr_ = sigTr
        self.siga_ = siga
        self.vsig_ = vsig
        self.sigS_ = sigS
        self.D_ = D
        self.diameter = diameter
        self.delR = delR

    def __str__(self):
        return self.delR
    
    def GetD(self,i,gp):
        if gp == "fast":
            D = self.D_[0]
        if gp == "thermal":
            D = self.D_[1]
        return D

    def GetsigTr(self,i,gp):
        sigTr = self.sigTr_[0]
        return sigTr

    def GetsigS(self,i,gp):
        sigS = self.sigS_[0]
        return sigS

    def Getsiga(self,i,gp):
        siga = self.siga_[0]
        return siga

    def Getvsig(self,i,gp):
        if gp == "fast":
            vsig = self.vsig_[0]
        elif gp == "thermal":
            vsig = self.vsig_[1]
        return vsig


#####################################################################
# Functions used in the progrm 
#####################################################################

# Input: the node location in the i direction (r direction). Index 
# at 0. 
# Output: The distance from the center at that node. AKA the length. 
def r(i):
    x = i*delR_1
    #x2 = i*delR_2
    if x < R_1:
        return x
    elif x >= R_1:
        return R_1 + (i-Ni+1)*delR_1

# Function that returns the material or materials given node i
# Uses the r function to return a system diameter value. 
def GetMaterial(i):
    r_ = [] # holdes values for the lenght of each fuel element
    temp = 0.0
    InterfaceValues = []
    NumberOfFuelElements = len(LoadingPattern)

    # Sets up the arrays to be used. r_ holds the length of each
    # fuel element. InterfaceValues holds the values for all 
    # interfaces in the system
    for k in xrange(NumberOfFuelElements-1):
        temp = temp + LoadingPattern[k].diameter
        r_.append(LoadingPattern[k].diameter)
        InterfaceValues.append(temp)

    # Returns the what fuel element you are in for a 
    # non interface value
    for n,k in enumerate(InterfaceValues):
        if r(i) < r_[0]:
            return LoadingPattern[0]
        if InterfaceValues[n-1] < r(i) < InterfaceValues[n]:
            return LoadingPattern[n]
        if r(i) > max(InterfaceValues):
            return Water

    # Returns two values at an interface node. The first is the 
    # material on the left, second is the material on the right. 
    for n,k in enumerate(InterfaceValues):
    	# Added to fix if the user only addes one fuel element
        if r(i) == InterfaceValues[0] and len(LoadingPattern) == 1:
            return LoadingPattern[0],Water
        if r(i) == InterfaceValues[0]:
            return LoadingPattern[0], LoadingPattern[1]
        if r(i) == InterfaceValues[n] and n != (len(InterfaceValues)-1):
            return LoadingPattern[n],LoadingPattern[n+1]
        if r(i) == max(InterfaceValues):
            return LoadingPattern[(len(InterfaceValues)-1)],Water

# Input: Node in the i direction. 
# Output: True or false logical to weither the node is at an 
# interface or not. 
def IsInterface(i):
    InterfaceValues = []
    Interface = False
    temp = 0.0
    NumberOfFuelElements = len(LoadingPattern)-1
    for k in xrange(NumberOfFuelElements):
        temp = temp + LoadingPattern[k].diameter
        InterfaceValues.append(temp)

    for k in InterfaceValues:
        if r(i) == k:
            Interface = True

    return Interface

# Calculates the total length of the fuel region. 
def FuelRegionLength(LoadingPattern):
    temp = 0.0
    NumberOfFuelElements = len(LoadingPattern)

    for k in xrange(NumberOfFuelElements):
        temp = temp + LoadingPattern[k].diameter
    return temp 

# Input is the number of mesh cells for each fuel element!!!
# returns the correct number of nodes. 
def MeshGeneration(NumberOfCellsFuelElementdiameter,NumberOfCellsWaterdiameter):
    cE = NumberOfCellsFuelElementdiameter
    cW = NumberOfCellsWaterdiameter

    Ni = len(LoadingPattern)*cE + 1 #Number of nodes in fuel region. Includes 
    # fuel water interface
    Nii = cW #Does not include fuel water interface

    Nitot = Ni + Nii

    return Ni,Nii,Nitot

# All the functions used in the A matrix generation. 
def a1(i,gp):
    if IsInterface(i) == False:
        material = GetMaterial(i)
        if r(i) == 0:
            return -material.GetD(i,gp)*(delR_1**3)/8
        elif (0 < r(i) < R_1 or r(i) > R_1):
            return -material.GetD(i,gp)*r(i)*delR_1**2

    elif IsInterface(i) == True:
        Leftmaterial,Rightmaterial = GetMaterial(i)
        D1 = Leftmaterial.GetD(i,gp)
        D2 = Rightmaterial.GetD(i,gp)
        return (-D1/delZ)*((r(i)*delR_1/2) + (delR_1**2/8)) + (-D2/delZ)*((r(i)*delR_1/2) + (delR_1**2/8))
    
def a2(i,gp):
    if IsInterface(i) == False:
        material = GetMaterial(i)
        if r(i) == 0:
            return -material.GetD(i,gp)*(delR_1**3)/8
        if 0 < r(i) < R_1 or r(i) > R_1:
            return -material.GetD(i,gp)*(r(i)-(delR_1/2))*delZ**2

    elif IsInterface(i) == True:
        Leftmaterial,Rightmaterial = GetMaterial(i)
        D1 = Leftmaterial.GetD(i,gp)
        D2 = Rightmaterial.GetD(i,gp)
        return (-D1/delR_1)*(r(i)-delR_1/2)*delZ
    
def a3(i,gp):
    if IsInterface(i) == False:
        material = GetMaterial(i)
        if gp == "fast":
            GG = material.GetsigTr(i,gp)
        else:
            GG = material.Getsiga(i,gp)
        if r(i) == 0:
            return GG*((delR_1**3)/8)*(delZ**2) + 2*material.GetD(i,gp)*(delR_1**3)/8 + material.GetD(i,gp)*(delR_1/2)*delZ**2
        elif 0 < r(i) < R_1 or r(i) > R_1:
            return GG*r(i)*delR_1**2*delZ**2 + 2*material.GetD(i,gp)*r(i)*delR_1**2 + 2*material.GetD(i,gp)*r(i)*delZ**2 
    elif IsInterface(i) == True:
        Leftmaterial,Rightmaterial = GetMaterial(i)
        if gp == "fast":
            GG1 = Rightmaterial.Getsiga(i,gp)
            GG2 = Leftmaterial.Getsiga(i,gp)
        else:
            GG1 = Rightmaterial.GetsigTr(i,gp)
            GG2 = Leftmaterial.GetsigTr(i,gp)
        D1 = Leftmaterial.GetD(i,gp)
        D2 = Rightmaterial.GetD(i,gp)
        return (D1/delR_1)*(r(i)-delR_1/2)*delZ + (2*D1/delZ)*((r(i)*delR_1/2) + (delR_1**2/8)) + GG1*delZ*((r(i)*delR_1/2) + (delR_1**2/8)) + (D2/delR_1)*(r(i)-delR_1/2)*delZ + (2*D2/delZ)*((r(i)*delR_1/2) + (delR_1**2/8)) + GG2*delZ*((r(i)*delR_1/2) + (delR_1**2/8))

def a4(i,gp):
    if IsInterface(i) == False:
        material = GetMaterial(i)
        if r(i) == 0:
            return -material.GetD(i,gp)*(delR_1/2)*delZ**2
        elif (0 < r(i) < R_1 or r(i) > R_1):
            return -material.GetD(i,gp)*(r(i) + delR_1/2)*delZ**2
    elif IsInterface(i) == True:
        Leftmaterial,Rightmaterial = GetMaterial(i)
        D1 = Leftmaterial.GetD(i,gp)
        D2 = Rightmaterial.GetD(i,gp)
        return (-D2/delR_1)*(r(i)-delR_1/2)*delZ

def a5(i,gp):
    if IsInterface(i) == False:
        material = GetMaterial(i)
        if r(i) == 0:
            return -material.GetD(i,gp)*(delR_1**3)/8
        elif 0 < r(i) < R_1 or r(i) > R_1:
            return -material.GetD(i,gp)*r(i)*delR_1**2
    elif IsInterface(i) == True:
        Leftmaterial,Rightmaterial = GetMaterial(i)
        D1 = Leftmaterial.GetD(i,gp)
        D2 = Rightmaterial.GetD(i,gp)
        return (-D1/delZ)*((r(i)*delR_1/2) + (delR_1**2/8)) + (-D2/delZ)*((r(i)*delR_1/2) + (delR_1**2/8))
      
def b_fis_term(i,gp):
    if IsInterface(i) == False:
        material = GetMaterial(i)
        if r(i) == 0:
            return material.Getvsig(i,gp)*((delR_1**3)/8)*delZ**2
        elif r(i) < R_1:
            return material.Getvsig(i,gp)*delZ**2*r(i)*delR_1**2
        else:
            return 0.0
    elif IsInterface(i) == True:
        print i,r(i)
        Leftmaterial,Rightmaterial = GetMaterial(i)
        vsig1 = Leftmaterial.Getvsig(i,gp)
        vsig2 = Rightmaterial.Getvsig(i,gp)
        return vsig1*delZ**2*r(i)*delR_1**2 
    
def b_down_term(i,gp):
    if IsInterface(i) == False:
        material = GetMaterial(i)
        if r(i) == 0:
            return material.GetsigS(i,gp)*((delR_1**3)/8)*delZ**2
        else: 
            return material.GetsigS(i,gp)*delZ**2*r(i)*delR_1**2
    elif IsInterface(i) == True:
        Leftmaterial,Rightmaterial = GetMaterial(i)
        sigTr1 = Leftmaterial.GetsigS(i,gp)
        sigTr2 = Rightmaterial.GetsigS(i,gp)
        return sigTr1*delZ**2*r(i)*delR_1**2 + sigTr2*delZ**2*r(i)*delR_1**2

# Takes in the A matrix and fills in its values. Note its very messy 
# but i think it works. 
def a_matrix_gen(A,i,k,gp):
	# Loops through the bottom nodes in the problem. Sweeps from 
	# center of reactor out. 
    for h in xrange (0,i-1):
        if h == 0:
            A[h][0] = a3(h,gp)
            A[h][1] = a4(h,gp)
            A[h][2+i0] = a5(h,gp)
        if h == i-2:
            A[h][i0] = a2(h,gp)
            A[h][i0+1] = a3(h,gp)
            A[h][2*i0+3] = a5(h,gp)
        if 0 < h < i-2:
            A[h][h-1] = a2(h,gp)
            A[h][h] = a3(h,gp)
            A[h][h+1] = a4(h,gp)
            A[h][h+2+i0] = a5(h,gp)
            
# Generates part 2 of the a matrix
# These are the rows that are sandwitched between the bottom row
# and the top row. 
    oldu = 0
    for u in xrange(0,k-3):
        for h in xrange(u*(i-1),u*(i-1) + (i-1)):

            if u !=0:
                if h == u*i-oldu-1:       
                    A[h][h-i0-2] = a1(h-(u)*(i-1),gp)
                    A[h][h] = a3(h-(u)*(i-1),gp)
                    A[h][h+1] = a4(h-(u)*(i-1),gp)
                    A[h][h+2+i0] =a5(h-(u)*(i-1),gp)
                if h == u*(i-1) + (i-2): #and u*Ni+Ni-1 != Ni-1:
                    A[h][h-i0-2] = a1(h-(u)*(i-1),gp)
                    A[h][h] = a3(h-(u)*(i-1),gp)
                    A[h][h-1] = a2(h-(u)*(i-1),gp)
                    A[h][h+i0+2] = a5(h-(u)*(i-1),gp)
            
                
                if u*(i-1) < h < u*(i-1) + (i-2):
                    A[h][h-i0-2] = a1(h-(u)*(i-1),gp)
                    A[h][h-1] = a2(h-(u)*(i-1),gp)
                    A[h][h] = a3(h-(u)*(i-1),gp)
                    A[h][h+1] = a4(h-(u)*(i-1),gp)
                    A[h][h+2+i0] = a5(h-(u)*(i-1),gp)
                
        oldu = u

# Generates part 3 of the matrix. The top row of values at Z = totalZ
    for h in xrange((i-1)*(k-2)-i+1,(i-1)*(k-2)):
        if h == (i-1)*(k-2)-i+1:
            A[h][h] = a3(h-(u+1)*(i-1),gp)
            A[h][h-i0-2] = a1(h-(u+1)*(i-1),gp)
            A[h][h+1] = a4(h-(u+1)*(i-1),gp)
        
        if h == (i-1)*(k-2)-1:
            A[h][h-2-i0] = a1(h-(u+1)*(i-1),gp)
            A[h][h-1] = a2(h-(u+1)*(i-1),gp)
            A[h][h] = a3(h-(u+1)*(i-1),gp)
        
        if (i-1)*(k-2)-i+1 < h < (i-1)*(k-2)-1:
            A[h][h+1] = a4(h-(u+1)*(i-1),gp)
            A[h][h] = a3(h-(u+1)*(i-1),gp)
            A[h][h-1] = a2(h-(u+1)*(i-1),gp)
            A[h][h-2-i0] = a1(h-(u+1)*(i-1),gp)
            
# Same thing as A matrix but for B. 
def b_matrix_gen(B,b,gp):
    for u in xrange(0,Nk-3):
        for h in xrange(0,(Nitot-1)*(Nk-2)):
            if 0 <= h < Nitot-1:
                B[h][h] = b(h,gp)
            if u*(Nitot-1) <= h <= u*(Nitot-1) + (Nitot-2):
                B[h][h] = b((h-(u)*(Nitot-1)),gp)
            if (Nitot-1)*(Nk-2)-Nitot+1 <= h < (Nitot-1)*(Nk-2):
                B[h][h] = b(h-(u+1)*(Nitot-1),gp)

#####################################################################
# Main program 
#####################################################################


# System properties. I put dummy variables in here.
# k = 1.0873
AsigTr_ = [0.0269]
Asiga_ = [0.1152]
Avsig_ = [0.00641,0.169]
AsigS_ = [0.015561] # 0.00155
AD_ = [1.2627,0.3543]

# k = 0.977
BsigTr_ = [0.0274]
Bsiga_ = [0.108]
Bvsig_ = [0.0055,0.146]
BsigS_ = [0.0157]
BD_ = [1.2427,0.3543]

# k = 1.21
CsigTr_ = [0.0263]
Csiga_ = [0.119]
Cvsig_ = [0.00746,0.186]
CsigS_ = [0.0156]
CD_ = [1.2627,0.3543]

# k = 1.44
DsigTr_ = [0.0254]
Dsiga_ = [0.1097]
Dvsig_ = [0.00873,0.198]
DsigS_ = [0.0155]
DD_ = [1.2627,0.3543]

WsigTr_ = [0.0494]
Wsiga_ = [0.0197]
Wvsig_ = [0.0,0.0]
WsigS_ = [0.0494]
WD_ = [1.13,0.16]


# Loading patter for fuel region. Note that this does include water
# if you want to only do one fuel element for the whole core
# LoadingPattern = [FuelX]. MeshGeneration(X,0). delR_2 = 0.0

FuelA = Material(AsigTr_,Asiga_,Avsig_,AsigS_,AD_,18.0,"FuelA")
FuelB = Material(BsigTr_,Bsiga_,Bvsig_,BsigS_,BD_,18.0,"FuelB")
FuelC = Material(CsigTr_,Csiga_,Cvsig_,CsigS_,CD_,18.0,"FuelC")
FuelD = Material(DsigTr_,Dsiga_,Dvsig_,DsigS_,DD_,18.0,"FuelD")
Water = Material(WsigTr_,Wsiga_,Wvsig_,WsigS_,WD_,15.0,"Water")

#LoadingPattern = [FuelD,FuelD,FuelD,FuelD,FuelD,FuelD,FuelD,FuelD,FuelD,FuelD,FuelD,FuelD,FuelB,FuelB,FuelB,FuelB,FuelB,Water,Water,Water,Water]
LoadingPattern = [FuelB,FuelA,FuelB,FuelC,Water]

Z = 2000.0  # height of the reactor
R_1 = FuelRegionLength  # radious of the reactor
#R_2 = 5.0 # radious of the reflector
Ni,Nii,Nitot = MeshGeneration(4,0)
Nk = 10
i0 = Nitot - 3 # This is the number of Zeros between the main
			   # triple diagional and the offset diagional

num_row_r1 = (Nitot-1)*(Nk-2)
num_col_r1 = (Nitot-1)*(Nk-2)

delR_1 = FuelRegionLength(LoadingPattern)/(Nitot)

delZ = Z/Nk
for i in xrange(Nitot):
    print i,r(i)
    print GetMaterial(i)
# A matrix
A_f = np.zeros(shape=(num_row_r1,num_col_r1))
A_t = np.zeros(shape=(num_row_r1,num_col_r1))

a_matrix_gen(A_f,Nitot,Nk,"fast")
a_matrix_gen(A_t,Nitot,Nk,"thermal")

# B matrix
B_fis_f = np.zeros(shape=(num_row_r1,num_col_r1))
B_fis_t = np.zeros(shape=(num_row_r1,num_col_r1))
B_down_t = np.zeros(shape=(num_row_r1,num_col_r1))

# Populates the B matrix
b_matrix_gen(B_fis_f,b_fis_term,'fast')
b_matrix_gen(B_fis_t,b_fis_term,'thermal')
b_matrix_gen(B_down_t,b_down_term,'thermal')


# Solves the actual problem

# Initial guess for flux
fluxT = np.ones(shape=(num_row_r1,1))
fluxF = np.ones(shape=(num_row_r1,1))
    
# Source term
k = 1
ST = (1/k)*(B_down_t.dot(fluxF))
SF = (1/k)*(B_fis_t.dot(fluxT) + B_fis_f.dot(fluxF))

# Iterations
fluxdiff = 1
j = 0
while(fluxdiff > 0.001):
    
    oldfluxF = fluxF
    oldfluxT = fluxT
    
    oldk = k
    oldSF = SF
    oldST = ST
    
    fluxF = np.linalg.inv(A_f).dot(oldSF)
    fluxT = np.linalg.inv(A_t).dot(oldST)
    SF = (1/k)*(B_fis_t.dot(fluxT) + B_fis_f.dot(fluxF))
    ST = (1/k)*(B_down_t.dot(fluxF))
    k = oldk*(np.sum(SF)/np.sum(oldSF))
    
    fluxdiff = 0
    
    for i in (0,Ni):
        fluxdiff = fluxdiff + ((fluxF[i][0] - (oldfluxF[i][0]))**2)
        fluxdiff = fluxdiff**0.5
    print fluxdiff
    j = j + 1


# Normalize
total = np.sum(fluxF)
fluxF = fluxF/total
print j
print k

fig = plt.figure()
ax = fig.gca(projection='3d')

xDatapoints = []
for i in xrange(Nitot-1):
    xDatapoints.append(r(i))
x = xDatapoints
y = np.arange(delZ,Z-delZ,delZ)
xx, yy = np.meshgrid(x,y)
#fig = plt.figure()
z = np.array((fluxF.reshape(xx.shape)))
#ax = fig.add_subplot(211, projection='3d')
ax.plot_surface(xx,yy,z,cmap=plt.cm.rainbow,linewidth=0, antialiased=False)
ax.set_xlabel('Radius (cm)')
ax.set_ylabel('Height (cm)')
ax.set_zlabel('$\phi (normalized)$')
#plt.draw()
plt.show()

"""
bx = fig.add_subplot(212)
bx.pcolor(x,y,fluxT.reshape(Nk-2,Nitot-1),cmap=plt.cm.rainbow)
bx.set_xlabel('Radius (cm)')
bx.set_ylabel('Height (cm)')
plt.draw()
plt.show()

fig = plt.figure()
ax = fig.gca(projection='3d')

"""

