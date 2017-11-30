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

    def __init__(self,sigTr,siga,vsig,D,radius,element):
        self.sigTr_ = sigTr
        self.siga_ = siga
        self.vsig_ = vsig
        self.D_ = D
        self.radius = radius
        self.element = element

    def __str__(self):
        return self.element

    def GetD(self,i,gp):
        if gp == "fast":
            D = self.D_[0]
        if gp == "thermal":
            D = self.D_[1]
        return D

    def GetsigTr(self,i,gp):
        if gp == "fast":
            sigTr = self.sigTr_[0]
        if gp == "thermal":
            sigTr = self.sigTr_[1]
        return sigTr

    def Getsiga(self,i,gp):
        if gp == "fast":
            siga = self.siga_[0]
        if gp == "thermal":
            siga = self.siga_[1]
        return siga

    def Getvsig(self,i,gp):
        if gp == "fast":
            vsig = self.vsig_[0]
        elif gp == "thermal":
            vsig = self.vsig_[1]


#####################################################################
# Functions used in the progrm 
#####################################################################
def r(i):
    x = i*delR_1
    x2 = i*delR_2
    if x < R_1:
        return x
    elif x >= R_1:
        return R_1 + (i-Ni+1)*delR_2

# Function that returns the material or materials given node i
# Uses the r function to return a system radius value. Heart of
# the program
def GetMaterial(i):
    r_ = []
    temp = 0.0
    InterfaceValues = []
    NumberOfFuelElements = len(LoadingPattern)

    for k in xrange(NumberOfFuelElements):
        temp = temp + LoadingPattern[k].radius
        r_.append(LoadingPattern[k].radius)
        InterfaceValues.append(temp)

    for n,k in enumerate(InterfaceValues):
        if r(i) < r_[0]:
            return LoadingPattern[0]
        if InterfaceValues[n-1] < r(i) < InterfaceValues[n]:
            return LoadingPattern[n]
        if r(i) > max(InterfaceValues):
            return Water

    for n,k in enumerate(InterfaceValues):
        if r(i) == InterfaceValues[0]:
            return LoadingPattern[0], LoadingPattern[1]
        if r(i) == InterfaceValues[n] and n != (len(InterfaceValues)-1):
            return LoadingPattern[n],LoadingPattern[n+1]
        if r(i) == max(InterfaceValues):
            return LoadingPattern[(len(InterfaceValues)-1)],Water

def IsInterface(i):
    InterfaceValues = []
    Interface = False
    temp = 0.0
    NumberOfFuelElements = len(LoadingPattern)
    for k in xrange(NumberOfFuelElements):
        temp = temp + LoadingPattern[k].radius
        InterfaceValues.append(temp)

    for k in InterfaceValues:
        if r(i) == k:
            Interface = True

    return Interface

def FuelRegionLength(LoadingPattern):
    temp = 0.0
    NumberOfFuelElements = len(LoadingPattern)

    for k in xrange(NumberOfFuelElements):
        temp = temp + LoadingPattern[k].radius
    return temp 

def MeshGeneration(NumberOfCellsFuelElementRadius,NumberOfCellsWaterRadius):
    cE = NumberOfCellsFuelElementRadius
    cW = NumberOfCellsWaterRadius

    Ni = len(LoadingPattern)*cE + 1 #Number of nodes in fuel region. Includes 
    # fuel water interface
    Nii = cW #Does not include fuel water interface

    Nitot = Ni + Nii

    return Ni,Nii,Nitot

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
        return (-D1/delZ)*((r(i)*delR_1/2) + (delR_1**2/8)) + (-D2/delZ)*((r(i)*delR_2/2) + (delR_2**2/8))
    
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
        if r(i) == 0:
            return material.Getsiga(i,gp)*((delR_1**3)/8)*(delZ**2) + 2*material.GetD(i,gp)*(delR_1**3)/8 + material.GetD(i,gp)*(delR_1/2)*delZ**2
        elif 0 < r(i) < R_1 or r(i) > R_1:
            return material.Getsiga(i,gp)*r(i)*delR_1**2*delZ**2 + 2*material.GetD(i,gp)*r(i)*delR_1**2 + 2*material.GetD(i,gp)*r(i)*delZ**2 
    elif IsInterface(i) == True:
        Leftmaterial,Rightmaterial = GetMaterial(i)
        D1 = Leftmaterial.GetD(i,gp)
        D2 = Rightmaterial.GetD(i,gp)
        siga1 = Leftmaterial.Getsiga(i,gp)
        siga2 = Rightmaterial.Getsiga(i,gp)
        sigTr1 = Leftmaterial.GetsigTr(i,gp)
        sigTr2 = Rightmaterial.GetsigTr(i,gp)
        return (D1/delR_1)*(r(i)-delR_1/2)*delZ + (2*D1/delZ)*((r(i)*delR_1/2) + (delR_1**2/8)) + siga1*delZ*((r(i)*delR_1/2) + (delR_1**2/8)) + (D2/delR_2)*(r(i)-delR_2/2)*delZ + (2*D2/delZ)*((r(i)*delR_2/2) + (delR_2**2/8)) + siga2*delZ*((r(i)*delR_2/2) + (delR_2**2/8))

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
        siga1 = Leftmaterial.Getsiga(i,gp)
        siga2 = Rightmaterial.Getsiga(i,gp)
        sigTr1 = Leftmaterial.GetsigTr(i,gp)
        sigTr2 = Rightmaterial.GetsigTr(i,gp)
        return (-D1/delZ)*((r(i)*delR_1/2) + (delR_1**2/8)) + (-D2/delZ)*((r(i)*delR_2/2) + (delR_2**2/8))
      
def b_fis_therm(i):
    material = GetMaterial(i)
    if r(i) == 0:
        return material.Getvsig(i,gp)*((delR_1**3)/8)*delZ**2
    elif r(i) <= R_1:
        return material.Getvsig(i,gp)*delZ**2*r(i)*delR_1**2
    else:
        return 0.0
    
def b_source_therm(i):
    if r(i) == 0:
        return sigTrR_1*((delR_1**3)/8)*delZ**2
    elif r(i) <= R_1:
        return sigTrR_1*delZ**2*r(i)*delR_1**2
    elif (r(i) > R_1):
        return sigTrR_2*delZ**2*r(i)*delR_2**2

def a_matrix_gen(A,i,k,rg,gp):
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

# Generates part 3 of the matrix
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
            
            
def b_matrix_gen(B,b):
    # Populates the B matrix
    for u in xrange(0,Nk-3):
        for h in xrange(0,(Ni-1)*(Nk-2)):
            if 0 <= h < Ni-1:
                B[h][h] = b(h)
            if u*(Ni-1) <= h <= u*(Ni-1) + (Ni-2):
                B[h][h] = b((h-(u)*(Ni-1)))
            if (Ni-1)*(Nk-2)-Ni+1 <= h < (Ni-1)*(Nk-2):
                B[h][h] = b(h-(u+1)*(Ni-1))

#####################################################################
# Main program 
#####################################################################

AsigTr_ = [3.62e-2,0.0494]
Asiga_ = [0.01207,0.0004]
Avsig_ = [0.008476,0.18514]
AD_ = [1.2627,1.13]
rad = 20.0

BsigTr_ = [3.62e-2,0.0494]
Bsiga_ = [0.01207,0.0004]
Bvsig_ = [0.008476,0.18514]
BD_ = [1.2627,1.13]
rad = 20.0

CsigTr_ = [3.62e-2,0.0494]
Csiga_ = [0.01207,0.0004]
Cvsig_ = [0.008476,0.18514]
CD_ = [1.2627,1.13]
rad = 20.0

FuelA = Material(AsigTr_,Asiga_,Avsig_,AD_,rad,"FuelA")
FuelB = Material(BsigTr_,Bsiga_,Bvsig_,BD_,rad,"FuelB")
FuelC = Material(CsigTr_,Csiga_,Cvsig_,CD_,rad,"FuelC")
Water = Material(CsigTr_,Csiga_,Cvsig_,CD_,rad,"Water")

LoadingPattern = [FuelA,Water]

Z = 100.0  # height of the reactor
R_1 = 60.0  # radious of the reactor
R_2 = 20.0 # radious of the reflector
Ni,Nii,Nitot = MeshGeneration(4,0)
Nk = 6
i0 = Nitot - 2

num_row_r1 = (Nitot-1)*(Nk-2)
num_col_r1 = (Nitot-1)*(Nk-2)

delR_1 = FuelRegionLength(LoadingPattern)/(Ni-1)
#delR_2 = R_2/(Nii)
delZ = Z/Nk


# A matrix
A_f = np.zeros(shape=(num_row_r1,num_col_r1))
A_t = np.zeros(shape=(num_row_r1,num_col_r1))

a_matrix_gen(A_f,Nitot,Nk,1,"fast")
a_matrix_gen(A_t,Nitot,Nk,1,"thermal")
"""

# B matrix
B_fis_f = np.zeros(shape=(num_row_r1,num_col_r1))
B_fis_t = np.zeros(shape=(num_row_r1,num_col_r1))
B_source_t = np.zeros(shape=(num_row_r1,num_col_r1))
# Populates the B matrix

b_matrix_gen(B_fis_f,b_fis_fast)
b_matrix_gen(B_fis_t,b_fis_therm)
b_matrix_gen(B_source_t,b_source_therm)


# Solves the actual problem

# Initial guess for flux
fluxT = np.ones(shape=(num_row_r1,1))
fluxF = np.ones(shape=(num_row_r1,1))
    
# Source term
k = 1
ST = (1/k)*(B_source_t.dot(fluxF))
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
    ST = (1/k)*(B_source_t.dot(fluxF))
    k = oldk*(np.sum(SF)/np.sum(oldSF))
    
    fluxdiff = 0
    
    for i in (0,Ni):
        fluxdiff = fluxdiff + ((fluxF[i][0] - (oldfluxF[i][0]))**2)
        fluxdiff = fluxdiff**0.5
    j = j + 1
    print fluxdiff


# Normalize
total = np.sum(fluxF)
fluxF = fluxF/total

print j
print k

x = np.arange(0,R_1 + R_2-delR_1-delR_2,delR_1+delR_2)
print x
y = np.arange(delZ,Z-delZ,delZ)
xx, yy = np.meshgrid(x,y)
z = np.array((flux.reshape(xx.shape)))
fig = plt.figure()
ax = fig.add_subplot(211, projection='3d')
ax.plot_surface(xx,yy,z,cmap=plt.cm.coolwarm)
ax.set_xlabel('Radius (cm)')
ax.set_ylabel('Height (cm)')
ax.set_zlabel('$\phi (normalized)$')
plt.draw()
plt.show()

bx = fig.add_subplot(212)
bx.pcolor(x,y,flux.reshape(Nk-2,Ni-1),cmap=plt.cm.rainbow)
bx.set_xlabel('Radius (cm)')
bx.set_ylabel('Height (cm)')
plt.draw()
plt.show()

fig = plt.figure()
ax = fig.gca(projection='3d')

"""


