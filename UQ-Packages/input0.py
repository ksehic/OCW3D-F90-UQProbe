# coding: utf-8
# UQpackage by kense, Aug 2018
# ###################
import numpy as np

# Inputs for the wave generation at x=0

TDur=600 # time duration [s]

dt=0.01 # time step [s]

scl = 50 #scaling factor

Hs=6.8 #significant wave height

HsL=Hs/float(scl) # significant wave height [m]

Tp=15 #time period

TpL=Tp/float(np.sqrt(scl)) # time period [s]

df=1/float(TDur) # frequ step

f0 = df # freq start

fHighCut=0.3 # freq cut

y = 0.5 #

ny = 1 # steps in the y-axis

################# Lets generate Beach space domain
import matlab.engine # to call a matlab function

eng = matlab.engine.start_matlab()

L0 = 15 # Flat bottom run-up length at the left end L_0
L1 =  5 # Sloping beach length L_1=?
h0 =  5 # Depth at x0
h1=   3 # Depth at x1
n1 =  20 # Number of grid points on the sloping beach part n_1?

eng.BeachGen(L0,L1,n1,h0,h1) # Lets run Matlab

################# Lets define input parameters for OceanWave3D.inp


################ Lets generate Wave Maker script
import wavegen

wavegen.x0(TDur,dt,HsL,TpL,f0,df,fHighCut,y,ny) #produce Wave Maker script

########### Create a folders and make a copy within
import os
import shutil

run = 1

loc = os.getcwd()

path = os.getcwd() + "/run%d" % (run) #run will be each new random realisation

runpath= os.mkdir(path)

shutil.copy(loc+"beach_%d" %(run), runpath)
shutil.copy(loc+"OceanWave3D.inp", runpath)
shutil.copy(loc+"wavegen.py", runpath)

######
