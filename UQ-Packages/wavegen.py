# coding: utf-8
# WaveMaker Code
# by kense
# ############
import numpy as np
import matplotlib.pyplot as plt

np.random.seed(seed=None) # random seed

def jonswap(f,Hs,Tp): # JONSWAP function

	fp=1/float(Tp)

	TH = Tp/float(np.sqrt(Hs))

	if TH<=3.6:

			gamma=5

	elif (3.6<TH) and (TH<=5):

			gamma=np.exp(5.75-1.15*TH)

	else:
			gamma=1

	if f<=fp:

    		sigma = 0.07

	else:

    		sigma = 0.09

	s=0.3125*Hs**2*Tp*(f/float(fp))**(-5)*np.exp(-1.25*(f/float(fp))**(-4))*(1-0.287*np.log(gamma))*gamma**(np.exp(-0.5*((f/float(fp)-1)/float(sigma))**2))

	return s

def x0(TDur,dt,Hs,Tp,f0,df,fHighCut,y,ny): # etaBC generation

	f=np.arange(f0,fHighCut+f0,df)

	omega=2*np.pi*f

	JS = np.zeros((1,f.size))

	for ij in range(0,f.size-1):

		fL=f[ij]

		JS1 = jonswap(fL,Hs,Tp)

		JS[0,ij] = JS1

	#tSpan=0:dt:TDur;
	tSpan=np.arange(0,TDur+dt,dt)

	etaBC=np.zeros((tSpan.size,1))

	A=np.zeros((1,omega.size))
	B=np.zeros((1,omega.size))
	AB=np.random.standard_normal((1,2*omega.size)) #random Fourie coefficients

	for j in range(0,omega.size):

	    A[0,j]=np.sqrt(JS[0,j]*df)*AB[0,j]
	    B[0,j]=np.sqrt(JS[0,j]*df)*AB[0,omega.size+j]

	for i in range(0,tSpan.size):

		etaBC[i,0]=np.sum(A*np.cos(2*np.pi*f*tSpan[i])+B*np.sin(2*np.pi*f*tSpan[i])); # Initial Wave Surface

	# Print Check
	#print omega.size

	#print etaBC
	#plt.plot(tSpan,etaBC,'b-*')
	#plt.xlabel('Time[s]')
	#plt.ylabel('Eta[m]')
	#plt.grid()
	#plt.show()

	n = tSpan.size

	# Lets print output the Wave Maker

	file = open("waveMakerSignalrandom.inp","w")

	file.write("# This is a Wave Generation surface response at x=0\n")
	file.write("%.7f %.0f %.0f\n" % (dt, n, ny))
	file.write("%.7f\n" % y)

	for i in range(0,tSpan.size):
	  file.write("%.7f %.7f\n" % (tSpan[i],etaBC[i,0]))

	file.close()

	## Which Fourie coefficients have been used

	file = open("randomAB.txt","w")

	file.write("# These are random Aj and Bj\n")

	for i in range(0,omega.size):
	  file.write("%.7s %.7s\n" % (AB[0,i],AB[0,omega.size+i]))

	file.close()
