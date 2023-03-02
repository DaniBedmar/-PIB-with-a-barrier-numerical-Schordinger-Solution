# ----------------------------
# Quantum Harmonic Oscillator
# ----------------------------
# Finite differences method as developed by Truhlar JCP 10 (1972) 123-132
#
# code by Jordi Faraudo
#  
#
import numpy as np
import matplotlib.pyplot as plt

#Potential as a function of position
def getV(x,Ebar):
    if (x==xlower) or (x==xupper):
        potvalue = 10000000000
    if (x>-0.05) and (x<0.05):
        potvalue = Ebar
    else:
        potvalue = 0
    return potvalue

#Discretized Schrodinger equation in n points (FROM 0 to n-1)
def Eq(n,h,x):
    F = np.zeros([n,n])
    for i in range(0,n):
        F[i,i] = -2*((h**2)*getV(x[i],Ebar) + 1)
        if i > 1:
           F[i,i-1] = 1
           if i < n-2:
              F[i,i+1] = 1
    return F

#-------------------------
# Main program
#-------------------------
#PARAMETERS DEFINITION
# Interval for calculating the wave function [-L/2,L/2]
L = 1
xlower = -L/2.0
xupper = L/2.0

#Energy barruer
Ebar=100

#Discretization options
h = 0.002  #discretization in space

#Create coordinates at which the solution will be calculated
x = np.arange(xlower,xupper,h)
#grid size (how many discrete points to use in the range [-L/2,L/2])
npoints=len(x) 

print("Using",npoints, "grid points.")
#Barrier=[]
#for i in range(npoints):
    #if (x>-0.05) and (x<0.05):
     #   Barrier.append(100)
    #else:
     #   Barrier.append(0)

        
#Calculation of discrete form of Schrodinger Equation
print("Calculating matrix...")
F=Eq(npoints,h,x)

#diagonalize the matrix F
print("Diagonalizing...")
eigenValues, eigenVectors = np.linalg.eig(F)

#Order results by eigenvalue
# w ordered eigenvalues and vs ordered eigenvectors
idx = eigenValues.argsort()[::-1]   
w = eigenValues[idx]
vs = eigenVectors[:,idx]

#Energy Level
E = - w/(2.0*h**2)

#Print Energy Values
print("RESULTS:")
for k in range(0,5):
	print("State ",k," Energy = %.2f" %E[k])

#Init Wavefunction (empty list with npoints elements)
psi = [None]*npoints

#Calculation of normalised Wave Functions
for k in range(0,len(w)):
	psi[k] = vs[:,k]
	integral = h*np.dot(psi[k],psi[k])
	psi[k] = psi[k]/integral**0.5

#Plot Wave functions
print("Plotting")

#v = int(input("\n Quantum Number (enter 0 for ground state):\n>"))
for v in range(0,5):
	plt.plot(x,psi[v],label=r'$\psi_v(x)$, k = ' + str(v))
	plt.axvline(x = 0.05, color = 'orange', label = 'axvline - full height')
	plt.axvline(x = -0.05, color = 'orange', label = 'axvline - full height')
	plt.title(r'$n=$'+ str(v) + r', $E_n$=' + '{:.2f}'.format(E[v]))
	#plt.axhline(x=-0.05, color = 'orange')
	#plt.axhline(x=0.05, color = 'orange')
	plt.legend()
	plt.xlabel(r'$x$(dimensionless)')
	plt.ylabel(r'$\psi(x)$')
	plt.show()

print("Bye")
