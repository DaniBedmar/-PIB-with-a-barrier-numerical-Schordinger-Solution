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
    if (x<0.5) and (x>-0.5):
        y=Ebar
    else:
        y=0
    return y


#Discretized Schrodinger equation in n points (FROM 0 to n-1)
def Eq(n,h,x):
    F = np.zeros([n,n])
    for i in range(0,n):
        F[i,i] = -2*((h**2)*getV(x[i],Ebar) + 1)
        #print(x[i],getV(x[i]))
        if i > 1:
           F[i,i-1] = 1
           if i < n-2:
              F[i,i+1] = 1
    return F

#-------------------------
# Main program
#-------------------------
# Parameters
L = 20              #Simulation spatial lenght
xlower = -L/2.0     
xupper = L/2.0
n=4                 #numbers of levels calculated
h = 0.02            #discretization in space#
lb=0.45             #left border (will be multiplied by the interval
rb=0.55             #right border (will be multiplied by the interval
Ebar=100               #Energy of the barrier

  

#Create coordinates at which the solution will be calculated
x = np.arange(xlower,xupper,h)
#grid size (how many discrete points to use in the range [-L/2,L/2])
npoints=len(x) 

print("Using",npoints, "grid points.")

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
for k in range(0,n):
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
for v in range(0,n):
	plt.plot(x,psi[v],label=r'$\psi_v(x)$, k = ' + str(v))
	plt.title(r'$n=$'+ str(v) + r', $E_n$=' + '{:.2f}'.format(E[v]))
	plt.legend()
	plt.xlabel(r'$x$(dimensionless)')
	plt.ylabel(r'$\psi(x)$')
	plt.show()

print("Bye")
