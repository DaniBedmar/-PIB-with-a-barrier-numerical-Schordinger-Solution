#-------------------------------------------------------------------------
#funcion para plotear el potencial de la barrera y los niveles energeticos
#-------------------------------------------------------------------------
#
#Importamos librerias
#
import numpy as np
import matplotlib.pyplot as plt

#Funcion para que plotee lo q queremos
def plotEnergy (L,Ebar,n,lb,rb,h):
    x=np.arange(0,L,h)
    Elvls=[]
    Barrier=[]
    for i in range(n):
        Elvls.append(((i+1)**2)*(np.pi**2)/2)
    for i in range(len(x)):
        if (x[i]<rb*L) and (x[i]>lb*L):
            Barrier.append(Ebar)
        else:
            Barrier.append(0)
    plt.plot(x,Barrier)
    for i in range(n):
        plt.axhline(y=((i+1)**2)*(np.pi**2)/2, color='orange')
    plt.show()
    return Barrier

    #print(x)
        

#print(Barrier)
L=20
h=0.001
plotEnergy(L,100,4,0.45,0.55,h)
#x=np.arange(0,L,h)
#plt.plot(x,Barrier)
#plt.show
