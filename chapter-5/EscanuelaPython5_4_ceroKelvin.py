#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#Importo librerias y defino variables
import numpy as np
import matplotlib.pyplot as plt

Ec = -4.7 #eV
Ef = -5.0 #eV
T = 1 #K
m0 = 0.5*9.1e-31/1.78266192162790e-36 #eV/c^2
m0SI=0.5*9.1e-31
C_G = 0.1e-15 #F
L = 40e-9 #m
W = 3*L
C_S = 0
C_D = 0
eta0 = C_G/(C_S+C_D+C_G)
q=1.602176565e-19 #C
h=4.135667696e-15 #eV⋅s
hSI=6.62607015e-34
c=3e8
C_Q=2.0*np.pi*q*q*m0SI*W*L/(hSI*hSI)
eta = C_G/(C_S+C_D+C_G+C_Q)
print(eta)
print(C_Q)


#ALGORITMO
def algoritmo(I:float) -> float:
    cte=(4.0e6*q*W/(h*h))*np.sqrt(8.0*c*c*m0/9.0)*((eta*q)**(3.0/2.0))
    i,j=0,0
    for Vgs in np.linspace(0.3,0.5,5):
        for Vds in np.linspace(0,0.5,50):
            mu_s=Ef
            #mu_d=mu_s-Vds
            VT = (Ec-mu_s)/eta0
            if Vds<eta*(Vgs-VT):
                I[i, j]= cte*((Vgs-VT)**(3.0/2.0)-(Vgs-VT-Vds/eta)**(3.0/2.0))
            else:
                I[i, j]= cte*(Vgs-VT)**(3.0/2.0)
            #print("I[{0},{1}]={2}\n".format(i,j,I[i,j]))
            i+=1
        i=0
        #print("\n")
        j+=1



#Defino I como un array de n valores
Vds=np.linspace(0,0.5,50) 
Vgs=np.linspace(0.3,0.5,5)
I=np.zeros(shape = (np.size(Vds), np.size(Vgs)))
algoritmo(I)


#PLOTEO

#Intensidad
#Ploteo
fig, axInt = plt.subplots()

Vds=np.linspace(0,0.5,50)#El de la intensidad
axInt.set_xlabel('$V_{DS} (V)$')
axInt.set_ylabel('$I_{DS} (\mu A)$')
for i in range(0,5):
    axInt.plot(Vds,I[:,i])

plt.savefig('Escanuela-T0K-Ej5_3.pdf')