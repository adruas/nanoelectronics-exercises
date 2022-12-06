#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#Importo librerias y defino variables
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

Ec = -4.7 #eV
Ef = -5.0 #eV
T = 298 #K
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


#Defino la distribución de Fermi-Dirac
def fermi (E:float, mu:float) ->float:
    k=8.617333262e-5 #eV/K
    return 1/(1+np.exp((E-mu)/(k*T)))


#ALGORITMO
def algoritmo(I:float) -> float:
    cte=4.0e6*q*W/(h*h)
    i,j=0,0
    for Vgs in np.linspace(0.3,0.5,5):
        for Vds in np.linspace(0,0.5,50):
            mu_s=Ef
            mu_d=mu_s-Vds
            VT = (Ec-mu_s)/eta0
            I[i, j] = cte*quad(lambda E:np.sqrt(abs(2.0*m0*c*c*(E+(eta*(Vgs-VT)))))*(fermi(E,mu_s)-fermi(E,mu_d)), mu_s-eta*(Vgs-VT), np.inf)[0]
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
for i in range(0,4):
    axInt.plot(Vds,I[:,i]) #NO LOG
plt.savefig('Escanuela-T298K-Ej5_3.pdf')

fig, axlogInt = plt.subplots()
Vds=np.linspace(0,0.1,10)
Ilog=np.zeros(shape = (np.size(Vds), np.size(Vgs)))
axlogInt.set_yscale('log')
for i in range(0,9): #9 primeras Vds
    for j in range(0,4):
        Ilog[i,j]=np.log(I[i,j]) #primeros valores de intensidad
axlogInt.set_xlabel('$V_{DS} (V)$')
axlogInt.set_ylabel('$\ln (I_{DS}) (\mu A)$')
for k in range(0,4):
    axlogInt.plot(Vds,Ilog[:,k]) #log
plt.tight_layout()
plt.savefig('Escanuela-logT298K-Ej5_3.pdf')

