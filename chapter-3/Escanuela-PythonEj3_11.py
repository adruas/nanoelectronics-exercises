#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#Importo librerias y defino variables
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

eta=0.5
homo=-5.5 #eV
gamma_s=0.1 #eV
gamma_d=0.1 #eV
T=298 #K
Ef=-5 #eV
n=100 #Numero de puntos
rangoV=4  
Vds=np.linspace(-rangoV,rangoV,n)

#Convierto de gamma a tau
hbar=6.582119569e-16 #eV⋅s
tau_s=hbar/gamma_s
tau_d=hbar/gamma_d

#Defino la distribución de Fermi-Dirac
k=8.617333262e-5 #eV/K
def fermi (E:float, mu:float, T:float) ->float:
    return 1/(1+np.exp((E-mu)/(k*T)))

#Defino la función densidad de estados (Lorentziana)
gamma=gamma_d+gamma_s
def DoS(E:float, gamma:float, homo:float, Ef:float) -> float:
    epsilon=homo-Ef
    v=E-epsilon
    u=gamma/2
    return gamma/(np.pi*((v*v)+(u*u)))

#ALGORITMO
def algoritmo(n:int,rangoV:int,q2_entre_C_ES:float,I1:float, I2:float ,Ef:float)->float:

    #Lorentziana
    #-----------
    infty=20*gamma #el límite de nuestras integrales
    q=1.602176565e-19
    N01=quad(lambda E: fermi(E, 0, T) * DoS(E, gamma, homo, Ef), -infty, infty)[0]
    i=0
    alpha=1e-2
    for Vds in np.arange(-rangoV,rangoV,2*rangoV/n):
        mu_s=eta*Vds
        mu_d=mu_s-Vds
        
        Uaux,U=0,1 #Valores para el U auxiliar y U iniciales
        eps=1e-5 #Epsilon: diferencia entre U y U auxiliar deseada
        while(abs(U-Uaux)>=eps):
            Uaux=U
            N1=quad(lambda E: DoS(E-U, gamma, homo, Ef) * ( (tau_d*fermi(E,mu_s,T)) 
                + (tau_s*fermi(E,mu_d,T)) ) / (tau_s+tau_d), -infty, infty)[0]
            U=q2_entre_C_ES*(N1-N01)
            U=Uaux+alpha*(U-Uaux)
        I1[i] =1e6*q*quad(lambda E: DoS(E-U, gamma, homo, Ef) * ( fermi(E,mu_s,T)
             - fermi(E,mu_d,T) ) / (tau_s+tau_d), -infty, infty)[0]
        i+=1

    #Dirac
    #-----
    N02=2*fermi(Ef,0,T)
    i=0
    for Vds in np.arange(-rangoV,rangoV,2*rangoV/n):
        mu_s=eta*Vds
        mu_d=mu_s-Vds

        Uaux,U=0,1 #Valores para el U auxiliar y U iniciales
        eps=1e-5 #Epsilon: diferencia entre U y U auxiliar deseada
        while(abs(U-Uaux)>=eps):
            Uaux=U
            N2=2*( (tau_d*fermi(U+homo-Ef,mu_s,T)) + (tau_s*fermi(U+homo-Ef,mu_d,T)) ) / (tau_s+tau_d)
            U=q2_entre_C_ES*(N2-N02)
            U=Uaux+alpha*(U-Uaux)
        I2[i]=1e6*2*q*( fermi(U+homo-Ef,mu_s,T) - fermi(U+homo-Ef,mu_d,T) ) / (tau_s+tau_d)
        i+=1



#Si q^2/C_ES = 1 eV
q2_entre_C_ES=1
#Defino I como un array de n valores
I1 = np.linspace(-rangoV,rangoV,n)
I2 = np.linspace(-rangoV,rangoV,n)
algoritmo(n,rangoV,q2_entre_C_ES,I1,I2,Ef)
conduct1=np.diff(I1)/np.diff(Vds)
conduct2=np.diff(I2)/np.diff(Vds)



#Ploteo
fig, (axInt, axCond) = plt.subplots(1,2) #Una al lado de la otra

Vds=np.linspace(-rangoV,rangoV,n)#El de la intensidad
axInt.set_xlabel('Tensión (V)')
axInt.set_ylabel('Corriente ($\mu$ A)')
axInt.plot(Vds,I1, label='Lorentziana')
axInt.plot(Vds,I2, label='Dirac')
axInt.legend(loc='upper left')

Vds=np.linspace(-rangoV,rangoV,n-1) #El de la conductancia
axCond.set_ylabel('dI/dV ($\mu A / V)$')
axCond.set_xlabel('Tensión (V)')
axCond.plot(Vds,conduct1)
axCond.plot(Vds,conduct2)

plt.tight_layout()
plt.savefig('Escanuela-FIGURA-Ej3_11.pdf')