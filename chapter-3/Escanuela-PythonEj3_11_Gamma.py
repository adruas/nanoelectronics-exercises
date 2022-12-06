#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#Importo librerias y defino variables
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

eta=0.5
homo=-5.5 #eV
T=298 #K
Ef=-5 #eV
n=100 #Numero de puntos
rangoV=4  
Vds=np.linspace(-rangoV,rangoV,n)


#Defino la distribución de Fermi-Dirac
k=8.617333262e-5 #eV/K
def fermi (E:float, mu:float, T:float) ->float:
    return 1/(1+np.exp((E-mu)/(k*T)))

#Defino la función densidad de estados (Lorentziana)
def DoS(E:float, gamma:float, homo:float, Ef:float) -> float:
    epsilon=homo-Ef
    v=E-epsilon
    u=gamma/2
    return gamma/(np.pi*((v*v)+(u*u)))

#ALGORITMO
def algoritmo(n:int,rangoV:int,q2_entre_C_ES:float,I1:float, I2:float ,Ef:float)->float:

    #Lorentziana GammaS=GammaD
    #-----------
    gamma_s1=0.1 #eV
    gamma_d1=0.1 #eV
    gamma1=gamma_d1+gamma_s1
    hbar=6.582119569e-16 #eV⋅s
    tau_s1=hbar/gamma_s1
    tau_d1=hbar/gamma_d1
    infty1=20*gamma1 #el límite de nuestras integrales
    q=1.602176565e-19
    N01=quad(lambda E: fermi(E, 0, T) * DoS(E, gamma1, homo, Ef), -infty1, infty1)[0]
    i=0
    alpha=1e-2
    for Vds in np.arange(-rangoV,rangoV,2*rangoV/n):
        mu_s=eta*Vds
        mu_d=mu_s-Vds
        
        Uaux,U=0,1 #Valores para el U auxiliar y U iniciales
        eps=1e-5 #Epsilon: diferencia entre U y U auxiliar deseada
        while(abs(U-Uaux)>=eps):
            Uaux=U
            N1=quad(lambda E: DoS(E-U, gamma1, homo, Ef) * ( (tau_d1*fermi(E,mu_s,T)) 
                + (tau_s1*fermi(E,mu_d,T)) ) / (tau_s1+tau_d1), -infty1, infty1)[0]
            U=q2_entre_C_ES*(N1-N01)
            U=Uaux+alpha*(U-Uaux)
        I1[i] =1e6*q*quad(lambda E: DoS(E-U, gamma1, homo, Ef) * ( fermi(E,mu_s,T)
             - fermi(E,mu_d,T) ) / (tau_s1+tau_d1), -infty1, infty1)[0]
        i+=1

    #Lorentziana GammaS=2GammaD
    #-----
    gamma_s2=0.2 #eV
    gamma_d2=0.1 #eV
    gamma2=gamma_d2+gamma_s2
    tau_s2=hbar/gamma_s2
    tau_d2=hbar/gamma_d2
    infty2=20*gamma2 #el límite de nuestras integrales
    N02=quad(lambda E: fermi(E, 0, T) * DoS(E, gamma2, homo, Ef), -infty2, infty2)[0]
    i=0
    alpha=1e-2
    for Vds in np.arange(-rangoV,rangoV,2*rangoV/n):
        mu_s=eta*Vds
        mu_d=mu_s-Vds
        
        Uaux,U=0,1 #Valores para el U auxiliar y U iniciales
        eps=1e-5 #Epsilon: diferencia entre U y U auxiliar deseada
        while(abs(U-Uaux)>=eps):
            Uaux=U
            N2=quad(lambda E: DoS(E-U, gamma2, homo, Ef) * ( (tau_d2*fermi(E,mu_s,T)) 
                + (tau_s2*fermi(E,mu_d,T)) ) / (tau_s2+tau_d2), -infty2, infty2)[0]
            U=q2_entre_C_ES*(N1-N01)
            U=Uaux+alpha*(U-Uaux)
        I2[i] =1e6*q*quad(lambda E: DoS(E-U, gamma2, homo, Ef) * ( fermi(E,mu_s,T)
             - fermi(E,mu_d,T) ) / (tau_s2+tau_d2), -infty2, infty2)[0]
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
axInt.plot(Vds,I1, label='$\Gamma_S=\Gamma_D$')
axInt.plot(Vds,I2, label='$\Gamma_S=2\Gamma_D$')
axInt.legend(loc='upper left')

Vds=np.linspace(-rangoV,rangoV,n-1) #El de la conductancia
axCond.set_ylabel('dI/dV ($\mu A / V)$')
axCond.set_xlabel('Tensión (V)')
axCond.plot(Vds,conduct1)
axCond.plot(Vds,conduct2)

plt.tight_layout()
plt.savefig('Escanuela-FIGURA-Ej3_11DifGamma.pdf')