#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#Importo librerias y defino variables
import numpy as np
import matplotlib.pyplot as plt

eta=0.5
Ef=-5 #eV
homo=-5.5 #eV
gamma_s=0.1 #eV
gamma_d=0.1 #eV
T=298 #K

#Convierto de gamma a tau
hbar=6.582119569e-16 #eV⋅s
tau_s=hbar/gamma_s
tau_d=hbar/gamma_d

#Defino la distribución de Fermi-Dirac
k=8.617333262e-5 #eV/K
def fermi (E:float, mu:float, T:float) ->float:
    return 1/(1+np.exp((E-mu)/(k*T)))

#Arrays de valores
n=100
rangoV=2
Vds=np.linspace(-rangoV,rangoV,n)
mu_s=eta*Vds
mu_d=mu_s-Vds

#ALGORITMO
q=1.602176565e-19
#Si C_ES -> oo
U=0
I1=2*q*(fermi(U+homo-Ef,mu_s,T)-fermi(U+homo-Ef,mu_d,T))/(tau_s+tau_d)
conduct1=np.diff(I1)/np.diff(Vds)
#Si q^2/C_ES = 1 eV
q2_entre_C_ES=1
N0=2*fermi(Ef,0,T) #=2

I2=np.linspace(-rangoV,rangoV,n) #Defino I2 como un array de n valores
contador=0
alpha=1e-2
for Vds in np.arange(-rangoV,rangoV,2*rangoV/n):
    mu_s=eta*Vds
    mu_d=mu_s-Vds

    Uaux,U=0,1 #Valores para el U auxiliar y U iniciales
    eps=1e-5 #Epsilon: diferencia entre U y U auxiliar deseada
    while(abs(U-Uaux)>=eps):
        Uaux=U
        N=2*((tau_d*fermi(U+homo-Ef,mu_s,T))+(tau_s*fermi(U+homo-Ef,mu_d,T)))/(tau_s+tau_d)
        U=q2_entre_C_ES*(N-N0)
        U=Uaux+alpha*(U-Uaux)
    I2[contador]=2*q*(fermi(U+homo-Ef,mu_s,T)-fermi(U+homo-Ef,mu_d,T))/(tau_s+tau_d)
    contador+=1
    
Vds=np.linspace(-rangoV,rangoV,n) #Reconvierto tipo de Vds a linspace
conduct2=np.diff(I2)/np.diff(Vds)
#Ploteo
fig, (axInt, axCond) = plt.subplots(1,2) #Una al lado de la otra

Vds=np.linspace(-rangoV,rangoV,n) #El de la intensidad
axInt.set_xlabel('Tensión (V)')
axInt.set_ylabel('Corriente ($\mu$ A)')
axInt.plot(Vds,I1*1e6, label='$q^2/C_{ES}=0\, eV$')
axInt.plot(Vds,I2*1e6, label='$q^2/C_{ES}=1\, eV$')
axInt.legend(loc='upper left')

Vds=np.linspace(-rangoV,rangoV,n-1) #El de la conductancia
axCond.set_ylabel('dI/dV ($\mu A / V)$')
axCond.set_xlabel('Tensión (V)')
axCond.plot(Vds,conduct1*1e6)
axCond.plot(Vds,conduct2*1e6)

plt.tight_layout()
plt.savefig('NanoEj3_9.pdf')
