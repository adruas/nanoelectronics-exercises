#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#Importo librerias y defino variables
import numpy as np
import matplotlib.pyplot as plt

eta=0.5
homo=-5.5 #eV
lumo=-1.5 #eV
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

#ALGORITMO
q=1.602176565e-19
#Defino una función para el algoritmo
def algoritmo(n:int,rangoV:int,q2_entre_C_ES:float,I1:float,Ef:float)->float:
    N0=2*fermi(Ef,0,T) #=2  
    contador=0
    alpha=1e-2
    for Vds in np.arange(-rangoV,rangoV,2*rangoV/n):
        mu_s=eta*Vds
        mu_d=mu_s-Vds
        
        Uaux,U=0,1 #Valores para el U auxiliar y U iniciales
        eps=1e-5 #Epsilon: diferencia entre U y U auxiliar deseada
        while(abs(U-Uaux)>=eps):
            Uaux=U
            N=2*((tau_d*fermi(U+homo-Ef,mu_s,T))+(tau_s*fermi(U+homo-Ef,mu_d,T))
                 +(tau_d*fermi(U+lumo-Ef,mu_s,T))+(tau_s*fermi(U+lumo-Ef,mu_d,T)))/(tau_s+tau_d)
            U=q2_entre_C_ES*(N-N0)
            U=Uaux+alpha*(U-Uaux)
        I1[contador]=2*q*((fermi(U+homo-Ef,mu_s,T)-fermi(U+homo-Ef,mu_d,T))
            +(fermi(U+lumo-Ef,mu_s,T)-fermi(U+lumo-Ef,mu_d,T)))/(tau_s+tau_d)
        contador+=1


n=100 #Numero de puntos
rangoV=6  
Vds=np.linspace(-rangoV,rangoV,n)

#Si q^2/C_ES = 1 eV y E_f=-2.5 eV
q2_entre_C_ES=1
Ef=-2.5 #eV
I1=np.linspace(-rangoV,rangoV,n) #Defino I como un array de n valores
algoritmo(n,rangoV,q2_entre_C_ES,I1,Ef)
conduct1=np.diff(I1)/np.diff(Vds)

#Si q^2/C_ES = 1 eV, E_F=-3.5
q2_entre_C_ES=1
Ef=-3.5 #eV
I2=np.linspace(-rangoV,rangoV,n) #Defino I como un array de n valores
algoritmo(n,rangoV,q2_entre_C_ES,I2,Ef)
conduct2=np.diff(I2)/np.diff(Vds)
#Si q^2/C_ES = 1 eV, E_F=-5
q2_entre_C_ES=1
Ef=-5 #eV
I3=np.linspace(-rangoV,rangoV,n) #Defino I como un array de n valores
algoritmo(n,rangoV,q2_entre_C_ES,I3,Ef)
conduct3=np.diff(I3)/np.diff(Vds)

#Ploteo
fig, (axInt, axCond) = plt.subplots(1,2) #Una al lado de la otra

Vds=np.linspace(-rangoV,rangoV,n)#El de la intensidad
axInt.set_xlabel('Tensión (V)')
axInt.set_ylabel('Corriente ($\mu$ A)')
axInt.plot(Vds,I1*1e6, label='$E_F=-2.5\, eV$')
axInt.plot(Vds,I2*1e6, label='$E_F=-3.5\, eV$')
axInt.plot(Vds,I3*1e6, label='$E_F=-5\, eV$')
axInt.legend(loc='upper left')

Vds=np.linspace(-rangoV,rangoV,n-1) #El de la conductancia
axCond.set_ylabel('dI/dV ($\mu A / V)$')
axCond.set_xlabel('Tensión (V)')
axCond.plot(Vds,conduct1*1e6)
axCond.plot(Vds,conduct2*1e6)
axCond.plot(Vds,conduct3*1e6)

plt.tight_layout()
plt.savefig('NanoEj3_10.pdf')
