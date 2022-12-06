#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#Importo librerias y defino variables
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#from scipy.integrate import quad

C_G=1e-18 #F
C_S=10e-18 #F
C_D=10e-18 #F
C_ES=C_S+C_D+C_G
lumo=-4.7 #eV
gamma=0.1 #eV
T=1.0 #K
Ef=-5.0 #eV


#Defino la distribución de Fermi-Dirac
k=8.617333262e-5 #eV/K
def fermi (E:float, mu:float) ->float:
    return 1/(1+np.exp((E-mu)/(k*T)))


#ALGORITMO
def algoritmo(I:float,cond:float) -> float:

    hbar=6.582119569e-16 #eV⋅s
    tau=hbar/gamma
    tau_s=tau
    tau_d=tau
    #infty=20*gamma #el límite de nuestras integrales
    q=1.602176565e-19 #C
    #Intensidad
    N0=2.0*fermi(Ef, 0)
    i,j=0,0
    alpha=1e-2
    for Vgs in np.arange(5.8,8,0.2):
        for Vds in np.arange(-0.2,0.2,0.01):
            mu_s=Ef
            mu_d=mu_s-Vds
            
            Uaux,U=0,1 #Valores para el U auxiliar y U iniciales
            eps=1e-5 #Epsilon: diferencia entre U y U auxiliar deseada
            while(abs(U-Uaux)>=eps):
                Uaux=U
                N=2.0*((tau_d*fermi(U+lumo,mu_s)) + (tau_s*fermi(U+lumo,mu_d))) / (tau_s+tau_d)
                U=-(Vgs*C_G/C_ES)-(Vds*C_D/C_ES)+(q*(N-N0)/C_ES)
                U=Uaux+alpha*(U-Uaux)
            I[i, j]=2.0*1e6*q*(fermi(U+lumo,mu_s)-fermi(U+lumo,mu_d)) / (tau_s+tau_d)
            i+=1
        i=0
        j+=1
        
    #Conductividad
    VdsCond=np.arange(-0.2,0.2,0.01) #Vds en [-0.2,0.2] V con paso 0.01 V
    VgsCond=np.arange(5.8,8,0.2) #Vgs en [5,8] V con paso 0.2 V
    for k in range(0,VgsCond.size):
        cond[:,k]=np.diff(I[:,k])/np.diff(VdsCond)



#Defino I como un array de n valores
Vds=np.arange(-0.2,0.2,0.01)#Vds en [-0.2,0.2] V con paso 0.01 V
Vgs=np.arange(5.8,8,0.2) #Vgs en [5,8] V con paso 0.2 V
I=np.zeros(shape = (np.size(Vds), np.size(Vgs)))
cond=np.zeros(shape = (np.size(Vds)-1, np.size(Vgs)))
algoritmo(I,cond)


#PLOTEO

#Intensidad
fig = plt.figure()
axInt1 = fig.gca(projection='3d')
axInt1.set_xlabel('$V_{GS} \ (V)$')
axInt1.set_ylabel('$V_{DS} \ (mV)$')
axInt1.set_zlabel('$I \ (\mu A)$')
VGS, VDS = np.meshgrid(Vgs,Vds*1e3)
axInt1.plot_surface(VGS,VDS,I,cmap=plt.cm.coolwarm,linewidth=0)
plt.savefig('Escanuela-IntSup-Ej5_1.pdf')

fig, axInt2 = plt.subplots()
axInt2.set_xlabel('$V_{GS} \ (V)$')
axInt2.set_ylabel('$V_{DS} \ (mV)$')
axInt2.set_title('$I \ (\mu A)$')
VGS, VDS = np.meshgrid(Vgs,Vds*1e3)
imInt2=axInt2.pcolormesh(VGS,VDS,I,cmap=plt.cm.coolwarm,linewidth=0,shading='auto')
fig.colorbar(imInt2, ax=axInt2)
fig.autofmt_xdate()
plt.savefig('Escanuela-IntMap-Ej5_1.pdf')


row=int((5.9-5.8)/0.01)

fig, axInt3 = plt.subplots()
axInt3.set_xlabel('$V_{DS} \ (V)$')
axInt3.set_ylabel('$I \ (\mu A)$')
axInt3.plot(Vds, I[:,row], label='$V_{GS}=5.9 \ V$')
fig.autofmt_xdate()
plt.savefig('Escanuela-IntVds-Ej5_1.pdf')

#----------------------------------

#Conductividad
Vds=np.arange(-0.2+0.01,0.2,0.01) #Para ploteo conductividad rango-0.01

fig = plt.figure()
axCond1 = fig.gca(projection='3d')
axCond1.set_xlabel('$V_{GS} \ (V)$')
axCond1.set_ylabel('$V_{DS} \ (mV)$')
axCond1.set_zlabel('$dI_{DS}/dV_{DS} \ (\mu A/V)$')
VGS, VDS = np.meshgrid(Vgs,Vds*1e3)
axCond1.plot_surface(VGS,VDS,cond,cmap=plt.cm.coolwarm,linewidth=0)
plt.savefig('Escanuela-CondSup-Ej5_1.pdf')

fig, axCond2 = plt.subplots()
axCond2.set_xlabel('$V_{GS} \ (V)$')
axCond2.set_ylabel('$V_{DS} \ (mV)$')
axCond2.set_title('$dI_{DS}/dV_{DS} \ (\mu A/V)$')
VGS, VDS = np.meshgrid(Vgs,Vds*1e3)
imCond2=axCond2.pcolormesh(VGS,VDS,cond,cmap=plt.cm.coolwarm,linewidth=0,shading='auto')
fig.colorbar(imCond2, ax=axCond2)
fig.autofmt_xdate()
plt.savefig('Escanuela-CondMap-Ej5_1.pdf')