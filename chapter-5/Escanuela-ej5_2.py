#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#Importo librerias y defino variables
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

Ec1 = -4.7 #eV
Ec2 = -4.6 #eV
Ef = -5.0 #eV
T = 50 #K
m0 = 9.1e-31/1.78266192162790e-36 #eV/c^2
C_G = 50e-18 #F
L = (100e-9)*(1e15) #fm
q2_entre_C_ES = 1 #eV



#Defino la distribución de Fermi-Dirac
def fermi (E:float, mu:float) ->float:
    k=8.617333262e-5 #eV/K
    return 1/(1+np.exp((E-mu)/(k*T)))


#ALGORITMO
def algoritmo(I:float) -> float:
    h=4.135667696e-15 #eV⋅s
    q=1.602176565e-19 #C
    hc = 1239.841984e6 #eV*fm

    C_ES=q/q2_entre_C_ES #F

    N01 = quad(lambda E: np.sqrt((2.0*m0)/abs(E-Ec1)) * fermi(E, 0.0) * np.heaviside(E-Ec1,0), -np.inf, np.inf)[0]
    print("\n\n\nN01={0}\n\n\n".format(N01))
    N02 = quad(lambda E: np.sqrt((2.0*m0)/abs(E-Ec2)) * fermi(E, 0.0) * np.heaviside(E-Ec2,0), -np.inf, np.inf)[0]
    print("\n\n\nN02={0}\n\n\n".format(N02))
    N0 = (L/hc)*(N01+N02)
    print("\n\n\nN0={0}\n\n\n".format(N0))

    i,j=0,0
    alpha=1e-2
    for Vgs in np.linspace(0.3,0.5,5):
        for Vds in np.linspace(0,0.5,50):
            mu_s=0
            mu_d=-Vds
            
            Uaux,U=0,1 #Valores para el U auxiliar y U iniciales
            eps=1e-5 #Epsilon: diferencia entre U y U auxiliar deseada
            while(abs(U-Uaux)>=eps):
                Uaux=U
                N1 = quad(lambda E:  np.sqrt((2.0*m0)/abs(E-Ec1-U+Ef))* np.heaviside(E-Ec1-U+Ef,0) * (fermi(E,mu_s)+fermi(E,mu_d)), -np.inf, np.inf)[0]
                #print("N1={0}\n".format(N1))
                N2 = quad(lambda E: np.sqrt((2.0*m0)/abs(E-Ec2-U+Ef)) * np.heaviside(E-Ec2-U+Ef,0) * (fermi(E,mu_s)+fermi(E,mu_d)), -np.inf, np.inf)[0]
                #print("N2={0}\n".format(N2))
                N = (L/hc) * (N1 + N2)
                #print("N={0}\n".format(N))
                #ppppp=input()
                #--------------
                #AUTOCONSISTENTE:
                #U = -Vgs*C_G/C_ES+q2_entre_C_ES*(N-N0) 

                #NO AUTOCONSISTENTE:
                U=0
                #--------------
                U = Uaux+alpha*(U-Uaux)
            #print("U=",U)
            I1 = quad(lambda E: np.heaviside(E-Ec1-U+Ef,0) * (fermi(E,mu_s)-fermi(E,mu_d)) , -np.inf, np.inf)[0]
            #print("I1={0}".format(I1))
            I2 = quad(lambda E: np.heaviside(E-Ec2-U+Ef,0) * (fermi(E,mu_s)-fermi(E,mu_d)) , -np.inf, np.inf)[0]
            #print("I2={0}".format(I2))
            I[i, j]= (2.0e6*q/h) * (I1 + I2)
            print("I[{0},{1}]={2}\n".format(i,j,I[i,j]))
            i+=1
        i=0
        print("\n")
        j+=1



#Defino I como un array de n valores
Vds=np.linspace(0,0.5,50) 
Vgs=np.linspace(0.3,0.5,5)
I=np.zeros(shape = (np.size(Vds), np.size(Vgs)))
algoritmo(I)


#PLOTEO

#Intensidad
#Ploteo
fig, axInt = plt.subplots() #Una al lado de la otra

Vds=np.linspace(0,0.5,50)#El de la intensidad
axInt.set_xlabel('$V_{DS} (V)$')
axInt.set_ylabel('$I_{DS} (\mu A)$')
axInt.plot(Vds,I[:,4])

plt.savefig('Escanuela-FIGURA-Ej5_2.pdf')