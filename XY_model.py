## IMPORT LIBRARY
import numpy as np
import scipy as sp
import matplotlib as mpl
import matplotlib.pyplot as plt
import pylab
import random
import matplotlib.image as mpimg
from matplotlib.legend_handler import HandlerLine2D
import math
from scipy.optimize import curve_fit
from numpy import pi





## applying Metropolis algorithm
# input: T/temperature
#        S/spins configuration(in 1d list)
#        H/ecternal field.default value=0



def sweep(T,S,H=None):
    N = int(np.size(S))
    L = int(np.sqrt(N))
    H = 0. if H==None else H
    delta = 0.5 #maximum cnange of theta
    nbr = {i : ((i // L) * L + (i + 1) % L, (i + L) % N,
                (i // L) * L + (i - 1) % L, (i - L) % N) \
                                        for i in list(range(N))}
    beta = 1.0 / T
    for step in list(range(N)):#one sweep in defined as N attempts of flip
        k = np.random.randint(0, N - 1)#randomly choose a spin
        for m in list(range(3)):
            energy_i = -sum(np.cos(S[k]-S[n]) for n in nbr[k]) 
            dtheta = delta * (np.random.random()*2-1)
            S_temp = S[k] + dtheta
            energy_f = -sum(np.cos(S_temp-S[n]) for n in nbr[k]) 
            delta_E = energy_f - energy_i
            if np.random.uniform(0.0, 1.0) < math.exp(-beta * delta_E):
                S[k] += dtheta
    return S

## calculate the energy through the configuration  
#  input: S/spin configuration in list
#         H/external field, defult 0
def energy(S,H=None):
    N=int(np.size(S))
    L = int(np.sqrt(N))
    E=0
    H = 0. if H==None else H
    nbr = {i : ((i // L) * L + (i + 1) % L, (i + L) % N,
                (i // L) * L + (i - 1) % L, (i - L) % N) \
                                        for i in list(range(N))}
    for kk in list(range(0,N)): #calculate energy per spin
        E += -sum(np.cos(S[kk]-S[n]) for n in nbr[kk])/N/2#nearst neighbor of kth spin
    return E

## thermodynamic quantities of 2D ising using Metropolis algorithm
# input: T/temperature
#        S/spin configuration
#        nsweeps/number of sweps
#        pr/just ingore it
def thermo(T,S,nsweeps,H=None,pr=None):
    H = 0. if H==None else H
    #size of system
    N = int(np.size(S)) 
    L = int(np.sqrt(N))
    beta=1.0/T
    list_M=[]
    list_E=[]
    list_Cv=[]
    list_chi=[]
    for k in list(range(nsweeps)):
        S = sweep(T,S,H)     
        #list_M.append(np.abs(np.sum(S)/N))
        list_E.append(energy(S))
    Energy=np.sum(list_E)/nsweeps
    #magnetization=np.sum(list_M)/nsweeps
    energy2=np.sum(np.power(list_E,2))/nsweeps
    #magetization2=np.sum(np.power(list_M,2))/nsweeps
    #chi=(magetization2-magnetization**2)*beta
    Cv=(energy2-Energy**2)*beta**2
    #return S,Energy,magnetization,chi,Cv
    return S,Energy,Cv

## To see thermoquantities evolve as we cooling the systems down
# input: L/size of grid
#        sample/'log' or 'lin',mean linear sampled T or log sampled( centered at critical point)
def cooling(L,T_init=2.5,T_final=0.01):
    N=L**2
    nsweeps=400
    # initialize spins. Orientations are taken from 0 - 2pi randomly.
    #initialize spin configuration  
    S = np.random.random((1,N))*2*pi
    S=S[0][:]
    list_T=list(np.linspace(T_init,T_final,50))
    list_E=[]
    list_M=[]
    list_Cv=[]
    list_chi=[]
    for t in list_T:
        #nsweep = 300 if t > 0.8 else 
        S,energy_t,Cv_t=thermo(t,S,nsweeps)
        #S,energy_t,magnetization_t,chi_t,Cv_t=thermo(t,S,nsweeps)
        list_E.append(energy_t)
        #list_M.append(magnetization_t)
        list_Cv.append(Cv_t)
        #list_chi.append(chi_t)
    #plt.plot(list_T,list_E)
    '''
    print(sample)
    plt.plot(list_T,list_chi,'.')
    plt.ylabel(r'$\chi$')
    plt.xlabel('T')
    plt.show()  
    plt.plot(list_T,list_M,'.')
    plt.ylabel(r'$\langle |m|\rangle$')
    plt.xlabel('T')
    plt.show()
    return list_E,list_M,list_Cv,list_chi,list_T
    '''
    plt.plot(list_T,list_Cv,'.')
    plt.ylabel(r'$C_v$')
    plt.xlabel('T')
    plt.show()
    plt.plot(list_T,list_E,'.')
    plt.ylabel(r'$\langle E \rangle$')
    plt.xlabel('T')
    plt.show()
    return list_E,list_Cv,list_T


## convert configuration inz list to matrix form
def list2matrix(S):
    N=int(np.size(S))
    L = int(np.sqrt(N))
    S=np.reshape(S,(L,L))
    return S

## visulize a configurtion
#  input：S/ spin configuration in list form
def visulize(S,T,i=None):
    N=int(np.size(S))
    L = int(np.sqrt(N))
    S = list2matrix(S)
    X, Y = np.meshgrid(np.arange(0,L ),np.arange(0, L))
    U = np.cos(S)
    V = np.sin(S)
    plt.figure(1,figsize=(8,8), dpi=100)
    #plt.title('Arrows scale with plot width, not view')
    Q = plt.quiver(X, Y, U, V, units='width')
    qk = plt.quiverkey(Q, 0.1, 0.1, 1, r'$2 \frac{m}{s}$', labelpos='E',
                   coordinates='figure')
    plt.title('T=%.2f'%T+',Size='+str(L)+',Sweeps=%i'%i)
    #if i != None:
    #    plt.savefig('C:/Users/liang/Desktop/Computating_Project/XY MODEL/20_0.5/'+str(i)+'.png')
    plt.show()

def visulizeXY(T=0.5,L=10):
    N = L**2
    S = np.random.random((1,N))*2*pi#initialize spin configuration
    S=S[0][:]
    nsweeps=300
    magnetization=[]
    stime=list(range(nsweeps))
    for k in stime:
        a=visulize(S,T,k)
        S=sweep(T,S,H=None)
        #plt.savefig(a,'Size='+str(L)+','+str(k)+'.png')
        
def visulizeXY_equ(T=0.5,L=5):
    N = L**2
    S = np.random.random((1,N))*2*pi#initialize spin configuration
    S=S[0][:]
    nsweeps=1000
    magnetization=[]
    stime=list(range(nsweeps))
    for k in stime:
        S=sweep(T,S,H=None)
    visulize(S,T)
    
    
if __name__ == '__main__':
    visulizeXY(T=0.5,L=5)