import numpy as np
import  matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.ticker as tkr
from mpl_toolkits.mplot3d import Axes3D


"Détermination du gain expérimental"

exp_H=[12.95334765,18.9696542,34.83128596,70.98679794,95.96800753,115.6914713,123.0637935,117.8336786,115.849218,114.9872316,106.3935156,101.1043794,
95.55003001,60.48814263,40.49944059,28.97113448,21.45267398,17.29245916]
exp_f=[240,250,260,265,266,267,267.5,268,268.2,268.5,269,269.5,270,275,280,285,290,295]

m=0.0025
k=7100

H=1
khi=0.009719626
Q=51.44230769
w0=2*np.pi*267.5
d=2*np.sqrt(k*m)*khi

n=len(exp_H)

def f(x):
    return(H/np.sqrt(1+(Q*(x/w0-w0/x))**2))

F=np.linspace(245,300,50)
W=[2*3.14*f for f in F]
mod=[f(w) for w in W]

plt.plot(exp_f,exp_H,'b',label='Expérimental')
plt.plot(F,mod,'r',label='Modèle')
plt.legend(loc='best')
plt.xlabel('Fréquence en Hz')
plt.ylabel('Gain H de la fonction de transfert')
plt.show()

"Mesure de l'écart théorie/expérimental"

mod2=[f(x) for x in exp_f]
Eccart=[(1-(exp_H[i]-mod2[i])/exp_H[i])*100 for i in range(n)]


plt.plot(exp_f,Eccart)
plt.xlabel('Fréquence en Hz')
plt.ylabel('Eccart relatif en pourcentage entre les valurs théoriques et expérimentales')
plt.show()

"Optimisation en énergie"
""" 2D """

"""en faisant varier ksi et kappa"""
ksi=[0.01,0.005]
kappa=[0.1,0.3]

fig, axs=plt.subplots(2,2)

Omega_b=np.linspace(0.9,1.1,1000)
Tau_e_b=np.logspace(-2,2,20)
for i in range(len(ksi)):
    for j in range(len(kappa)):
        Omega_co=np.sqrt(1+kappa[j]**2)
        
        for T in Tau_e_b:
            y=kappa[j]**2*T*Omega_b**2/(2*abs(1-(1+2*T*ksi[i])*Omega_b**2+1j*Omega_b*(T*(1+kappa[i]**2)+2*ksi[i]-T*Omega_b**2))**2)
            axs[i,j].plot(Omega_b,y)

plt.show()

"en faisant varier kappa"
fig, axs=plt.subplots(2,1)

Omega_b=np.linspace(0.9,1.1,1000)
Tau_e_b=np.logspace(-2,2,20)
for j in range(len(kappa)):
    Omega_co=np.sqrt(1+kappa[j]**2)
    ksi = 0.0097    
    for T in Tau_e_b:
        y=kappa[j]**2*T*Omega_b**2/(2*abs(1-(1+2*T*ksi)*Omega_b**2+1j*Omega_b*(T*(1+kappa[i]**2)+2*ksi-T*Omega_b**2))**2)
        axs[j].plot(Omega_b,y)


axs[0].set_xlabel('Omega_barre')
axs[0].set_ylabel('P_barre')
axs[1].set_xlabel('Omega_barre')
axs[1].set_ylabel('P_barre')
plt.show()


""" 3D """

fig=plt.figure()
ax = fig.gca(projection='3d')

Omega_b=np.linspace(0.98,1.04,100)
Tau_e_b=np.logspace(-2,2,100)
X,Y = np.meshgrid(Omega_b,Tau_e_b)

def f(Om,T,kappa_v,ksi_v):
            return kappa_v**2*T*Om**2/(2*abs(1-(1+2*T*ksi_v)*Om**2+1j*Om*(T*(1+kappa_v**2)+2*ksi_v-T*Om**2))**2)
        

P=f(X,Y,0.285,0.0015)
surf=ax.plot_surface(X,np.log10(Y),P,
                rstride=1,
                cstride=1,
                cmap=cm.nipy_spectral_r,
                linewidth=0,
                antialiased=False)

ax.set_ylabel(r'Tau_e_barre',labelpad=7)
ax.set_xlabel(r'Omega_barre',labelpad=7)
ax.set_zlabel('P_barre',labelpad=7)

ax.xaxis.set_major_locator(tkr.AutoLocator())
ax.yaxis.set_major_locator(tkr.AutoLocator())
ax.zaxis.set_major_locator(tkr.AutoLocator())

ax.get_xaxis().get_major_formatter().set_useOffset(True)
ax.get_xaxis().get_major_formatter().set_useOffset(True)
ax.get_xaxis().get_major_formatter().set_useOffset(True)

fig.colorbar(surf, shrink=0.7, aspect=20, pad=0.12)

plt.show()

""" recherche du max sur la pulsation pour tracer puissance/R """

Omega_b=np.linspace(0.98,1.04,100)
Tau_e_b=np.logspace(-1.5,1,100)
Resistance = [x/(2.8*10**(-4)) for x in Tau_e_b]

"""valeurs théoriques"""

Puissance_th=[0.000384213,0.000422605,0.000451205,0.000475174,0.000493158,0.000508879,0.000511884,0.000535736,0.000548261,0.000557299,0.000555107,
0.000560729,0.000556118,0.000567219,0.000572083,0.000577296,0.000570538,0.000550421,0.000541323,0.000503653]
Résistance_th=[985,1195,1504,1792,2176,2675,3307,3884,4658,5585,6787,8295,9940,11970,14980,17970,21560,27180,32810,39090]

def f(Om,T,kappa_v,ksi_v):
            return kappa_v**2*T*Om**2/(2*abs(1-(1+2*T*ksi_v)*Om**2+1j*Om*(T*(1+kappa_v**2)+2*ksi_v-T*Om**2))**2)

L=[ ]
for tau in Tau_e_b:
    H = [ f(om,tau,0.285,0.00973) for om in Omega_b ]
    L.append( max (H)*4.75*10**(-5) )
    
plt.plot(Resistance,L,label='Modélisation théorique')
plt.plot(Résistance_th,Puissance_th,'-r+',label='Relevés expérimentaux')
plt.xlabel('Résistance en Ohm')
plt.ylabel('Puissance en W')
plt.legend(loc='best')
plt.grid()
plt.show()


""" recherche du max """

fig=plt.figure()
ax = fig.gca(projection='3d')

L_ksi=np.linspace(0.1,1,100)
L_kappa=np.linspace(0.005,0.5,100)

n=len(L_ksi)
p=len(L_kappa)
Pmax=np.zeros((n,p))

for ks in range(n):
    for kp in range(p):
        Pt=f(X,Y,L_kappa[kp],L_ksi[ks])
        ind = np.unravel_index(np.argmax(Pt, axis=None), Pt.shape)
        Pmax[ks,kp]=P[ind]

X1,Y1=np.meshgrid(L_ksi,L_kappa)

sf=ax.plot_surface(X1,Y1,Pmax,
                rstride=1,
                cstride=1,
                cmap=cm.nipy_spectral_r,
                linewidth=0,
                antialiased=False)

ax.set_ylabel(r'Kappa',labelpad=7)
ax.set_xlabel(r'Ksi',labelpad=7)
ax.set_zlabel('Puissance maximum',labelpad=7)

ax.xaxis.set_major_locator(tkr.AutoLocator())
ax.yaxis.set_major_locator(tkr.AutoLocator())
ax.zaxis.set_major_locator(tkr.AutoLocator())

ax.get_xaxis().get_major_formatter().set_useOffset(True)
ax.get_xaxis().get_major_formatter().set_useOffset(True)
ax.get_xaxis().get_major_formatter().set_useOffset(True)

fig.colorbar(sf, shrink=0.7, aspect=20, pad=0.12)

plt.show()
        

    
