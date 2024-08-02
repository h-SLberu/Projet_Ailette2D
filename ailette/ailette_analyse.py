# Importation des modules
import numpy as np
import matplotlib.pyplot as plt
import pytest
import ailette_corr
try:
    from ailette_fct import *
except:
    pass

#------------------------------------------------------------------------------
# Code principal pour l'analyse des résultats
# Il faudra faire appel aux fonctions programmées dans ailette_fct.py afin de
# modéliser les pertes de chaleurs selon différentes combinaisons géométriques.
# Ensuite, il faudra identifier le meilleur rendement.
#------------------------------------------------------------------------------

# Assignation des paramètres
# ATTENTION! Ne pas changer le nom des attributs
class parametres():
    L = 0.1       # [m] Longueur
    D = 0.0025       # [m] Diamètre
    k = 45       # [W/m*K] Conductivité thermique
    T_a = 25+273     # [K] Température de l'air ambiant
    T_w = 125+273     # [K] Température de la paroi
    h = 150       # [W/m^2*K] Coefficient de convection
    N = 1000      # [-] Nombre de points en z

prm = parametres()

# Appel des fonctions pour le calcul du profil de température
Tz, z = mdf(prm)

#analytique
L = prm.L
D = prm.D
k = prm.k
T_a = prm.T_a
T_w = prm.T_w
h = prm.h
N = prm.N

m = np.sqrt( (4*h) / (k*D) )
T = (np.cosh(m * (L-z))*(T_w - T_a)) / np.cosh(m * L) + T_a

# Graphique
plt.figure()
plt.plot(z, Tz, 'r', label='solution obtenue')
plt.plot(z, T, 'k--', label='solution analytique')

plt.title('Profil de température ')
plt.xlabel('Disrance z (m)')
plt.ylabel('Température (K)')
plt.legend()
plt.savefig('ailette.png',dpi=300)


#plt.show()

# Calcul de la dissipation pour chaque géométrie

#dissipation = inte(Tz, z, prm)

eps = 1
prm.N = 15

while True:
    
    p1, p2 = mdf(prm)
    p2 = np.array(p2)
    
    sln = (prm.T_w - prm.T_a)*np.cosh(m*(prm.L-p2))/np.cosh(m*prm.L) + prm.T_a
    
    q = inte(p1, p2, prm)
    qref = inte(sln, p2, prm)
    eps = np.abs((q - qref)/qref)
    
    if eps < 1e-2:
        break
    else:
        prm.N += 1
print(prm.N)

npts = 100
Dim_min_found = False
colors = ['-r', '-b', '-g', '-m', '-c', '-y']
aL = np.array([0.005, 0.0075, 0.01, 0.0125, 0.015])
aD = np.linspace(0.001, 0.02, npts)
aQ_e = np.empty(npts)
q_lim = np.linspace(10, 10, npts)

for L in aL:
    aq = np.array([])
    
    for D in aD:
        prm.L = L
        prm.D = D

        p1, p2 = mdf(prm)
        q = inte(p1, p2, prm)
        
        # Che
        if q >= 10 and Dim_min_found is False:
            
            D_min = D
            L_min = L
            Dim_min_found = True
            
        aq = np.append(aq, q)
        
    aQ_e = np.vstack((aQ_e, aq))
    
    # Enlever 1a ligne crée par np.empty()
    aQ = np.delete(aQ_e, 0, axis=0)

# Graphique
plt.figure()
for i in range(5):
    # Graphique
    plt.plot(aD, aQ[i], colors[i])

plt.plot(aD, q_lim, colors[5])
plt.title('Chaleur dissipée en fonction du diamètre pour 5 longueurs')
plt.xlabel('Diamètre (m)')
plt.ylabel('Chaleur dissipée (K)')
plt.legend([f'L = {aL[i]} m' if i < 5 else f'q = 10 W' for i in range(6)])
plt.savefig('chaleur.png',dpi=300)
plt.show()

print(f'D minimum = {D_min:.4f} m et L minimum = {L_min} m, donnant un V minimum de {(L_min*D_min**2)/4:.4} m^3')

# Correction
pytest.main(['-q', '--tb=long', 'ailette_corr.py'])
