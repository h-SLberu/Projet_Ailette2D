# Importation des modules
import numpy as np
import matplotlib.pyplot as plt

try:
    from ailette_fct_v31         import *
except:
    pass

#------------------------------------------------------------------------------
# Code principal pour l'analyse des résultats
# Il faudra faire appel aux fonctions programmées dans ailette_fct.py afin de
# modéliser les pertes de chaleurs selon différentes combinaisons géométriques.
#------------------------------------------------------------------------------



# Assignation des paramètres
class parametres():
    
        L = 0.2      # [m] Longueur
        R = 0.2       # [m] Rayon
        k = 10       # [W/m*K] Conductivité thermique
        T_a = 0 + 273.15     # [K] Température de l'air ambiant
        T_w = 100 + 273.15     # [K] Température de la base

        
        
        
prm=parametres()        

Z = [0,prm.L]
R = [0,prm.R]

nr = 50
nz = 50

Bi = [0.1, 1, 10, 20, 100]

z,r = position(Z, R, nz, nr)


# #%% PREMIÈRE ANALYSE : Détermination du profil de température à r = 0 et à r = R



######################################## ICI ##########################


fig, axes = plt.subplots(nrows=1, ncols=len(Bi), figsize=(17,6))
fig.suptitle("Cartes de chaleur pour différents nombres de Biot selon la condition frontière d’isolation", fontsize=16)

for i, biot in enumerate(Bi):
    # Obtenez la matrice A et le vecteur b pour un Bi donné
    A, b = mdf_CF_isolation(z,r,nz,nr,biot)

    # Résolvez le système linéaire
    T_profile = np.linalg.solve(A, b)
    # Affichez chaque heatmap dans une sous-figure
    T_profile_reshaped = T_profile.reshape(nz, nr).transpose()
    im = axes[i].imshow(T_profile_reshaped, cmap='hot', interpolation='nearest', extent=[0, 0.2, 0, 0.2])
    axes[i].set_title(f"Biot = {biot}")
    axes[i].set_xlabel("Position en z")
    axes[i].set_ylabel("Position en r")
    cbar = fig.colorbar(im, ax=axes[i])
    cbar.set_label("Températures en kelvin")

plt.tight_layout()
plt.show()
###############################################################################################

# Résolution pour la condition frontière de convection
fig, axes = plt.subplots(nrows=1, ncols=len(Bi), figsize=(17,6))
fig.suptitle("Cartes de chaleur pour différents nombres de Biot selon la condition frontière de convection", fontsize=16)


for i, biot in enumerate(Bi):
    # Obtenez la matrice A et le vecteur b pour un Bi donné
    A_convection, b_convection = mdf_CF_convection(z, r, nr, nz, biot)

    # Résolvez le système linéaire
    T_profile = np.linalg.solve(A_convection, b_convection)

    # Affichez chaque heatmap dans une sous-figure
    im = axes[i].imshow(T_profile.reshape((nz, nr)), cmap='hot', interpolation='nearest', extent=[0, 0.2, 0, 0.2])
    axes[i].set_title(f"Biot = {biot}")
    axes[i].set_xlabel("Position en z")
    axes[i].set_ylabel("Position en r")
    cbar = fig.colorbar(im, ax=axes[i])
    cbar.set_label("Températures en kelvin")

plt.tight_layout()
plt.show()




#%% GÉNÉRATION D'UN GRAPHIQUE T(z) selon le nombre de Biot 

T_w = prm.T_w
T_a = prm.T_a
k = prm.k
L = prm.L
R = prm.R
z1= np.linspace(0,L,nz)


# # Obtain the matrix A and vector b
# A, b = mdf_CF_isolation(z, r, nr, nz, biot)
# #Test pour voir profil de T Bi=1
# A1, b1 = mdf_CF_isolation(z, r, nr, nz, 1)
# T_profile_Bi1=np.linalg.solve(A1, b1).reshape((nz, nr))
# print(T_profile_Bi1) 

# Solve the linear system
#T_profile = np.linalg.solve(A, b).reshape((nz, nr))
#T_profile_Bi1=np.linalg.solve(A1, b1).reshape((nz, nr))
#print(T_profile_Bi1) 

plt.figure(figsize=(8, 6))
for biot in Bi:
    # Résolution du système linéaire
    A, b = mdf_CF_isolation(z,r,nr,nz,biot)
    T = np.linalg.solve(A, b).reshape((nz, nr))
    
    
    # Sélectionner la température au centre pour r=0
    T_centre = T[:,0] #en r=0

    # Tracer la température en fonction de la position z
    plt.plot(z1, T_centre, label=f'Bi = {biot}')

plt.xlabel('Position z (m)')
plt.ylabel('Température (K)')

plt.title('Température en r=0 en selon z pour 5 nombres de Bi (condition isolée)')
plt.legend(('Bi = 0.1', 'Bi = 1', 'Bi = 10', 'Bi = 20' ,'Bi = 100'))
plt.grid(True)
plt.show()


plt.figure(figsize=(8, 6))
for biot in Bi:
    # Résolution du système linéaire
    A, b = mdf_CF_isolation(z,r,nr,nz,biot)
    T = np.linalg.solve(A, b).reshape((nz, nr))
    #print('les temperatures sont:',T)
    
    
    # Sélectionner la température au centre pour r=R
    T_ext = T[:, -1] #en r=R

    # Tracer la température en fonction de la position z
    plt.plot(z1, T_ext, label=f'Bi = {biot}')

plt.xlabel('Position z (m)')
plt.ylabel('Température (K)')

plt.title('Température en r=R en selon z pour 5 nombres de Bi (condition isolée)')
plt.legend(('Bi = 0.1', 'Bi = 1', 'Bi = 10', 'Bi = 20' ,'Bi = 100'))
plt.grid(True)
plt.show()
#print(A)

#%% GÉNÉRATION D'UN GRAPHIQUE T(z) pour la solution analytique

# Calcul de la solution analytique pour Bi < 0.1
T_w = prm.T_w
L = prm.L

z1 = np.linspace(0, L, nz)
Bi_analytique = [0.000001, 0.00001, 0.0001, 0.001, 0.01]

# Initialize a list to store the results for each Bi_analytique
results = []

for i in Bi_analytique:
    m = np.sqrt((2 * ((i * prm.k)/ (prm.R * 2))) / (prm.k * prm.R))
    T_analytique = prm.T_a + ((T_w - prm.T_a) * (np.cosh(m * (L - z1)) / np.cosh(m * L)))
    
    # Append the results to the list
    results.append(T_analytique)

# Tracer la solution analytique
plt.figure(figsize=(8, 6))

# Plot each solution from the list
for i, result in enumerate(results):
    plt.plot(z1, result, label=f'Bi = {Bi_analytique[i]}')

plt.xlabel('Position z')
plt.ylabel('Température(K)')
plt.title('Solution analytique de la température en fonction de la position z (Bi < 0.1)')
plt.legend()
plt.grid(True)
plt.show()



    
#%% GÉNÉRATION DES VALEURS DE FLUX ET DÉTERMINATION DE L'ERREUR COMMISE (AVEC RATIO ENTRE SIMULÉ ET ANALYTIQUE)

T_w = prm.T_w
T_a = prm.T_a
k = prm.k
L = prm.L
R = prm.R


# Calcul des flux de chaleur et des ratios pour chaque Bi

for biot in Bi:
    # Obtention des températures simulées et la solution analytique pour chaque Bi
    A, b = mdf_CF_convection(z,r,nr,nz,biot)
    temperature_simule = np.linalg.solve(A, b)
    temperature_analyt = T_a + (T_w - T_a) * np.cosh(m * (L - z)) / np.cosh(m * L)
    
    # Calcul des flux de chaleur
    q_simule = inte_eq10(temperature_simule, z, prm, biot) 
    q_analyt = inte_eq10(temperature_analyt, z, prm, biot)
    print(q_analyt)
    print(q_simule)

    
    ratios = [q_simule[i] / q_analyt[i] for i in range(len(Bi))]

# Impression des ratios correspondants aux 5 nombres de Bi
print("Ratios entre les flux de chaleur simulés et analytiques pour chaque Bi:")
for i, biot in enumerate(Bi):
	print(f"Bi = {biot}: Ratio = {ratios[i]}") 





