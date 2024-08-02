# Importation des modules
import numpy as np

def mdf(prm):
    """Fonction qui calcule le profil de température le long de l'ailette

    Entrées:
        - prm : Objet class parametres()
            - L : Longueur
            - D : Diamètre
            - k : Conductivité thermique
            - T_a : Température ambiante
            - T_w : Température du mur
            - h : Coefficient de convection
            - N : Nombre de points utilisés pour la méthode

    Sorties (dans l'ordre énuméré ci-bas):
        - Vecteur (np.array) donnant la température tout au long de l'ailette en Kelvin
        - Vecteur (np.array) donnant la position tout au long de l'ailette (axe z) en mètre
    """

    L = prm.L
    D = prm.D
    k = prm.k
    T_a = prm.T_a
    T_w = prm.T_w
    h = prm.h
    N = prm.N
    
    z  = np.linspace(0,L,N)
    dz = L / ( N - 1 )
    b = np.zeros(N)
    A = np.zeros([N,N])
    
    A[0,0] = 1
    
    for i in range(1,N-1):
        A[i, i-1] = D * k
        A[i, i]   = - ( 2 * D * k ) - (4 * h * ( dz**2 ) )
        A[i, i+1] = D * k
        b[i] = - 4 * h * (dz**2) * T_a 
        
    A[N-1,N-1] = 3
    A[N-1,N-2] = -4
    A[N-1,N-3] = 1
    
    b[0] = T_w
    
    Tz = np.linalg.solve(A, b)

    return  Tz, z

def inte(T,z,prm):
    """Fonction qui intègre la convection sur la surface de l'ailette.

    Entrées:
        - T : Vecteur comprenant les températures  en Kelvin sur la longueur de l'ailette
                pour une combinaison géométrique donnée
        - z : Vecteur comprenant les points sur la longueur en mètre
        - prm : Objet class parametres()
            - k : Conductivité thermique
            - T_a : Température ambiante
            - T_w : Température du mur
            - h : Coefficient de convection
            - N : Nombre de points utilisés pour la méthode

    Sortie:
        - Valeur numérique de l'intégrale résultante (perte en W)
    """
    
    D = prm.D
    k = prm.k
    T_a = prm.T_a
    T_w = prm.T_w
    h = prm.h
    N = prm.N
    
    I=0
    for i in range(0,len(T)-1):
        B = h * np.pi * D * (T[i] - T_a)  
        b = h * np.pi * D * (T[i+1] - T_a) 
        dz = z[i+1] - z[i]
        I = I + (B + b) * dz /2

    return I