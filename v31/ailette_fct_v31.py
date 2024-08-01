# Importation des modules
import numpy as np


# Assignation des paramètres
class parametres():
    
        L = 0.2      # [m] Longueur
        R = 0.2       # [m] Rayon
        k = 10       # [W/m*K] Conductivité thermique
        T_a = 0 + 273.15     # [K] Température de l'air ambiant
        T_w = 100 + 273.15     # [K] Température de la base
        
        nr = 50        #Nombre de points en r
        nz = 50        #Nombre de points en z
prm=parametres()        


def position(Z,R,nz,nr):
    """ Fonction générant deux matrices de discrétisation de l'espace

    Entrées:
        - R : Bornes du domaine en r, r = [r_min, r_max]
        - Z : Bornes du domaine en z, z = [z_min, z_max]
        - nr : Discrétisation de l'espace en x (nombre de points)
        - nz : Discrétisation de l'espace en y (nombre de points)

    Sorties (dans l'ordre énuméré ci-bas):
        - r : Matrice (array) de dimension (nr x nz) qui contient la position en r
        - z : Matrice (array) de dimension (nr x nz) qui contient la position en z
            * Exemple d'une matrice position :
            * Si r = [-1, 1] et z = [0, 1]
            * Avec nr = 3 et nz = 3
                r = [-1    0    1]
                    [-1    0    1]
                    [-1    0    1]

                z = [1    1    1  ]
                    [0.5  0.5  0.5]
                    [0    0    0  ]
    """

    # Fonction à écrire
    z=np.zeros([nz,nr])
    r = np.zeros([nz,nr])
    
    Z = [0,prm.L]
    R = [0,0.2]
    
    Lz=(abs(Z[1])-abs(Z[0]))/(nz-1)
    Lr=(abs(R[0])+abs(R[1]))/ (nr-1)

    for i in range(0,nr):
        z[i,0] = Z[0]

    for i in range(0,nz):
        r[0,i] = R[1]

    for i in range(0,nz):
        for j in range(1,nr):
            r[j,i]= r[j-1,i] - Lr
        
    for i in range(0,nr):
        for j in range(1,nz):
            z[i,j]= z[i,j-1] + Lz
    
    #for j in range(nz):
    #    z[:,j]=Z[0] + j*Lz
        
    #for i in range(nr):
    #    r[i,:]=R[1] - i*Lr

    return z,r
    
#TESTTTTTTTTTTTT

def mdf_CF_isolation(z,r,nz,nr,Bi): 
    """ Fonction assemblant la matrice A et le vecteur b
    Entrées:
        - r : Bornes du domaine en r, R = [r_min, r_max]
        - z : Bornes du domaine en z, Z = [z_min, z_max]
        - nr : Discrétisation de l'espace en r (nombre de points)
        - nz : Discrétisation de l'espace en z (nombre de points)
        - Bi : Nombre de Biot
    Sorties (dans l'ordre énuméré ci-bas):
        - A : Matrice (array)
        - b : Vecteur (array)
    """
    L = prm.L

    T_a = prm.T_a
    T_w = prm.T_w
    nr=prm.nr
    nz=prm.nz

    # Fonction à écrire
    def get_k(i,j):
        k=i*nr+j
        return k
    
    k=get_k
    R = [0,1]
    Z = [0,prm.L]
    Lr=abs((R[0]-R[1]) / (nr-1))
    Lr=R[1]-R[0] 
    Lz=Z[1]-Z[0] 
    dr=Lr/(nr-1) 
    dz=Lz/(nz-1) 
    
    #z,r=position(z,r,nz,nr)  
    N=nr*nz

    A=np.zeros([N,N]) 
    b=np.zeros(N) 

    #Au coeur du domaine OK
    for i in range(1,nz-1): 
        for j in range(1,nr-1):          
            k = i * nr + j
            A[k,k]= -((2/(dr**2))+(2/(dz**2)))
            A[k,k-1]= (1/(dr**2)) - (1/(r[i,j]*2*dr))
            A[k, k+1]= (1/(dr**2))+(1/(dr*r[i,j]*2))
            A[k,k-nr]= (1/(dz**2)) 
            A[k,k+nr]= (1/(dz**2))
            b[k]=0


    #Bout de l'ailette (condition -isolée- à droite)
    i=nz-1 #derniere colonne
    for j in range(0,nr-1): #parcours verticalement
        #methode gear arriere
        k=i*nr+j
        A[k,k]=3
        A[k,k-nr]=-4
        A[k,k-2*nr]= 1
        b[k] = 0

  
    #Base de l'ailette (condition à gauche) OK
    i=0 #premiere colonne
    for j in range(0,nr):
        k=i*nr+j
        A[k,k]=1
        b[k]=T_w
    
    #Surface ext. (condition en haut)
    j=0
    for i in range(0,nz):
        k=i*nr+j
        A[k,k]=3/(2*dr) + (Bi*prm.k/(prm.R*2))/prm.k #noter que h=(Bi*prm.k/(prm.R*2))
        A[k,k+1]= -4/(2*dr)
        A[k,k+2]= 1/(2*dr)
        b[k] = ((Bi*prm.k/(prm.R*2)))*T_a/prm.k 

    #Condition de symétrie en r=0 
    j=nr-1
    for i in range(0,nz):
        k=i*nr+j
        A[k,k]=3
        A[k,k-nr]=-4
        A[k,k-2*nr]= 1
        b[k] = 0
   
    return A,b


def mdf_CF_convection(z,r,nr,nz,Bi): 
    """ Fonction assemblant la matrice A et le vecteur b
    Entrées:
        - r : Bornes du domaine en r, R = [r_min, r_max]
        - z : Bornes du domaine en z, Z = [z_min, z_max]
        - nr : Discrétisation de l'espace en r (nombre de points)
        - nz : Discrétisation de l'espace en z (nombre de points)
        - Bi : Nombre de Biot
    Sorties (dans l'ordre énuméré ci-bas):
        - A : Matrice (array)
        - b : Vecteur (array)
    """
    L = prm.L

    T_a = prm.T_a
    T_w = prm.T_w
    nr=prm.nr
    nz=prm.nz

    # Fonction à écrire
    def get_k(i,j):
        k=i*nr+j
        return k
    
    k=get_k
    R = [0,1]
    Z = [0,prm.L]
    Lr=abs((R[0]-R[1]) / (nr-1))
    Lr=R[1]-R[0] 
    Lz=Z[1]-Z[0] 
    dr=Lr/(nr-1) 
    dz=Lz/(nz-1) 
    
    #z,r=position(z,r,nz,nr)  
    N=nr*nz

    A=np.zeros([N,N]) 
    b=np.zeros(N) 

    #Au coeur du domaine OK
    for i in range(1,nz-1): 
        for j in range(1,nr-1):          
            k = i * nr + j
            A[k,k]= -((2/(dr**2))+(2/(dz**2)))
            A[k,k-1]= (1/(dr**2)) - (1/(r[i,j]*2*dr))
            A[k, k+1]= (1/(dr**2))+(1/(dr*r[i,j]*2))
            A[k,k-nr]= (1/(dz**2)) 
            A[k,k+nr]= (1/(dz**2))
            b[k]=0


    #Bout de l'ailette (condition - de convection - à droite)
    i=nz-1 #derniere colonne
    for j in range(0,nr-1): #parcours verticalement
        #methode gear arriere
        k=i*nr+j
        A[k,k]=3
        A[k,k-nr]=-4
        A[k,k-2*nr]= 1
        b[k] = (Bi*prm.k/(prm.R*2))/prm.k #noter que h=(Bi*prm.k/(prm.R*2))

 
  
    #Base de l'ailette (condition à gauche) OK
    i=0 #premiere colonne
    for j in range(0,nr):
        k=i*nr+j
        A[k,k]=1
        b[k]=T_w
    
    #Surface ext. (condition en haut)
    j=0
    for i in range(0,nz):
        k=i*nr+j
        A[k,k]=3/(2*dr) + (Bi*prm.k/(prm.R*2))/prm.k #noter que h=(Bi*prm.k/(prm.R*2))
        A[k,k+1]= -4/(2*dr)
        A[k,k+2]= 1/(2*dr)
        b[k] = ((Bi*prm.k/(prm.R*2)))*T_a/prm.k 

    #Condition de symétrie en r=0 
    j=nr-1
    for i in range(0,nz):
        k=i*nr+j
        A[k,k]=3
        A[k,k-nr]=-4
        A[k,k-2*nr]= 1
        b[k] = 0
   
    return A,b

def inte_eq9(T, z, prm):
    """Fonction qui intègre le flux de chaleur sur une section transversale de l'ailette.

    Entrées:
        - T : Vecteur comprenant les températures en Kelvin sur la longueur de l'ailette
              pour une combinaison géométrique donnée
        - z : Matrice comprenant les positions en r et z sur le domaine
        - prm : Objet class parametres()
            - k : Conductivité thermique
            - h : Coefficient de convection
            - R : Rayon

    Sortie:
        - Valeur numérique de l'intégrale résultante (perte en W)
    """
    k = prm.k
    R = prm.R

    # Fonction à intégrer : -k * (∂T/∂z) * r
    dz_dr = np.gradient(z[:, 1], z[:, 0])  # Calcul du gradient dz/dr
    fonc = lambda z: -k * dz_dr * R

    delta_z = z[1, 1] - z[0, 1]
    I = (delta_z / 2) * np.sum(fonc(z[1:, :]) + fonc(z[:-1, :]))

    return I

def inte_eq10(T,z,prm,Bi):
    """Fonction qui intègre le flux de chaleur sur le contour de l'ailette.

    Entrées:
        - T : Vecteur comprenant les températures  en Kelvin sur la longueur de l'ailette
                pour une combinaison géométrique donnée
        - z : Vecteur comprenant les points sur la longueur en mètre
        - prm : Objet class parametres()
            - k : Conductivité thermique
            - T_a : Température ambiante
            - T_w : Température du mur
            - R : Rayon

    Sortie:
        - Valeur numérique de l'intégrale résultante (perte en W)
    """
    T_a = prm.T_a
    R = prm.R
    k = prm.k
    T_w = prm.T_w
    
    h = (Bi * k)/(2*R)
    fonc = lambda T : 2*h*np.pi*R*(T- T_a)   
    
    delta = z[1]-z[0]
    I = (delta/2) * sum(fonc(T[1:])+fonc(T[:-1]))

    return I
