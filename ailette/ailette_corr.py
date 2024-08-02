import numpy as np
try:
    from ailette_fct import *
except:
    print("ERREUR! Il y a une erreur fatale dans le fichier ailette_fct.py")


class prm():
    k = 180
    T_a = 293
    T_w = 373
    h = 120
    N = 5
    L = 0.2
    D = 0.07


class Test:

    def test_mdf(self):
        T = np.array([373.,355.7094,344.3911,337.9672,335.8259],dtype=float)
        z = np.array([0.,0.05,0.1,0.15,0.2],dtype=float)
        rep_mdf = mdf(prm())
        err_mdf = abs(np.asarray(rep_mdf[0]) - np.array([373.0000000,355.7093954,344.3911143,
                                                             337.9672250,335.8259286]))
        err_z = abs(np.asarray(rep_mdf[1] - np.array([0.,0.05,0.1,0.15,0.2])))

        assert (all(err_mdf < 1e-06))
        assert (all(err_z < 1e-06))

    def test_inte(self):
        T = np.array([373.,355.7094,344.3911,337.9672,335.8259],dtype=float)
        z = np.array([0.,0.05,0.1,0.15,0.2],dtype=float)
        rep_inte = inte(T,z,prm())

        print(rep_inte)

        assert (abs(rep_inte - 290.91736392545505) < 1e-04)

