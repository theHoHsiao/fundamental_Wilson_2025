from . import fitting
import numpy as np
from scipy import linalg


def fold_correlators(C):
    return (C + np.roll(np.flip(C, axis=1), 1, axis=1)) / 2


def fold_correlators_cross(C):
    C_fold = (C - np.roll(np.flip(C, axis=1), 1, axis=1)) / 2

    C_fold[:, 0] = C[:, 0]

    return C_fold


def meson_mass_sample(C_tmp, TI, TF):
    C_fold = fold_correlators(C_tmp)

    E_fit, A_fit, X2 = fitting.fit_cosh_booerr(C_fold, TI, TF)

    return E_fit, A_fit / np.sqrt(E_fit), round(X2, 2)


def meson_decay_sample(Css, Csp, TI, TF):
    # load the ensamble info
    GLB_T = np.shape(Css)[1]

    (
        E_fit,
        b_fit,
        X2,
    ) = fitting.fit_cosh_simultaneous(Css, Csp, TI, TF, GLB_T)

    return (
        E_fit,
        b_fit / np.sqrt(E_fit),
        round(X2, 2),
    )


def GEVP_fixT(Cmat, t0, ti, tf):
    Mshape = Cmat.shape

    Lambda_n = np.zeros(shape=(Mshape[0], Mshape[1], Mshape[2]))
    Vector_n = np.zeros(shape=Mshape, dtype=complex)

    T_dot = np.arange(ti, tf, 1, dtype=int)
    for t in T_dot:
        for N in range(Mshape[0]):
            value, vector = linalg.eig(
                Cmat[N, t], Cmat[N, t0], overwrite_a=True, overwrite_b=True
            )

            Vector_n[N, t] = vector

            for n in range(Mshape[2]):
                Lambda_n[N, t, n] = np.real(value[n])

    Lambda_n.sort(axis=2)

    return Lambda_n[:, :, ::-1], Vector_n
