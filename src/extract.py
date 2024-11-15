from . import fitting
import numpy as np
from scipy import linalg
from .bootstrap import BootstrapSampleSet


def fold_correlators(C):
    return (C + np.roll(np.flip(C, axis=1), 1, axis=1)) / 2


def fold_correlators_cross(C):
    C_fold = (C - np.roll(np.flip(C, axis=1), 1, axis=1)) / 2

    C_fold[:, 0] = C[:, 0]

    return C_fold


def meson_mass_sample(C_tmp, plateau_start, plateau_end):
    E_mean, A_mean, X2, E_samples, A_samples = fitting.fit_cosh_booerr(
        C_tmp, plateau_start, plateau_end
    )

    E_fit = BootstrapSampleSet(E_mean, E_samples)
    A_fit = BootstrapSampleSet(A_mean / np.sqrt(E_mean), A_samples / np.sqrt(E_samples))

    return E_fit, A_fit, round(X2, 2)


def meson_decay_sample(Css, Csp, plateau_start, plateau_end):
    # load the ensamble info
    GLB_T = np.shape(Css.mean)[1]

    (E_mean, A_mean, X2, E_samples, A_samples) = fitting.fit_cosh_simultaneous(
        Css, Csp, plateau_start, plateau_end, GLB_T
    )

    E_fit = BootstrapSampleSet(E_mean, E_samples)
    A_fit = BootstrapSampleSet(A_mean / np.sqrt(E_mean), A_samples / np.sqrt(E_samples))

    return (
        E_fit,
        A_fit,
        round(X2, 2),
    )


def GEVP_fixT(Cmat_mean, Cmat, t0, ti, tf):
    Mshape = Cmat.shape

    Lambda_n = np.zeros(shape=(Mshape[0], Mshape[1], Mshape[2]))
    Lambda_n_mean = np.zeros(shape=(1, Mshape[1], Mshape[2]))

    # Vector_n = np.zeros(shape=Mshape, dtype=complex)

    T_dot = np.arange(ti, tf, 1, dtype=int)
    for t in T_dot:
        for N in range(Mshape[0]):
            value, vector = linalg.eig(
                Cmat[N, t], Cmat[N, t0], overwrite_a=True, overwrite_b=True
            )

            # Vector_n[N, t] = vector

            for n in range(Mshape[2]):
                Lambda_n[N, t, n] = np.real(value[n])

        #### mean set ####
        value_m, vector_m = linalg.eig(
            Cmat_mean[0, t], Cmat_mean[0, t0], overwrite_a=True, overwrite_b=True
        )

        for n in range(Mshape[2]):
            Lambda_n_mean[0, t, n] = np.real(value_m[n])

    Lambda_n.sort(axis=2)
    Lambda_n_mean.sort(axis=2)

    samples = Lambda_n[:, :, ::-1]
    mean = Lambda_n_mean[:, :, ::-1]
    eigenvalues = []
    for n in range(Mshape[2]):
        eigenvalues.append(BootstrapSampleSet(mean[:, :, n], samples[:, :, n]))

    return eigenvalues  # , Vector_n
