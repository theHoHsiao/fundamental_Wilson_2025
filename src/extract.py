import matplotlib.pyplot as plt
import numpy as np
from scipy import linalg

from . import fitting
from .bootstrap import BootstrapSampleSet, BOOTSTRAP_SAMPLE_COUNT
from .plots_common import plot_mass_eff_cosh


def extract_meson_mass(C_tmp, plateau_start, plateau_end):
    E_fit, A_fit, chisquare = fitting.fit_cosh_bootstrap(
        C_tmp, plateau_start, plateau_end
    )

    return E_fit, A_fit, round(chisquare, 2)


def meson_decay_constant(Css, Csp, plateau_start, plateau_end):
    # load the ensamble info
    lattice_t = np.shape(Css.mean)[1]

    E_fit, A_fit, chisquare = fitting.fit_coshsinh_simultaneous(
        Css, Csp, plateau_start, plateau_end, lattice_t
    )

    return E_fit, A_fit, round(chisquare, 2)


def gevp_fixT(Cmat_mean, Cmat, t0, ti, tf):
    Mshape = Cmat.shape

    Lambda_n_sample = np.zeros(shape=(Mshape[0], Mshape[1], Mshape[2]))
    Lambda_n_mean = np.zeros(shape=(1, Mshape[1], Mshape[2]))

    T_dot = np.arange(ti, tf, 1, dtype=int)
    for t in T_dot:
        for N in range(Mshape[0]):
            value, vector = linalg.eig(
                Cmat[N, t], Cmat[N, t0], overwrite_a=True, overwrite_b=True
            )

            for n in range(Mshape[2]):
                Lambda_n_sample[N, t, n] = np.real(value[n])

        #### mean set ####
        value_m, vector_m = linalg.eig(
            Cmat_mean[0, t], Cmat_mean[0, t0], overwrite_a=True, overwrite_b=True
        )

        for n in range(Mshape[2]):
            Lambda_n_mean[0, t, n] = np.real(value_m[n])

    Lambda_n_sample.sort(axis=2)
    Lambda_n_mean.sort(axis=2)

    samples = Lambda_n_sample[:, :, ::-1]
    mean = Lambda_n_mean[:, :, ::-1]
    eigenvalues = []
    for n in range(Mshape[2]):
        eigenvalues.append(BootstrapSampleSet(mean[:, :, n], samples[:, :, n]))

    return eigenvalues


def extract_energy_states(eigenvalues, args):

    masses=[]
    chisquares_dof=[]

    for n in range(len(eigenvalues)):
        eigenvalue_n = eigenvalues[n]
        plateau_start = getattr(args, f"E{n}_plateau_start")
        plateau_end = getattr(args, f"E{n}_plateau_end")
        
        if plateau_start == 0 or plateau_end == 0:
            E_fit = BootstrapSampleSet(np.nan, np.nan * np.zeros(BOOTSTRAP_SAMPLE_COUNT))
            chisquare = np.nan

        else:
            try:
                E_fit, A_fit, chisquare = fitting.fit_cosh_bootstrap(eigenvalue_n, plateau_start, plateau_end)
            except:
                print(f"Error fitting {n}th energy state {plateau_start} {plateau_end}")
                fig, ax = plt.subplots(layout="constrained")
                plot_mass_eff_cosh(ax, eigenvalue_n, 2, args.Nt/2 + 1, "$E_"+f"{n}"+"$")
                plt.show()


        masses.append(E_fit)
        chisquares_dof.append(chisquare)
    
    return masses, chisquares_dof