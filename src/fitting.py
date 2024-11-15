import numpy as np
from scipy.optimize import curve_fit
import corrfitter as cf
import gvar as gv
from scipy.optimize import minimize


def make_models(T, tmin, tmax, tp):
    """Create corrfitter model for G(t)."""
    return [cf.Corr2(datatag="Gab", tp=tp, tmin=tmin, tmax=tmax, a="a", b="a", dE="dE")]


def make_prior(N):
    prior = gv.BufferDict()
    # setting the sdev of the prioir to infinity amounts to turning off the prior contribution to chi2
    prior["log(a)"] = gv.log(gv.gvar(N * [0.1], N * [np.inf]))
    prior["log(dE)"] = gv.log(gv.gvar(N * [0.1], N * [np.inf]))
    return prior


def first_fit_parameters(fit):
    p = fit.p
    E = np.cumsum(p["dE"])
    a = p["a"]
    chi2 = fit.chi2
    dof = fit.dof
    return E, a, chi2, dof


def bootstrap_fit(fitter, dset, T, tmin, tmax, tp):
    gv.ranseed(12)

    pdatalist = (
        cf.process_dataset(ds, make_models(T, tmin, tmax, tp))
        for ds in gv.dataset.bootstrap_iter(dset, n=10)
    )
    bs = gv.dataset.Dataset()
    for bsfit in fitter.bootstrapped_fit_iter(pdatalist=pdatalist):
        bs.append(
            E=np.cumsum(bsfit.pmean["dE"][:2]),
            a=bsfit.pmean["a"][:2],
        )
    bs = gv.dataset.avg_data(bs, bstrap=True)
    E = bs["E"]
    a = bs["a"]
    print("{:2}  {:15}  {:15}".format("E", E[0], E[1]))
    print("{:2}  {:15}  {:15}".format("a", a[0], a[1]))


def fit_correlator_without_bootstrap(
    data_corr, t_lattice, tmin, tmax, Nmax, tp, p0, plotting=False, printing=False
):
    t_lattice = abs(t_lattice)

    fitter = cf.CorrFitter(models=make_models(t_lattice, tmin, tmax, tp))

    for N in range(1, Nmax + 1):
        prior = make_prior(N)
        fit = fitter.lsqfit(data=data_corr, prior=prior, p0=p0)
        p0 = fit.pmean

        if printing:
            print("nterm =", N, 30 * "=")
            print(fit)

    E, a, chi2, dof = first_fit_parameters(fit)
    if plotting:
        fit.show_plots(view="log")

    return E, a, chi2, dof


def fit_exp_std(C_boot, plateau_start, plateau_end):
    # This function fits the mean correlators with a exp function
    # the error is estimated by standard deviation with a covariance matrix

    cov = np.cov(C_boot[0:-1].T)

    if np.isnan(cov.sum()):
        print("cov contain nan")

    if np.isnan(cov.T.sum()):
        print("cov contain nan")

    correlator_set = dict(Gab=gv.gvar(C_boot[-1], cov))

    E, a, chi2, dof = fit_correlator_without_bootstrap(
        correlator_set,
        0,
        plateau_start,
        plateau_end,
        1,
        None,
        plotting=False,
        printing=True,
    )

    return gv.mean(E[0]), gv.sdev(E[0]), chi2 / dof


def fit_cosh_std(C_boot, plateau_start, plateau_end, lattice_t):
    # This function fits the mean correlators with a cosh function
    # the error is estimated by standard deviation with a covariance matrix

    def func(t, a, M):
        return a * a * M * (np.exp(-M * t) + np.exp(-M * (lattice_t - t))) / 2

    x0, pcov = curve_fit(
        func,
        np.arange(plateau_start, plateau_end),
        C_boot[-1, plateau_start:plateau_end],
    )

    p0 = dict(
        {"log(a)": np.array([np.log(abs(x0[0]))]), "log(dE)": np.array([np.log(x0[1])])}
    )

    print(
        "A* ( exp[-mt] + exp[-m(T-t)] ) fitting time region ",
        plateau_start,
        "to",
        plateau_end,
        ": ",
    )

    cov = np.cov(C_boot[0:-1].T)

    if np.isnan(cov.sum()):
        print("cov contain nan")

    if np.isnan(cov.T.sum()):
        print("cov contain nan")

    correlator_set = dict(Gab=gv.gvar(C_boot[-1], cov))

    E, a, chi2, dof = fit_correlator_without_bootstrap(
        correlator_set,
        0,
        plateau_start,
        plateau_end,
        1,
        lattice_t,
        p0,
        plotting=False,
        printing=True,
    )

    return gv.mean(E[0]), gv.sdev(E[0]), chi2 / dof


def fit_cosh_booerr(C, plateau_start, plateau_end):
    # This function fits the correlators with a cosh function

    C_boot = C.samples

    num_sample = C_boot.shape[0]
    lattice_t = C_boot.shape[1]

    def func(t, a, M):
        return a * a * M * (np.exp(-M * t) + np.exp(-M * (lattice_t - t))) / 2

    x0, pcov = curve_fit(
        func,
        np.arange(plateau_start, plateau_end),
        C.mean[0, plateau_start:plateau_end],
    )

    p0 = dict(
        {"log(a)": np.array([np.log(abs(x0[0]))]), "log(dE)": np.array([np.log(x0[1])])}
    )

    E_sample = np.zeros(num_sample)
    a_sample = np.zeros(num_sample)
    chi2_dof = np.zeros(num_sample)

    cov = np.cov(C_boot.T)

    for n in range(num_sample):
        correlator_set = dict(Gab=gv.gvar(C_boot[n], cov))

        E, a, chi2, dof = fit_correlator_without_bootstrap(
            correlator_set,
            0,
            plateau_start,
            plateau_end,
            1,
            lattice_t,
            p0,
            plotting=False,
            printing=False,
        )
        E_sample[n] = gv.mean(E[0])
        a_sample[n] = gv.mean(a[0])
        chi2_dof[n] = chi2 / dof

    correlator_set = dict(Gab=gv.gvar(C.mean[0], cov))
    E_mean, a_mean, chi2, dof = fit_correlator_without_bootstrap(
        correlator_set,
        0,
        plateau_start,
        plateau_end,
        1,
        lattice_t,
        p0,
        plotting=False,
        printing=False,
    )

    return (
        gv.mean(E_mean[0]),
        gv.mean(a_mean[0]),
        chi2 / dof,
        E_sample,
        a_sample,
    )


def fit_exp_booerr(C, plateau_start, plateau_end):
    # This function fits the correlators with a exp function

    C_boot = C.samples

    def func(t, a, M):
        return a * a * M * (np.exp(-M * t)) / 2

    x0, pcov = curve_fit(
        func,
        np.arange(plateau_start, plateau_end),
        C.mean[0, plateau_start:plateau_end],
    )

    p0 = dict(
        {"log(a)": np.array([np.log(abs(x0[0]))]), "log(dE)": np.array([np.log(x0[1])])}
    )

    # load the ensamble info
    num_sample = C_boot.shape[0]

    E_sample = np.zeros(num_sample)
    a_sample = np.zeros(num_sample)
    chi2_dof = np.zeros(num_sample)

    cov = np.cov(C_boot.T)

    for n in range(num_sample):
        correlator_set = dict(Gab=gv.gvar(C_boot[n], cov))

        E, a, chi2, dof = fit_correlator_without_bootstrap(
            correlator_set,
            0,
            plateau_start,
            plateau_end,
            1,
            None,
            p0,
            plotting=False,
            printing=False,
        )
        E_sample[n] = gv.mean(E[0])
        a_sample[n] = gv.mean(a[0])
        chi2_dof[n] = chi2 / dof

    correlator_set = dict(Gab=gv.gvar(C.mean[0], cov))
    E_mean, a_mean, chi2, dof = fit_correlator_without_bootstrap(
        correlator_set,
        0,
        plateau_start,
        plateau_end,
        1,
        None,
        p0,
        plotting=False,
        printing=False,
    )

    return (
        gv.mean(E_mean[0]),
        gv.mean(a_mean[0]),
        chi2 / dof,
        E_sample,
        a_sample,
    )


def sim_model(T, tmin, tmax, tp):
    """Create corrfitter model for Gaa(t) and Gab(t)."""
    model = [
        cf.Corr2(datatag="Gaa", tp=tp, tmin=tmin, tmax=tmax, a="a", b="a", dE="dE"),
        cf.Corr2(datatag="Gab", tp=-tp, tmin=tmin, tmax=tmax, a="a", b="b", dE="dE"),
    ]
    return model


def sim_prior(N):
    prior = gv.BufferDict()
    # setting the sdev of the prioir to infinity amounts to turning off the prior contribution to chi2
    prior["log(a)"] = gv.log(gv.gvar(N * [0.1], N * [np.inf]))
    prior["log(b)"] = gv.log(gv.gvar(N * [0.1], N * [np.inf]))
    prior["log(dE)"] = gv.log(gv.gvar(N * [0.1], N * [np.inf]))
    return prior


def sim_fit_parameters(fit):
    p = fit.p
    E = np.cumsum(p["dE"])
    a = p["a"]
    b = p["b"]
    chi2 = fit.chi2
    dof = fit.dof
    return E, a, b, chi2, dof


def fit_correlator_simultaneous(
    data_corrs, T, tmin, tmax, Nmax, tp, p0, plotting=False, printing=False
):
    T = abs(T)

    fitter = cf.CorrFitter(models=sim_model(T, tmin, tmax, tp))
    for N in range(1, Nmax + 1):
        prior = sim_prior(N)
        fit = fitter.lsqfit(data=data_corrs, prior=prior, p0=p0)
        p0 = fit.pmean

        if printing:
            print("nterm =", N, 30 * "=")
            print(fit)

    E, a, b, chi2, dof = sim_fit_parameters(fit)

    if plotting:
        # fit.show_plots(view='ratio')
        fit.show_plots(view="log")

    return E, a, b, chi2, dof


def fit_cosh_simultaneous(Corr_ss, Corr_sp, plateau_start, plateau_end, lattice_t):
    # This function fits the correlators with cosh x sinh functions

    x0 = sim_coshsinh_fit(
        Corr_sp.mean[0], Corr_ss.mean[0], lattice_t, plateau_start, plateau_end
    )

    p0 = dict(
        {
            "log(a)": np.array([np.log(abs(x0[0]))]),
            "log(b)": np.array([np.log(x0[1])]),
            "log(dE)": np.array([np.log(x0[2])]),
        }
    )

    # load the ensamble info
    Css = Corr_ss.samples
    Csp = Corr_sp.samples

    num_sample = Css.shape[0]

    E_sample = np.zeros(num_sample)
    a_sample = np.zeros(num_sample)
    b_sample = np.zeros(num_sample)
    chi2_dof = np.zeros(num_sample)

    cov_ss = np.cov(Css.T)
    cov_sp = np.cov(Csp.T)

    for n in range(num_sample):
        correlator_set = dict(Gab=gv.gvar(Csp[n], cov_sp), Gaa=gv.gvar(Css[n], cov_ss))

        E, a, b, chi2, dof = fit_correlator_simultaneous(
            correlator_set,
            0,
            plateau_start,
            plateau_end,
            1,
            lattice_t,
            p0,
            plotting=False,
            printing=False,
        )
        E_sample[n] = gv.mean(E[0])
        a_sample[n] = gv.mean(a[0])
        b_sample[n] = gv.mean(b[0])
        chi2_dof[n] = chi2 / dof

    correlator_set = dict(
        Gab=gv.gvar(Corr_sp.mean[0], cov_sp), Gaa=gv.gvar(Corr_ss.mean[0], cov_ss)
    )

    E_mean, a_mean, b_mean, chi2, dof = fit_correlator_simultaneous(
        correlator_set,
        0,
        plateau_start,
        plateau_end,
        1,
        lattice_t,
        p0,
        plotting=False,
        printing=False,
    )

    return gv.mean(E_mean[0]), gv.mean(b_mean[0]), chi2 / dof, E_sample, b_sample


def sim_coshsinh_fit(C1, C2, T, ti, tf):
    y1 = C1[ti:tf]
    y2 = C2[ti:tf]
    comboY = np.append(y1, y2)
    t_slice = np.arange(ti, tf)
    comboX = np.append(t_slice, t_slice)

    def func_sp(t, As, f, M):
        return As * f * (np.exp(-M * t) - np.exp(-M * (T - t))) / 2.0

    def func_ss(t, As, f, M):
        return As**2 * (np.exp(-M * t) + np.exp(-M * (T - t))) / (2.0 * M)

    def comboFunc(
        comboData, a, b, c
    ):  # single data set passed in, extract separate data
        extract1 = comboData[: len(y1)]  # first data
        extract2 = comboData[len(y2) :]  # second data

        result1 = func_sp(extract1, a, b, c)
        result2 = func_ss(extract2, a, b, c)

        return np.append(result1, result2)

    # curve fit the combined data to the combined function
    fittedParameters, pcov = curve_fit(comboFunc, comboX, comboY)

    return fittedParameters


def fit_exp_simultaneous(Css, Csp, plateau_start, plateau_end, lattice_t):
    # This function fits the correlators with a exp function

    x0, y_fit_1, y_fit_2 = sim_coshsinh_fit(
        Csp[-1], Css[-1], lattice_t, plateau_start, plateau_end
    )

    p0 = dict(
        {
            "log(a)": np.array([np.log(abs(x0[0]))]),
            "log(b)": np.array([np.log(x0[1])]),
            "log(dE)": np.array([np.log(x0[2])]),
        }
    )

    print("p0 =", x0)

    # load the ensamble info
    num_sample = Css.shape[0]

    print(
        "A* ( exp[-mt] ) fitting time region ", plateau_start, "to", plateau_end, ": "
    )

    E_sample = np.zeros(num_sample)
    a_sample = np.zeros(num_sample)
    b_sample = np.zeros(num_sample)
    chi2_dof = np.zeros(num_sample)

    cov_ss = np.cov(Css.T)
    cov_sp = np.cov(Csp.T)

    for n in range(num_sample):
        correlator_set = dict(Gab=gv.gvar(Csp[n], cov_sp), Gaa=gv.gvar(Css[n], cov_ss))

        pp = False

        if n == num_sample - 1:
            pp = True

        E, a, b, chi2, dof = fit_correlator_simultaneous(
            correlator_set,
            0,
            plateau_start,
            plateau_end,
            1,
            None,
            p0,
            plotting=False,
            printing=pp,
        )
        E_sample[n] = gv.mean(E[0])
        a_sample[n] = gv.mean(a[0])
        b_sample[n] = gv.mean(b[0])
        chi2_dof[n] = chi2 / dof

    E_err = E_sample[0:-1].std()
    b_err = b_sample[0:-1].std()

    return E_sample, E_err, chi2_dof.mean(), b_sample, b_err


def meson_M2(X, LAT_A, Y):
    num_sample = np.shape(X)[1]
    num_pars = 3

    print("meson M2 fitting... M = a**2 * ( 1 + b * m**2 ) + c * lat_a ")

    def func(V, a, b, c):
        m, lat_a = V
        return a * (1 + b * m) + c * lat_a

    def Cov(i, j):
        return np.mean((Y[i, :] - np.mean(Y[i, :])) * (Y[j, :] - np.mean(Y[j, :])))

    x0, pcov = curve_fit(func, (X[:, -1], LAT_A[:, -1]), Y[:, -1])

    size = np.shape(X)[0]

    M = np.asmatrix(np.zeros(shape=(size, size)))

    for a in range(size):
        M[a, a] = Cov(a, a)

    M_I = M.I

    def X2_boot_const(pars):
        def Cf_vector():
            V = np.zeros(size)
            for i in range(size):
                V[i] = Y[i, N] - func((X[i, N], LAT_A[i, -1]), a, b, c)

            return V

        (a, b, c) = pars

        V = np.asmatrix(Cf_vector())

        return V * M_I * V.T

    def nor_X2(a, b, c):
        V = np.zeros(size)
        for i in range(size):
            V[i] = Y[i, -1] - func((X[i, -1], LAT_A[i, -1]), a, b, c)

        V = np.asmatrix(V)

        chisqr = (V * M_I * V.T)[0, 0]
        return chisqr

    fit_val = np.zeros(shape=(num_pars, num_sample))

    for N in range(num_sample):
        res = minimize(X2_boot_const, x0, method="Nelder-Mead", tol=10**-16)

        for n_pars in range(num_pars):
            fit_val[n_pars, N] = res.x[n_pars]

    X2 = nor_X2(fit_val[0].mean(), fit_val[1].mean(), fit_val[2].mean())

    return fit_val, X2 / (size - num_pars - 1)


def meson_beta(X, Y):
    num_sample = np.shape(X)[1]
    num_pars = 2

    print("meson M2 fitting... Y = A + BX^2 ")

    def func(m, a, b):
        return a + b * m

    def Cov(i, j):
        return np.mean((Y[i, :] - np.mean(Y[i, :])) * (Y[j, :] - np.mean(Y[j, :])))

    x0, pcov = curve_fit(func, X[:, -1], Y[:, -1])

    size = np.shape(X)[0]

    M = np.asmatrix(np.zeros(shape=(size, size)))

    for a in range(size):
        M[a, a] = Cov(a, a)

    M_I = M.I

    def X2_boot_const(pars):
        def Cf_vector():
            V = np.zeros(size)
            for i in range(size):
                V[i] = Y[i, N] - func(X[i, N], a, b)

            return V

        (a, b) = pars

        V = np.asmatrix(Cf_vector())

        return V * M_I * V.T

    def nor_X2(a, b):
        V = np.zeros(size)
        for i in range(size):
            V[i] = Y[i, -1] - func(X[i, -1], a, b)

        V = np.asmatrix(V)

        chisqr = (V * M_I * V.T)[0, 0]
        return chisqr

    fit_val = np.zeros(shape=(num_pars, num_sample))

    for N in range(num_sample):
        res = minimize(X2_boot_const, x0, method="Nelder-Mead", tol=10**-16)

        for n_pars in range(num_pars):
            fit_val[n_pars, N] = res.x[n_pars]

    X2 = nor_X2(fit_val[0].mean(), fit_val[1].mean())

    return fit_val, X2 / (size - num_pars)


def meson_beta_quad(X, Y):
    num_sample = np.shape(X)[1]
    num_pars = 3

    def func(m, a, b, c):
        return a + b * m + c * m**2

    def Cov(i, j):
        return np.mean((Y[i, :] - np.mean(Y[i, :])) * (Y[j, :] - np.mean(Y[j, :])))

    x0, pcov = curve_fit(func, X[:, -1], Y[:, -1])

    size = np.shape(X)[0]

    M = np.asmatrix(np.zeros(shape=(size, size)))

    for a in range(size):
        M[a, a] = Cov(a, a)

    M_I = M.I

    def X2_boot_const(pars):
        def Cf_vector():
            V = np.zeros(size)
            for i in range(size):
                V[i] = Y[i, N] - func(X[i, N], a, b, c)

            return V

        (a, b, c) = pars

        V = np.asmatrix(Cf_vector())

        return V * M_I * V.T

    def nor_X2(a, b, c):
        V = np.zeros(size)
        for i in range(size):
            V[i] = Y[i, -1] - func(X[i, -1], a, b, c)

        V = np.asmatrix(V)

        chisqr = (V * M_I * V.T)[0, 0]
        return chisqr

    fit_val = np.zeros(shape=(num_pars, num_sample))

    for N in range(num_sample):
        res = minimize(X2_boot_const, x0, method="Nelder-Mead", tol=10**-16)

        for n_pars in range(num_pars):
            fit_val[n_pars, N] = res.x[n_pars]

    X2 = nor_X2(fit_val[0].mean(), fit_val[1].mean(), fit_val[2].mean())

    return fit_val, X2 / (size - num_pars)
