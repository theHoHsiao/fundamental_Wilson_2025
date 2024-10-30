#!/usr/bin/env python3

import re
import numpy as np

from .bootstrap import get_rng, sample_bootstrap_1d, BootstrapSampleSet
from . import extract


def C_R(ch):
    return {
        "v": -20.57,
        "av": -15.82,
        "ps": -15.82,
    }.get(ch, ch)


def get_correlators(
    ensembles, beta=None, mAS=None, Nt=None, Ns=None, num_source=None, epsilon=None
):
    candidate_ensembles = []
    for ensemble in ensembles.values():
        if beta is not None and ensemble.get("beta", {(): None})[()] != beta:
            continue
        if mAS is not None and (
            len(masses := ensemble.get("quarkmasses", [])) != 1 or masses[0] != mAS
        ):
            continue
        if Nt is not None and ensemble.get("lattice", [None])[0] != Nt:
            continue
        if Ns is not None and tuple(ensemble.get("lattice", [None])[-3:]) != (
            Ns,
            Ns,
            Ns,
        ):
            continue
        if epsilon is not None and ensemble.get("Wuppertal_eps_anti", [])[0] != epsilon:
            continue
        candidate_ensembles.append(ensemble)
    if num_source is not None and len(candidate_ensembles) != num_source:
        raise ValueError("Did not uniquely identify one ensemble.")
    elif len(candidate_ensembles) == 0:
        raise ValueError("No ensembles found.")
    else:
        return candidate_ensembles


def get_correlator_samples(
    ensemble,
    measurement,
    min_trajectory=None,
    max_trajectory=None,
    trajectory_step=1,
):
    indices = np.asarray(
        [
            int(re.match(".*n([0-9]+)$", filename.decode()).groups()[0])
            for filename in ensemble["configurations"]
        ]
    )
    filtered_indices = (
        ((indices >= min_trajectory) if min_trajectory is not None else True)
        & ((indices <= max_trajectory) if max_trajectory is not None else True)
        & ((indices - (min_trajectory or 0)) % trajectory_step == 0)
    )

    C = ensemble[measurement][:, filtered_indices]

    return BootstrapSampleSet(
        C.mean(axis=1), sample_bootstrap_1d(C.T, get_rng(ensemble.name))
    )


def channel_tags(ch):
    return {
        "ps": ["g5"],
        "v": ["g1", "g2", "g3"],
        "t": ["g0g1", "g0g2", "g0g3"],
        "av": ["g5g1", "g5g2", "g5g3"],
        "at": ["g0g5g1", "g0g5g2", "g0g5g3"],
        "s": ["id"],
    }.get(ch, ch)


def fold_correlators(C):
    return (C + np.roll(np.flip(C, axis=1), 1, axis=1)) / 2


def ps_extraction(ensemble, args):
    corr_aa = get_correlator_samples(
        ensemble,
        "g5",
        args.min_trajectory,
        args.max_trajectory,
        args.trajectory_step,
    )

    aa_mean = np.zeros(shape=(1, args.Nt))
    aa_mean[0] = corr_aa.mean * args.Ns**3
    C_aa = BootstrapSampleSet(aa_mean, corr_aa.samples * args.Ns**3)

    corr_ab = get_correlator_samples(
        ensemble,
        "g5_g0g5_re",
        args.min_trajectory,
        args.max_trajectory,
        args.trajectory_step,
    )

    ab_mean = np.zeros(shape=(1, args.Nt))
    ab_mean[0] = corr_ab.mean * args.Ns**3
    C_ab = BootstrapSampleSet(ab_mean, corr_ab.samples * args.Ns**3)

    m_tmp, a_tmp, chi2 = extract.meson_decay_sample(
        C_aa, C_ab, args.plateau_start, args.plateau_end
    )
    return m_tmp, a_tmp, chi2


def ch_extraction(ensemble, args):
    CHs = channel_tags(args.channel)

    tmp_bin = []
    tmp_bin_mean = []
    for j in range(len(CHs)):
        tmp_bin.append(
            get_correlator_samples(
                ensemble,
                CHs[j],
                args.min_trajectory,
                args.max_trajectory,
                args.trajectory_step,
            ).samples
            * args.Ns**3
        )
        tmp_bin_mean.append(
            get_correlator_samples(
                ensemble,
                CHs[j],
                args.min_trajectory,
                args.max_trajectory,
                args.trajectory_step,
            ).mean
            * args.Ns**3
        )

    mean = np.zeros(shape=(1, args.Nt))
    mean[0] = np.array(tmp_bin_mean).mean(axis=0)

    corr = BootstrapSampleSet(mean, np.array(tmp_bin).mean(axis=0))

    m_tmp, a_tmp, chi2 = extract.meson_mass_sample(
        corr, args.plateau_start, args.plateau_end
    )

    return m_tmp, a_tmp, chi2
