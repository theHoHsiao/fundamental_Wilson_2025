#!/usr/bin/env python3
import numpy as np
import re


def get_ensemble(
    ensembles, beta=None, mAS=None, Nt=None, Ns=None, num_source=1, epsilon=None
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
    if len(candidate_ensembles) != num_source:
        raise ValueError("Did not uniquely identify one ensemble.")
    elif len(candidate_ensembles) == 0:
        raise ValueError("No ensembles found.")
    else:
        return candidate_ensembles


def filter_configurations(
    ensemble, min_trajectory=None, max_trajectory=None, trajectory_step=1
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

    return filtered_indices
