#!/usr/bin/env python3
from .dump import read_sample_files
from .fitting import meson_beta
from argparse import ArgumentParser, FileType
import numpy as np
from .dump import dump_dict, dump_samples
from uncertainties import ufloat


def get_args():
    parser = ArgumentParser()

    parser.add_argument(
        "data_filenames",
        nargs="+",
        metavar="sample_filename",
        help="Filenames of sample files used for extrapolation",
    )
    parser.add_argument(
        "--beta",
        type=float,
        choices=[6.6, 6.65, 6.7, 6.75, 6.8],
        default=None,
        help="Performing on the ensembles with the beta value",
    )
    parser.add_argument("--output_file_mean", type=FileType("w"), default="-")

    parser.add_argument(
        "--output_file_samples",
        type=FileType("w"),
        default=None,
        help="Where to output the bootstrap samples for fitting results",
    )

    return parser.parse_args()


def prepare_data(data, args):
    X = []
    Y = []
    bare_mass = []

    for datum in data:
        if "ps_mass_samples" not in datum:
            continue
        if datum["beta"] != args.beta:
            continue

        bare_mass.append(datum["mAS"])

        m_ps = np.append(
            datum["ps_mass_samples"].samples, datum["ps_mass_samples"].mean
        )

        f_ps = np.append(
            (datum["ps_decay_constant_samples"].samples),
            (datum["ps_decay_constant_samples"].mean),
        )

        mpcac = np.append(datum["mPCAC_samples"].samples, datum["mPCAC_samples"].mean)

        X.append(f_ps**2)
        Y.append(m_ps**2 / mpcac)

    return np.log(np.array(X)), np.log(np.array(Y)), np.array(bare_mass)


def main():
    args = get_args()
    data = read_sample_files(args.data_filenames)

    X, Y, bare_mass = prepare_data(data, args)

    fit_val, chisquare = meson_beta(X, Y)

    y = fit_val[1] * 2 + 1

    fit_A = ufloat(fit_val[0, -1], fit_val[0, 0:-1].std())
    fit_B = ufloat(fit_val[1, -1], fit_val[1, 0:-1].std())
    fit_y = ufloat(y[-1], y.std())

    dump_dict(
        {
            "beta": f"{args.beta}",
            "bare_mass_range": f"[{bare_mass.min()}, {bare_mass.max()}]",
            "chi_sqr_dof": chisquare,
            "A": fit_A,
            "B": fit_B,
            "y": fit_y,
        },
        args.output_file_mean,
    )

    if args.output_file_samples:
        dump_samples(
            {
                "beta": f"{args.beta}",
                f"A_b{args.beta}_samples": fit_val[0, 0:-1],
                f"A_b{args.beta}_value": fit_val[0, -1],
                f"B_b{args.beta}_samples": fit_val[1, 0:-1],
                f"B_b{args.beta}_value": fit_val[1, -1],
            },
            args.output_file_samples,
        )


if __name__ == "__main__":
    main()
