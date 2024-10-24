#!/usr/bin/env python3
from .dump import read_sample_files
from .fitting import meson_beta
from argparse import ArgumentParser, FileType
import numpy as np
from .dump import dump_dict, dump_samples
from uncertainties import ufloat
from .mass import C_R


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
    m_ps_sqr = []
    m_ch_sqr = []

    for datum in data:
        if "ps_mass_samples" not in datum:
            continue
        if datum["beta"] != args.beta:
            continue

        m_ps = np.append(
            datum["ps_mass_samples"].samples, datum["ps_mass_samples"].mean
        )

        Z_factor = 1 + 2 * (C_R("ps")) * (8 / datum["beta"]) / (
            16 * 3.141592653589793**2 * datum["plaquette_samples"].samples
        )

        Z_factor_mean = 1 + 2 * (C_R("ps")) * (8 / datum["beta"]) / (
            16 * 3.141592653589793**2 * datum["plaquette_samples"].mean
        )

        f_ps = np.append(
            (datum["ps_matrix_element_samples"].samples * Z_factor),
            (datum["ps_matrix_element_samples"].mean * Z_factor_mean),
        )

        m_ps_sqr.append(m_ps**2)
        m_ch_sqr.append(f_ps**2)

    return np.array(m_ps_sqr), np.array(m_ch_sqr)


def main():
    args = get_args()
    data = read_sample_files(args.data_filenames)

    m_ps, f_ps = prepare_data(data, args)

    fit_val, X2 = meson_beta(m_ps, f_ps)

    fit_A = ufloat(fit_val[0, -1], fit_val[0, 0:-1].std())
    fit_B = ufloat(fit_val[1, -1], fit_val[1, 0:-1].std())

    dump_dict(
        {
            "beta": f"{args.beta}",
            "chi_sqr_dof": X2,
            "A": fit_A,
            "B": fit_B,
        },
        args.output_file_mean,
    )

    if args.output_file_samples:
        dump_samples(
            {
                "beta": f"{args.beta}",
                f"A_{args.beta}_samples": fit_val[0, 0:-1],
                f"A_{args.beta}_value": fit_val[0, -1],
                f"B_{args.beta}_samples": fit_val[1, 0:-1],
                f"B_{args.beta}_value": fit_val[1, -1],
            },
            args.output_file_samples,
        )


if __name__ == "__main__":
    main()
