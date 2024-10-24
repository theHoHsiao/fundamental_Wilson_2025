#!/usr/bin/env python3
from .dump import read_sample_files
from . import fitting
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
        "--channel",
        choices=["ps", "v", "av"],
        default=None,
        help="Measuring channel",
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
    lat_a = []
    R = []

    for datum in data:
        if "w0_samples" not in datum:
            continue
        if "ps_mass_samples" not in datum:
            continue
        if "ps_matrix_element_samples" not in datum:
            continue
        if np.isnan((datum[f"{args.channel}_mass_samples"].samples).mean()):
            continue

        Z_factor = 1 + 2 * (C_R("ps")) * (8 / datum["beta"]) / (
            16 * 3.141592653589793**2 * datum["plaquette_samples"].samples
        )

        Z_factor_mean = 1 + 2 * (C_R("ps")) * (8 / datum["beta"]) / (
            16 * 3.141592653589793**2 * datum["plaquette_samples"].mean
        )

        w0 = np.append(datum["w0_samples"].samples, datum["w0_samples"].mean)

        m_ps = np.append(
            datum["ps_mass_samples"].samples, datum["ps_mass_samples"].mean
        )

        m_ch = np.append(
            datum[f"{args.channel}_mass_samples"].samples,
            datum[f"{args.channel}_mass_samples"].mean,
        )

        f_ps = np.append(
            (datum["ps_matrix_element_samples"].samples * Z_factor),
            (datum["ps_matrix_element_samples"].mean * Z_factor_mean),
        )

        # print(m_ps.shape)

        m_ps_sqr.append((w0 * m_ps) ** 2)

        R.append(m_ch / f_ps)

        lat_a.append(1 / w0)

    return np.array(m_ps_sqr), np.array(lat_a), np.array(R)


def main():
    args = get_args()
    data = read_sample_files(args.data_filenames)

    m_ps_sqr, lat_a, m_ch_sqr = prepare_data(data, args)

    fit_val, X2 = fitting.meson_M2(m_ps_sqr, lat_a, m_ch_sqr)

    # print(fit_val[0, 0:-1])
    fit_M = ufloat(fit_val[0, -1], fit_val[0, 0:-1].std())
    fit_L = ufloat(fit_val[1, -1], fit_val[1, 0:-1].std())
    fit_W = ufloat(fit_val[2, -1], fit_val[2, 0:-1].std())

    dump_dict(
        {
            "channel": f"R_m{args.channel}dfps",
            "chi_sqr_dof": X2,
            f"R_{args.channel}dfps": fit_M,
            f"L_{args.channel}dfps": fit_L,
            f"W_{args.channel}dfps": fit_W,
        },
        args.output_file_mean,
    )

    if args.output_file_samples:
        dump_samples(
            {
                "channel": f"R_m{args.channel}dfps",
                f"R_{args.channel}dfps_samples": fit_val[0, 0:-1],
                f"R_{args.channel}dfps_value": fit_val[0, -1],
                f"L_{args.channel}dfps_samples": fit_val[1, 0:-1],
                f"L_{args.channel}dfps_value": fit_val[1, -1],
                f"W_{args.channel}dfps_samples": fit_val[2, 0:-1],
                f"W_{args.channel}dfps_value": fit_val[2, -1],
            },
            args.output_file_samples,
        )


if __name__ == "__main__":
    main()
