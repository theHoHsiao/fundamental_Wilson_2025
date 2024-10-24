#!/usr/bin/env python3
from .dump import read_sample_files
from . import fitting
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
        "--channel",
        choices=["ps", "v", "t", "av", "at", "s", "rhoE1"],
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
    m_ch_sqr = []

    for datum in data:
        if "w0_samples" not in datum:
            continue
        if "smear_ps_mass_samples" not in datum:
            continue
        if f"smear_{args.channel}_mass_samples" not in datum:
            continue
        if np.isnan((datum[f"smear_{args.channel}_mass_samples"].samples).mean()):
            continue

        w0 = np.append(datum["w0_samples"].samples, datum["w0_samples"].mean)
        m_ps = np.append(
            datum["smear_ps_mass_samples"].samples, datum["smear_ps_mass_samples"].mean
        )
        m_ch = np.append(
            datum[f"smear_{args.channel}_mass_samples"].samples,
            datum[f"smear_{args.channel}_mass_samples"].mean,
        )

        print(m_ps.shape)

        m_ps_sqr.append((w0 * m_ps) ** 2)

        m_ch_sqr.append((w0 * m_ch) ** 2)

        lat_a.append(1 / w0)

    return np.array(m_ps_sqr), np.array(lat_a), np.array(m_ch_sqr)


def main():
    args = get_args()
    data = read_sample_files(args.data_filenames)

    m_ps_sqr, lat_a, m_ch_sqr = prepare_data(data, args)

    fit_val, X2 = fitting.meson_M2(m_ps_sqr, lat_a, m_ch_sqr)

    print(fit_val[0, 0:-1])
    fit_M = ufloat(fit_val[0, -1], fit_val[0, 0:-1].std())
    fit_L = ufloat(fit_val[1, -1], fit_val[1, 0:-1].std())
    fit_W = ufloat(fit_val[2, -1], fit_val[2, 0:-1].std())

    dump_dict(
        {
            "channel": f"m_{args.channel}",
            "chi_sqr_dof": X2,
            f"M_{args.channel}": fit_M,
            f"L_{args.channel}": fit_L,
            f"W_{args.channel}": fit_W,
        },
        args.output_file_mean,
    )

    if args.output_file_samples:
        dump_samples(
            {
                "channel": f"m_{args.channel}",
                f"M_{args.channel}_samples": fit_val[0, 0:-1],
                f"M_{args.channel}_value": fit_val[0, -1],
                f"L_{args.channel}_samples": fit_val[1, 0:-1],
                f"L_{args.channel}_value": fit_val[1, -1],
                f"W_{args.channel}_samples": fit_val[2, 0:-1],
                f"W_{args.channel}_value": fit_val[2, -1],
            },
            args.output_file_samples,
        )


if __name__ == "__main__":
    main()
