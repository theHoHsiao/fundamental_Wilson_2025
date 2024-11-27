#!/usr/bin/env python3
from argparse import ArgumentParser, FileType
from uncertainties import ufloat

from .dump import read_sample_files, dump_dict, dump_samples


def get_args():
    parser = ArgumentParser()

    parser.add_argument(
        "data_filename",
        nargs="+",
        metavar="sample_filename",
        help="Filename of a sample file",
    )
    parser.add_argument(
        "--channel",
        choices=["rhoE1"],
        default=None,
        help="Measuring channel",
    )
    parser.add_argument(
        "--smear",
        default="",
        help="Measuring channel",
    )
    parser.add_argument(
        "--output_file_mean",
        type=FileType("w"),
        default="-",
        help="Where to output the mean and uncertainty",
    )
    parser.add_argument(
        "--output_file_samples",
        type=FileType("w"),
        default=None,
        help="Where to output the bootstrap samples",
    )

    return parser.parse_args()


def main():
    args = get_args()
    if args.smear:
        file_head = "smear_"
    else:
        file_head = ""

    datum = read_sample_files(args.data_filename)[0]

    mv_sample = datum[f"{file_head}v_mass_samples"].samples
    mv_mean = datum[f"{file_head}v_mass_samples"].mean

    m_rhoE1_sample = datum[f"{file_head}{args.channel}_mass_samples"].samples
    m_rhoE1_mean = datum[f"{file_head}{args.channel}_mass_samples"].mean

    R_sample = m_rhoE1_sample / mv_sample
    R_mean = m_rhoE1_mean / mv_mean

    metadata = {
        "ensemble_name": datum["ensemble_name"],
        "beta": datum["beta"],
        "mAS": datum["mAS"],
        "Nt": datum["Nt"],
        "Ns": datum["Ns"],
    }
    dump_dict(
        {
            **metadata,
            f"{file_head}{args.channel}_Rmv": ufloat(R_mean, R_sample.std()),
        },
        args.output_file_mean,
    )
    if args.output_file_samples:
        dump_samples(
            {
                **metadata,
                f"{file_head}{args.channel}_Rmv_samples": R_sample,
                f"{file_head}{args.channel}_Rmv_value": R_mean,
            },
            args.output_file_samples,
        )


if __name__ == "__main__":
    main()
