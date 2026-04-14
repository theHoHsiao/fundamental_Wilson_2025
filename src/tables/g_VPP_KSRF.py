from argparse import ArgumentParser, FileType
from uncertainties import ufloat
from ..dump import read_sample_files
from numpy import sqrt
from ..provenance import text_metadata, get_basic_metadata


def get_args():
    parser = ArgumentParser()

    parser.add_argument(
        "data_filenames",
        nargs="+",
        metavar="sample_filename",
        help="Filenames of sample files containing data to plot",
    )
    parser.add_argument(
        "--output_file",
        type=FileType("w"),
        default="-",
        help="Where to place the table",
    )
    return parser.parse_args()


def g_VPP_KSRF(fit_pars):

    for parameter in fit_pars:
        if parameter["channel"] == "v":
            mv = parameter["M_a2_samples"]
        if parameter["channel"] == "ps":
            fps = parameter["F_a2_samples"]
    
    g_mean = sqrt(mv.mean / ( 2 * fps.mean))
    g_samples = sqrt(mv.samples / ( 2 * fps.samples))
    std = g_samples.std()


    return  ufloat(g_mean, std)


def format_table_main(results):
    content = []
    content.append(
        (
            "{:.02uSL}".format(results)
        )
    )
    return "".join(content)




def main():
    args = get_args()
    fit_pars = read_sample_files(args.data_filenames, group_key="channel")
    result = g_VPP_KSRF(fit_pars)
    
    print(text_metadata(get_basic_metadata(), comment_char="%"), file=args.output_file)
    print(format_table_main(result), file=args.output_file)
    


if __name__ == "__main__":
    main()