#!/usr/bin/env python3


import numpy as np

from ..definitions_common import format_definitions
from ..tables_common import common_table_main, get_header, get_footer


def get_similar_value(data, similarity_threshold_sigma=2):
    means = data.map(lambda x: x.nominal_value)
    std_devs = data.map(lambda x: x.std_dev)
    weighted_mean = (means / std_devs**2).sum() / (1 / std_devs**2).sum()
    normalised_differences = abs(means - weighted_mean) / std_devs.mean()
    if (normalised_differences > similarity_threshold_sigma).any():
        raise ValueError("These data aren't actually that close.")
    return weighted_mean


def format_table(df):
    header = get_header(
        [
            r"$\beta$",
            r"fit range ($am_0$)",
            r"$\chi^2/N_{\rm d.o.f.}$",
            r"$\tilde{C}$ ",
            r"$Y$",
            r"$y$",
        ],
        "1em",
    )
    footer = get_footer()
    content = []

    for row in df.itertuples():
        if np.isinf(row.chisquare):
            formatted_chisquare = r"$\infty$"
        else:
            formatted_chisquare = f"{row.chisquare:.02f}"
        content.append(
            (
                r"{} & $[{}, {}]$ & {} & ${:.02uSL}$ & ${:.02uSL}$ & ${:.02uSL}$ \\"
                "\n"
            ).format(
                row.beta,
                row.mAS_min,
                row.mAS_max,
                formatted_chisquare,
                row.A,
                row.B,
                row.y,
            )
        )

    mean_y = get_similar_value(df.y)
    definitions = format_definitions({"DEFTCommonYValue": f"{mean_y:.1f}"})
    return header + "".join(content) + footer, definitions


if __name__ == "__main__":
    common_table_main(format_table, index_name="beta", definitions=True)
