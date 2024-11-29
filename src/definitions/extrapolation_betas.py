#!/usr/bin/env python3

from ..definitions_common import common_definitions_main


def get_definition(data):
    betas = sorted(set([datum["beta"] for datum in data]))
    formatted_betas = (
        r"$\beta = "
        + "$, $".join([f"{beta}" for beta in betas[:-1]])
        + f"$, and ${betas[-1]}$"
    )
    return {"ChiPTBetaValues": formatted_betas}


if __name__ == "__main__":
    common_definitions_main(get_definition, group_key="ensemble_name")
