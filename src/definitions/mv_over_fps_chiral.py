#!/usr/bin/env python3


from ..definitions_common import common_definitions_main


def get_definition(data):
    result = {}
    for datum in data:
        channel = datum["channel"]
        for latex_suffix, json_prefix in [("Ratio", "R"), ("CL", "L"), ("CW", "W")]:
            result[f"M{channel.upper()}FPSExtrapolation{latex_suffix}"] = (
                "{:.02uSL}".format(
                    datum[f"{json_prefix}_{channel}dfps_samples"].to_ufloat()
                )
            )

    return result


if __name__ == "__main__":
    common_definitions_main(get_definition, group_key="channel")
