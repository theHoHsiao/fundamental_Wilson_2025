#!/usr/bin/env python3


from ..definitions_common import common_definitions_main


def get_definition(data):
    mhat_v = None
    fhat_ps = None
    for datum in data:
        if datum["channel"] == "v":
            mhat_v = datum["M_samples"]
        elif datum["channel"] == "ps":
            fhat_ps = datum["F_samples"]

    ratio = (mhat_v / fhat_ps) ** 0.5
    return {"MHatVFHatPSRatio": f"{ratio.to_ufloat():.02uSL}"}


if __name__ == "__main__":
    common_definitions_main(get_definition, group_key="channel")
