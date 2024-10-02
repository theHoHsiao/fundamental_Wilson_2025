#!/usr/bin/env python3

import hashlib
import numpy as np

from flow_analysis.stats.bootstrap import (
    basic_bootstrap,
    sample_bootstrap_1d,
    BOOTSTRAP_SAMPLE_COUNT,
)


def get_rng(name):
    filename = name.strip("/")
    filename_hash = hashlib.md5(filename.encode("utf8")).digest()
    seed = abs(int.from_bytes(filename_hash, "big"))
    return np.random.default_rng(seed)


__all__ = [
    "basic_bootstrap",
    "sample_bootstrap_1d",
    "BOOTSTRAP_SAMPLE_COUNT",
    "get_rng",
]
