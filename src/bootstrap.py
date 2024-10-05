#!/usr/bin/env python3

import hashlib
import numpy as np

from flow_analysis.stats.bootstrap import (
    basic_bootstrap,
    sample_bootstrap_1d,
    bootstrap_finalize,
    BOOTSTRAP_SAMPLE_COUNT,
)

from uncertainties import ufloat


def get_rng(name):
    filename = name.strip("/")
    filename_hash = hashlib.md5(filename.encode("utf8")).digest()
    seed = abs(int.from_bytes(filename_hash, "big"))
    return np.random.default_rng(seed)


class BootstrapSampleSet:
    def __init__(self, samples):
        assert len(samples) == BOOTSTRAP_SAMPLE_COUNT
        self.samples = np.asarray(samples)

    def __mul__(self, other):
        if isinstance(other, BootstrapSampleSet):
            return BootstrapSampleSet(self.samples * other.samples)
        else:
            return BootstrapSampleSet(self.samples * other)

    def __truediv__(self, other):
        if isinstance(other, BootstrapSampleSet):
            return BootstrapSampleSet(self.samples / other.samples)
        else:
            return BootstrapSampleSet(self.samples / other)

    def __repr__(self):
        return f"BootstrapSampleSet[mean={self.mean()}, std={self.std()}]"

    def __format__(self, format_spec):
        return f"{{:{format_spec}}}".format(self.to_ufloat())

    def mean(self):
        return self.samples.mean()

    def std(self):
        return self.samples.std()

    def to_ufloat(self):
        return ufloat(self.mean(), self.std())


__all__ = [
    "basic_bootstrap",
    "sample_bootstrap_1d",
    "bootstrap_finalize",
    "BOOTSTRAP_SAMPLE_COUNT",
    "get_rng",
]
