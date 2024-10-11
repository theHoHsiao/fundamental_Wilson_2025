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
    def __init__(self, mean, samples):
        self.mean = mean

        if samples is None or len(samples) == 0:
            self.samples = np.ones(BOOTSTRAP_SAMPLE_COUNT) * np.nan
        elif len(samples) != BOOTSTRAP_SAMPLE_COUNT:
            raise ValueError("Bootstrap sample count mismatch")
        else:
            self.samples = np.asarray(samples)

    def __mul__(self, other):
        if isinstance(other, BootstrapSampleSet):
            return BootstrapSampleSet(
                self.mean * other.mean, self.samples * other.samples
            )
        else:
            return BootstrapSampleSet(self.mean * other, self.samples * other)

    def __truediv__(self, other):
        if isinstance(other, BootstrapSampleSet):
            return BootstrapSampleSet(
                self.mean / other.mean, self.samples / other.samples
            )
        else:
            return BootstrapSampleSet(self.mean / other, self.samples / other)

    def __repr__(self):
        return f"BootstrapSampleSet[mean={self.mean}, std={self.std()}]"

    def __format__(self, format_spec):
        return f"{{:{format_spec}}}".format(self.to_ufloat())

    def std(self):
        return self.samples.std()

    def to_ufloat(self):
        return ufloat(self.mean, self.std())


__all__ = [
    "basic_bootstrap",
    "sample_bootstrap_1d",
    "bootstrap_finalize",
    "BOOTSTRAP_SAMPLE_COUNT",
    "get_rng",
]
