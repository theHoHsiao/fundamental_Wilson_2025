#!/usr/bin/env python3

import hashlib
import numpy as np

from flow_analysis.stats.bootstrap import (
    basic_bootstrap,
    sample_bootstrap_0d as _sample_bootstrap_0d,
    sample_bootstrap_1d as _sample_bootstrap_1d,
    bootstrap_finalize as _bootstrap_finalize,
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

    def __add__(self, other):
        if isinstance(other, BootstrapSampleSet):
            return BootstrapSampleSet(
                self.mean + other.mean, self.samples + other.samples
            )
        else:
            return BootstrapSampleSet(self.mean + other, self.samples + other)

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        if isinstance(other, BootstrapSampleSet):
            return BootstrapSampleSet(
                self.mean - other.mean, self.samples - other.samples
            )
        else:
            return BootstrapSampleSet(self.mean - other, self.samples - other)

    def __rsub__(self, other):
        return self - other

    def __mul__(self, other):
        if isinstance(other, BootstrapSampleSet):
            return BootstrapSampleSet(
                self.mean * other.mean, self.samples * other.samples
            )
        else:
            return BootstrapSampleSet(self.mean * other, self.samples * other)

    def __rmul__(self, other):
        return self * other

    def __truediv__(self, other):
        if isinstance(other, BootstrapSampleSet):
            return BootstrapSampleSet(
                self.mean / other.mean, self.samples / other.samples
            )
        else:
            return BootstrapSampleSet(self.mean / other, self.samples / other)

    def __rtruediv__(self, other):
        return BootstrapSampleSet(other / self.mean, other / self.samples)

    def __repr__(self):
        return f"BootstrapSampleSet[mean={self.mean}, std={self.std()}]"

    def __format__(self, format_spec):
        return f"{{:{format_spec}}}".format(self.to_ufloat())

    def __getitem__(self, key):
        return BootstrapSampleSet(self.mean[key], self.samples[:, key])

    def arccosh(self):
        return BootstrapSampleSet(np.arccosh(self.mean), np.arccosh(self.samples))

    def std(self):
        return self.samples.std(axis=0)

    def weighted_mean(self):
        """
        Compute the mean along the non-bootstrap dimension of a multidimensional
        BootstrapSampleSet, weighted by the uncertainties.
        """
        return BootstrapSampleSet(self.samples.mean(), self.samples.mean(axis=1))

    def to_ufloat(self):
        if isinstance(self.mean, np.ndarray):
            return [ufloat(mean, std) for mean, std in zip(self.mean, self.std())]
        else:
            return ufloat(self.mean, self.std())


def sample_bootstrap_0d(values, *args, **kwargs):
    values_array = np.asarray(values)
    return BootstrapSampleSet(
        values_array.mean(), _sample_bootstrap_0d(values_array, *args, **kwargs)
    )


def sample_bootstrap_1d(values, *args, **kwargs):
    values_array = np.asarray(values)
    return BootstrapSampleSet(
        values_array.mean(axis=0), _sample_bootstrap_1d(values_array, *args, **kwargs)
    )


def bootstrap_finalize(samples):
    if isinstance(samples, BootstrapSampleSet):
        return samples.to_ufloat()
    else:
        return _bootstrap_finalize(samples)


__all__ = [
    "basic_bootstrap",
    "sample_bootstrap_0d",
    "sample_bootstrap_1d",
    "bootstrap_finalize",
    "BOOTSTRAP_SAMPLE_COUNT",
    "get_rng",
]
