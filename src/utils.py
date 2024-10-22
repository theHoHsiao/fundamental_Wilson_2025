#!/usr/bin/env python3

from more_itertools import pairwise


def get_index_separation(indices):
    separation = indices[1] - indices[0]
    for idx1, idx2 in pairwise(indices):
        if idx2 - idx1 != separation:
            raise NotImplementedError(
                "Configurations have non-uniform separation or are out of order."
            )
    return separation
