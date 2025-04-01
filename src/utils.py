#!/usr/bin/env python3

from more_itertools import pairwise


def get_index_separation(indices):
    indices.sort()
    separation = indices[1] - indices[0]
    for idx1, idx2 in pairwise(indices):
        if idx2 - idx1 != separation:
            raise NotImplementedError(
                f"Configurations have non-uniform separation or are out of order. {idx1}, {idx2}, {separation}\n {indices}"    
            )
    return separation
