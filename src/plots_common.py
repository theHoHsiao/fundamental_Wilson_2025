#!/usr/bin/env python3

import matplotlib.pyplot as plt


def save_or_show(fig, filename=None):
    if filename == "/dev/null":
        plt.close(fig)
    elif filename is not None:
        fig.savefig(filename)
        plt.close(fig)
    else:
        plt.show()
