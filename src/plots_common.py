#!/usr/bin/env python3

import matplotlib.pyplot as plt
from uncertainties import UFloat


def save_or_show(fig, filename=None):
    if filename == "/dev/null":
        plt.close(fig)
    elif filename is not None:
        fig.savefig(filename)
        plt.close(fig)
    else:
        plt.show()


def errorbar_ufloat(ax, x, y, *args, **kwargs):
    if (
        isinstance(x[0], UFloat)
        or hasattr(x, "values")
        and isinstance(x.values[0], UFloat)
    ):
        x_values = [xi.nominal_value for xi in x]
        x_errors = [xi.std_dev for xi in x]
    else:
        x_values = x
        x_errors = None

    if (
        isinstance(y[0], UFloat)
        or hasattr(y, "values")
        and isinstance(y.values[0], UFloat)
    ):
        y_values = [yi.nominal_value for yi in y]
        y_errors = [yi.std_dev for yi in y]
    else:
        y_values = y
        y_errors = None

    ax.errorbar(
        x_values,
        y_values,
        xerr=x_errors,
        yerr=y_errors,
        ls="none",
        *args,
        **kwargs,
    )
