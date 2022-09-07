#!/usr/bin/env python3

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

"""
Several functions have been copied verbatim or with minor modifications from scikit-optimize:
https://github.com/scikit-optimize/scikit-optimize/blob/0.9.X/skopt/plots.py

BSD 3-Clause License

Copyright (c) 2016-2020 The scikit-optimize developers.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

import argparse
import os
import pickle
from itertools import count

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import LogLocator
from skopt.plots import plot_convergence, plot_evaluations, expected_minimum, expected_minimum_random_sampling
from skopt.space import Categorical


def _evenly_sample(dim, n_points):
    """
    Source:
    https://github.com/scikit-optimize/scikit-optimize/blob/a2369ddbc332d16d8ff173b12404b03fea472492/skopt/plots.py#L1298
    """
    cats = np.array(getattr(dim, 'categories', []), dtype=object)
    if len(cats):  # Sample categoricals while maintaining order
        xi = np.linspace(0, len(cats) - 1, min(len(cats), n_points),
                         dtype=int)
        xi_transformed = dim.transform(cats[xi])
    else:
        bounds = dim.bounds
        # XXX use linspace(*bounds, n_points) after python2 support ends
        xi = np.linspace(bounds[0], bounds[1], n_points)
        xi_transformed = dim.transform(xi)
    return xi, xi_transformed


def partial_dependence_2D(space, model, i, j, samples,
                          n_points=40):
    """
    Source:
    https://github.com/scikit-optimize/scikit-optimize/blob/a2369ddbc332d16d8ff173b12404b03fea472492/skopt/plots.py#L974
    """
    # The idea is to step through one dimension, evaluating the model with
    # that dimension fixed and averaging either over random values or over
    # the given ones in x_val in all other dimensions.
    # (Or step through 2 dimensions when i and j are given.)
    # Categorical dimensions make this interesting, because they are one-
    # hot-encoded, so there is a one-to-many mapping of input dimensions
    # to transformed (model) dimensions.

    # dim_locs[i] is the (column index of the) start of dim i in
    # sample_points.
    # This is usefull when we are using one hot encoding, i.e using
    # categorical values
    dim_locs = np.cumsum([0] + [d.transformed_size for d in space.dimensions])

    def _calc(x, y):
        """
        Helper-function to calculate the average predicted
        objective value for the given model, when setting
        the index1'th dimension of the search-space to the value x
        and setting the index2'th dimension to the value y,
        and then averaging over all samples.
        """
        rvs_ = np.array(samples)  # copy
        rvs_[:, dim_locs[j]:dim_locs[j + 1]] = x
        rvs_[:, dim_locs[i]:dim_locs[i + 1]] = y
        return np.mean(model.predict(rvs_))

    xi, xi_transformed = _evenly_sample(space.dimensions[j], n_points)
    yi, yi_transformed = _evenly_sample(space.dimensions[i], n_points)
    # Calculate the partial dependence for all combinations of these points.
    zi = [[_calc(x, y) for x in xi_transformed] for y in yi_transformed]

    # Convert list-of-list to a numpy array.
    zi = np.array(zi)

    return xi, yi, zi


def _map_categories(space, points, minimum):
    """
    Source:
    https://github.com/scikit-optimize/scikit-optimize/blob/a2369ddbc332d16d8ff173b12404b03fea472492/skopt/plots.py#L1264
    """
    points = np.asarray(points, dtype=object)  # Allow slicing, preserve cats
    iscat = np.repeat(False, space.n_dims)
    min_ = np.zeros(space.n_dims)
    pts_ = np.zeros(points.shape)
    for i, dim in enumerate(space.dimensions):
        if isinstance(dim, Categorical):
            iscat[i] = True
            catmap = dict(zip(dim.categories, count()))
            pts_[:, i] = [catmap[cat] for cat in points[:, i]]
            min_[i] = catmap[minimum[i]]
        else:
            pts_[:, i] = points[:, i]
            min_[i] = minimum[i]
    return pts_, min_, iscat


def _evaluate_min_params(result, params='result',
                         n_minimum_search=None,
                         random_state=None):
    """
    Source:
    https://github.com/scikit-optimize/scikit-optimize/blob/a2369ddbc332d16d8ff173b12404b03fea472492/skopt/plots.py#L1340
    """
    if isinstance(params, str):
        if params == 'result':
            # Using the best observed result
            x_vals = result.x
        elif params == 'expected_minimum':
            if result.space.is_partly_categorical:
                # space is also categorical
                raise ValueError('expected_minimum does not support any'
                                 'categorical values')
            # Do a gradient based minimum search using scipys own minimizer
            if n_minimum_search:
                # If a value for
                # expected_minimum_samples has been parsed
                x_vals, _ = expected_minimum(
                    result,
                    n_random_starts=n_minimum_search,
                    random_state=random_state)
            else:  # Use standard of 20 random starting points
                x_vals, _ = expected_minimum(result,
                                             n_random_starts=20,
                                             random_state=random_state)
        elif params == 'expected_minimum_random':
            # Do a minimum search by evaluating the function with
            # n_samples sample values
            if n_minimum_search is not None:
                # If a value for
                # n_minimum_samples has been parsed
                x_vals, _ = expected_minimum_random_sampling(
                    result,
                    n_random_starts=n_minimum_search,
                    random_state=random_state)
            else:
                # Use standard of 10^n_parameters. Note this
                # becomes very slow for many parameters
                n_minimum_search = 10 ** len(result.x)
                x_vals, _ = expected_minimum_random_sampling(
                    result,
                    n_random_starts=n_minimum_search,
                    random_state=random_state)
        else:
            raise ValueError('Argument ´eval_min_params´ must be a valid'
                             'string (´result´)')
    elif isinstance(params, list):
        assert len(params) == len(result.x), 'Argument' \
                                             '´eval_min_params´ of type list must have same length as' \
                                             'number of features'
        # Using defined x_values
        x_vals = params
    else:
        raise ValueError('Argument ´eval_min_params´ must'
                         'be a string or a list')
    return x_vals


def plot_objective(result, levels=10, n_points=40, n_samples=250, size=2,
                   zscale='linear', dimensions=None, sample_source='random',
                   minimum='result', minimum_validation=None, n_minimum_search=None,
                   plot_dims=None, show_points=True, cmap='viridis_r', vmin=None, vmax=None):
    """
    Source:
    https://github.com/scikit-optimize/scikit-optimize/blob/a2369ddbc332d16d8ff173b12404b03fea472492/skopt/plots.py#L542
    """
    # Here we define the values for which to plot the red dot (2d plot) and
    # the red dotted line (1d plot).
    # These same values will be used for evaluating the plots when
    # calculating dependence. (Unless partial
    # dependence is to be used instead).
    space = result.space
    # Get the relevant search-space dimensions.
    if plot_dims is None:
        # Get all dimensions.
        plot_dims = []
        for row in range(space.n_dims):
            if space.dimensions[row].is_constant:
                continue
            plot_dims.append((row, space.dimensions[row]))
    else:
        plot_dims = space[plot_dims]
    # Number of search-space dimensions we are using.
    n_dims = len(plot_dims)
    if dimensions is not None:
        assert len(dimensions) == n_dims
    x_vals = _evaluate_min_params(result, minimum, n_minimum_search)
    if sample_source == "random":
        x_eval = None
        samples = space.transform(space.rvs(n_samples=n_samples))
    else:
        x_eval = _evaluate_min_params(result, sample_source,
                                      n_minimum_search)
        samples = space.transform([x_eval])
    x_samples, minimum, _ = _map_categories(space, result.x_iters, x_vals)

    if zscale == 'log':
        locator = LogLocator()
    elif zscale == 'linear':
        locator = None
    else:
        raise ValueError("Valid values for zscale are 'linear' and 'log',"
                         " not '%s'." % zscale)

    fig, ax = plt.subplots(n_dims, n_dims,
                           figsize=(size * n_dims, size * n_dims))

    # fig.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.95,
    #                    hspace=0.1, wspace=0.1)

    kwargs = {}
    if vmin is not None:
        kwargs["vmin"] = vmin

    if vmax is not None:
        kwargs["vmax"] = vmax

    for i in range(n_dims):
        for j in range(n_dims):
            if i > j:
                index1, dim1 = plot_dims[i]
                index2, dim2 = plot_dims[j]
                ax_ = ax[i, j]
                xi, yi, zi = partial_dependence_2D(space, result.models[-1],
                                                   index1, index2,
                                                   samples, n_points)
                cs = ax_.contourf(xi, yi, zi, levels,
                                  locator=locator, cmap=cmap, **kwargs)
                fig.colorbar(cs, label="score", use_gridspec=True)  # , ax=ax_)
                if show_points:
                    ax_.scatter(x_samples[:, index2], x_samples[:, index1],
                                c='k', s=10, lw=0.)
                ax_.scatter(minimum[index2], minimum[index1],
                            c=['r'], s=100, lw=0., marker='.')
                if minimum_validation is not None:
                    ax_.scatter(minimum_validation[0], minimum_validation[1], c=["r"], s=75, lw=0., marker="X")
    return ax


def make_cli():
    class SplitFloatArgs(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            setattr(namespace, self.dest, [float(v) for v in values.split(",")])
    cli = argparse.ArgumentParser()
    cli.add_argument("pickle",
                     type=str,
                     help="Path to the .pickle file generated by the param. optimization script.")

    cli.add_argument("-o", "--output-prefix",
                     type=str,
                     required=True,
                     help="Path to output prefix.")

    cli.add_argument("--vmin", type=float)
    cli.add_argument("--vmax", type=float)

    cli.add_argument("--optimal-params-validation",
                     action=SplitFloatArgs,
                     help="Coordinates of the lowest score on the validation datasets.")

    cli.add_argument("--gradient-levels",
                     default=10,
                     help="Number of color steps to use to plot the gradient of the obj. function.")
    return cli


def import_pickled_data(path):
    with open(path, "rb") as f:
        return pickle.load(f)


def plot_wrapper(data, fx):
    plot = fx(data)
    fig = plot.get_figure()
    return fig


if __name__ == "__main__":
    args = make_cli().parse_args()

    data = import_pickled_data(args.pickle)

    out_prefix = args.output_prefix
    output_dir = os.path.dirname(out_prefix)

    optimal_point_validation = args.optimal_params_validation
    assert optimal_point_validation is None or len(optimal_point_validation) == 2

    if output_dir != "" and output_dir != ".":
        os.makedirs(output_dir, exist_ok=True)

    plot_convergence(data)
    plt.tight_layout()
    plt.savefig(f"{out_prefix}_convergence.png", dpi=1200)
    plt.savefig(f"{out_prefix}_convergence.svg")

    plot_evaluations(data, dimensions=["occupancy", "pnooc->pnooc"])
    plt.tight_layout()
    plt.savefig(f"{out_prefix}_evaluations.png", dpi=1200)
    plt.savefig(f"{out_prefix}_evaluations.svg")

    plot_objective(data,
                   levels=int(args.gradient_levels),
                   size=3,
                   dimensions=["occupancy", "pnooc->pnooc"],
                   minimum_validation=optimal_point_validation,
                   vmin=args.vmin,
                   vmax=args.vmax)
    plt.tight_layout()
    plt.savefig(f"{out_prefix}_objective.png", dpi=1200)
    plt.savefig(f"{out_prefix}_objective.svg")
