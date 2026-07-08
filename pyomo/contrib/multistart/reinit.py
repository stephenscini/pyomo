# ____________________________________________________________________________________
#
# Pyomo: Python Optimization Modeling Objects
# Copyright (c) 2008-2026 National Technology and Engineering Solutions of Sandia, LLC
# Under the terms of Contract DE-NA0003525 with National Technology and Engineering
# Solutions of Sandia, LLC, the U.S. Government retains certain rights in this
# software.  This software is distributed under the 3-clause BSD License.
# ____________________________________________________________________________________

"""Helper functions for variable reinitialization."""

import logging
import random
from pyomo.common.dependencies.scipy import stats
from pyomo.core.expr.visitor import (
    identify_variables,
)

from pyomo.core import Var

logger = logging.getLogger('pyomo.contrib.multistart')


def rand(val, lb, ub):
    return random.uniform(lb, ub)  # uniform distribution between lb and ub

def latin_hypercube(val, lb, ub, sampler):
    sample = sampler.random(n=1)
    sample = stats.qmc.scale(sample, lb, ub)
    return sample

def _generate_lhs_sample(vlist, config):
    n_vars = len(vlist)
    bnds_list = []
    for v in vlist:
        # the bounds should not be None because we
        # set the bounds to default_bound in
        # bound_all_nonlinear_variables
        lb = v.lb
        ub = v.ub
        bnds_list.append((lb, ub))
    sampler = stats.qmc.LatinHypercube(d=n_vars, seed=config.seed)
    sample = sampler.random(n=config.seed)
    l_bounds = [i[0] for i in bnds_list]
    u_bounds = [i[1] for i in bnds_list]
    sample = stats.qmc.scale(sample, l_bounds, u_bounds)

def midpoint_guess_and_bound(val, lb, ub):
    """Midpoint between current value and farthest bound."""
    far_bound = ub if ((ub - val) >= (val - lb)) else lb  # farther bound
    return (far_bound + val) / 2


def rand_guess_and_bound(val, lb, ub):
    """Random choice between current value and farthest bound."""
    far_bound = ub if ((ub - val) >= (val - lb)) else lb  # farther bound
    return random.uniform(val, far_bound)


def rand_distributed(val, lb, ub, divisions=9):
    """Random choice among evenly distributed set of values between bounds."""
    set_distributed_vals = linspace(lb, ub, divisions)
    return random.choice(set_distributed_vals)


def simple_midpoint(val, lb, ub):
    return (lb + ub) * 0.5


def linspace(lower, upper, n):
    """Linearly spaced range."""
    return [lower + x * (upper - lower) / (n - 1) for x in range(n)]


strategies = {
    "rand": rand,
    "midpoint_guess_and_bound": midpoint_guess_and_bound,
    "rand_guess_and_bound": rand_guess_and_bound,
    "rand_distributed": rand_distributed,
    "midpoint": simple_midpoint,
    "latin_hypercube": latin_hypercube
}


def reinitialize_variables(model, config):
    """Reinitializes all variable values in the model.

    Excludes fixed, noncontinuous, and unbounded variables.

    """
    if config.strategy is "latin_hypercube":
        vlist = list(identify_variables(model, include_fixed=False))
        _
        for v in vlist:

        
    else:
    for var in model.component_data_objects(ctype=Var, descend_into=True):
        if var.is_fixed() or not var.is_continuous():
            continue
        if var.lb is None or var.ub is None:
            if not config.suppress_unbounded_warning:
                logger.warning(
                    'Skipping reinitialization of unbounded variable '
                    '%s with bounds (%s, %s). '
                    'To suppress this message, set the '
                    'suppress_unbounded_warning flag.' % (var.name, var.lb, var.ub)
                )
            continue
        val = var.value if var.value is not None else (var.lb + var.ub) / 2
        # apply reinitialization strategy to variable
        var.set_value(
            strategies[config.strategy](val, var.lb, var.ub), skip_validation=True
        )
