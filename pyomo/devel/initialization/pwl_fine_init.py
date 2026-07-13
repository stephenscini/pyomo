# ____________________________________________________________________________________
#
# Pyomo: Python Optimization Modeling Objects
# Copyright (c) 2008-2026 National Technology and Engineering Solutions of Sandia, LLC
# Under the terms of Contract DE-NA0003525 with National Technology and Engineering
# Solutions of Sandia, LLC, the U.S. Government retains certain rights in this
# software.  This software is distributed under the 3-clause BSD License.
# ____________________________________________________________________________________

from pyomo.contrib.piecewise.piecewise_linear_expression import (
    PiecewiseLinearExpression,
)
from typing import List, MutableMapping, Sequence
from pyomo.core.base.constraint import ConstraintData
from pyomo.common.collections import ComponentMap, ComponentSet
from pyomo.core.base.block import BlockData
from pyomo.contrib.solver.common.base import SolverBase
from pyomo.contrib.solver.common.results import SolutionStatus
from pyomo.contrib.solver.common.results import Results
from pyomo.contrib.solver.solvers.scip.scip_direct import ScipDirect
from pyomo.contrib.solver.solvers.scip.scip_persistent import ScipPersistent
from pyomo.contrib.solver.solvers.gurobi.gurobi_direct_minlp import GurobiDirectMINLP
from pyomo.contrib.solver.solvers.highs import Highs
import logging
from pyomo.common.modeling import unique_component_name
from pyomo.core.expr.visitor import StreamBasedExpressionVisitor, identify_components
from pyomo.devel.initialization.bounds.bound_variables import (
    bound_all_nonlinear_variables,
)
from pyomo.devel.initialization.utils import (
    fix_vars_with_equal_bounds,
    get_vars,
    shallow_clone,
)
import pyomo.environ as pyo

logger = logging.getLogger(__name__)


def _minimize_infeasibility(m):
    trans = pyo.TransformationFactory('core.add_slack_variables')
    trans.apply_to(m, add_slack_objective=False)

    obj_expr = 0

    found_obj = False
    for obj in m.component_data_objects(pyo.Objective, active=True, descend_into=True):
        if found_obj:
            raise RuntimeError(
                'initialization module currently only supports models '
                'with zero or one active objectives'
            )
        if obj.sense == pyo.minimize:
            obj_expr += 0.1 * obj.expr
        else:
            obj_expr -= 0.1 * obj.expr
        obj.deactivate()
        found_obj = True

    obj_name = unique_component_name(m, 'slack_obj')
    new_obj = 10 * trans.get_summed_slacks_expr(m) + obj_expr
    setattr(m, obj_name, pyo.Objective(expr=new_obj))



def _get_pwl_constraints(
    m: BlockData,
) -> MutableMapping[PiecewiseLinearExpression, List[ConstraintData]]:
    comp_types = set()
    comp_types.add(PiecewiseLinearExpression)
    pwl_expr_to_con_map = ComponentMap()
    con_list = list(
        m.component_data_objects(pyo.Constraint, active=True, descend_into=True)
    )
    obj_list = list(
        m.component_data_objects(pyo.Objective, active=True, descend_into=True)
    )
    for comp in con_list + obj_list:
        pwl_exprs = list(identify_components(comp.expr, comp_types))
        if not pwl_exprs:
            continue
        assert len(pwl_exprs) == 1
        e = pwl_exprs[0]
        if e not in pwl_expr_to_con_map:
            pwl_expr_to_con_map[e] = []
        pwl_expr_to_con_map[e].append(comp)
    return pwl_expr_to_con_map


def _initialize_with_piecewise_linear_approx_fine(
    nlp: BlockData,
    nlp_solver: SolverBase,
    mip_solver: SolverBase,
    default_bound=1.0e8,
    aggressive_substitution=True,
    bounds_tol: float = 1e-6,
) -> Results:
    if isinstance(mip_solver, (ScipDirect, ScipPersistent)):
        opts = {'limits/solutions': 1}
    elif isinstance(mip_solver, (GurobiDirectMINLP, Highs)):
        opts = {'SolutionLimit': 1}
    else:
        raise NotImplementedError(
            'Currently, the initialization module only works with new solver '
            'interfaces, so the global solvers are limited to ScipDirect, '
            'ScipPersistent, and GurobiDirectMINLP.'
        )
    # Check if time limit is provided for global solver
    if mip_solver.config.time_limit is None:
        logger.warning(
            'No time limit set for global optimizer. '
            'For a large model, this may take a long time. '
            'Consider setting a time limit using global_solver.config.time_limit.'
        )

    logger.info('Starting initialization using a piecewise linear approximation')
    pwl = shallow_clone(nlp)
    logger.info('created a shallow clone of the model')

    # first introduce auxiliary variables so that we don't try to
    # approximate any functions of more than two variables
    trans = pyo.TransformationFactory(
        'contrib.piecewise.univariate_nonlinear_decomposition'
    )
    trans.apply_to(pwl, aggressive_substitution=aggressive_substitution)
    logger.info('applied the univariate_nonlinear_decomposition transformation')

    # now we need to try to get bounds on all of the nonlinear variables
    bound_all_nonlinear_variables(pwl, default_bound=default_bound)
    logger.info('bounded nonlinear variables')

    # Now, we need to fix variables with equal (or nearly equal) bounds.
    # Otherwise, the PWL transformation complains
    fix_vars_with_equal_bounds(pwl)
    logger.info('fixed variables with equal bounds')

    # now we modify the model by introducing slacks to make sure the PWL
    # approximation is feasible
    # all of the slacks appear linearly, so we don't need to worry about
    # upper bounds for them
    _minimize_infeasibility(pwl)
    logger.info('reformulated model to minimize infeasibility')

    # build the PWL approximation
    trans = pyo.TransformationFactory('contrib.piecewise.nonlinear_to_pwl')
    trans.apply_to(pwl, num_points=25, additively_decompose=False)
    logger.info('replaced nonlinear expressions with piecewise linear expressions')

    """
    Now we want to 
    1. solve the PWL approximation
    2. Initialize the NLP to the solution
    3. Try solving the NLP
    """

    # PWL transformation (and map the variables)
    orig_vars = list(get_vars(pwl))
    pwl.orig_vars = orig_vars
    trans = pyo.TransformationFactory('contrib.piecewise.disaggregated_logarithmic')
    _pwl = trans.create_using(pwl)
    new_vars = _pwl.orig_vars
    del pwl.orig_vars
    del _pwl.orig_vars
    logger.info('applied the disaggregated logarithmic transformation')

    # solve the MILP
    res = mip_solver.solve(
        _pwl, load_solutions=True, raise_exception_on_nonoptimal_result=False
    )
    logger.info(f'solved MILP: {res.solution_status}, {res.termination_condition}')

    # load the variable values back into orig_vars
    if res.solution_status in {SolutionStatus.feasible, SolutionStatus.optimal}:
        for ov, nv in zip(orig_vars, new_vars):
            ov.set_value(nv.value, skip_validation=True)

    return res