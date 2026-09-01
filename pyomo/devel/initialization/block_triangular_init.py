# ____________________________________________________________________________________
#
# Pyomo: Python Optimization Modeling Objects
# Copyright (c) 2008-2026 National Technology and Engineering Solutions of Sandia, LLC
# Under the terms of Contract DE-NA0003525 with National Technology and Engineering
# Solutions of Sandia, LLC, the U.S. Government retains certain rights in this
# software.  This software is distributed under the 3-clause BSD License.
# ____________________________________________________________________________________

from pyomo.core.base.block import BlockData
from pyomo.contrib.solver.common.base import SolverBase
from pyomo.contrib.solver.common.results import SolutionStatus
from pyomo.contrib.solver.solvers.scip.scip_direct import ScipDirect
from pyomo.contrib.solver.solvers.scip.scip_persistent import ScipPersistent
from pyomo.contrib.solver.solvers.gurobi.gurobi_direct_minlp import GurobiDirectMINLP
import logging

# ===========IMPORTS FROM IDAES BLOCK TRIANGULARIZATION====================
from pyomo.environ import check_optimal_termination
from pyomo.common.config import Bool, ConfigDict, ConfigValue
from pyomo.contrib.incidence_analysis import (
    IncidenceGraphInterface,
    solve_strongly_connected_components,
)

# from idaes.core.initialization.initializer_base import (
#     InitializerBase,
#     InitializationStatus,
# )
# from idaes.core.solvers import get_solver
# from idaes.core.util.exceptions import InitializationError
# from idaes.core.util.model_statistics import number_activated_constraints
# =========================================================================


logger = logging.getLogger(__name__)

# Heavily influenced by BlockTriangularizationInitializer in IDAES-PSE
# Meant to get degrees of freedom to zero automatically, and then initialize with
# strongly_connected_components like in above framework.

# Needs graph analysis
def _check_matching_and_report_dof(nlp):

    igraph = IncidenceGraphInterface(nlp, include_inequality=False)
    max_matching = igraph.maximum_matching()
    degrees_of_freedom = len(max_matching) - len(igraph.variables)
    print(f"The provided model has {degrees_of_freedom} degrees of freedom.")
    return degrees_of_freedom

def _solve_blocks_with_scc(nlp, nlp_solver, solver_kw=None):

    calc_var_opt = {"linesearch": False}

    res_list = solve_strongly_connected_components(
        nlp,
        solver=nlp_solver,
        solve_kwds=solver_kw,
        calc_var_kwds=calc_var_opt
        # use_calc_var=False,
    )
    return res_list



# Needs dof reduction and value setting (most unknowns)
# Needs matching check for maximum/perfect matching
# Needs connection to solve_strongly_connected_components
# Needs try and retry nlp solve outside.


def _initialize_with_block_triangularization(nlp: BlockData, nlp_solver: SolverBase):
 
    dof = _check_matching_and_report_dof(nlp)
    # If dof != 0:
        # find_perfect_matching


    # Make keywords to share into solver for ssc
    scc_opts = {
        "load_solutions": False,
        "raise_exception_on_nonoptimal_result": False,
    }
    # Once degrees of freedom is 0, solve nlp with scc
    res_list = _solve_blocks_with_scc(nlp, 
                                      nlp_solver=nlp_solver,
                                      solver_kw=scc_opts,

                                      )
    feasible_res_list = []
    for res in res_list:
        if res.solution_status in {SolutionStatus.feasible, SolutionStatus.optimal}:
            res.solution_loader.load_vars()
            feasible_res_list.append(res)

    print(f"Strongly connected components finished with {len(feasible_res_list)}/{len(res_list)} feasible components.")

    # Might need to expand more, but for now leaving it here and returning res_list
    return res_list

