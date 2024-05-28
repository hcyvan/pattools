import cvxpy as cp
import numpy as np
from scipy.optimize import nnls


def opt_nnls(A, b):
    x, rnorm = nnls(A, b)
    if np.all(x >= 0):
        return x
    return None


def opt_qp(A, b):
    x = cp.Variable(14)
    objective = cp.Minimize(cp.sum_squares(A @ x - b))
    constraints = [x >= 0, cp.sum(x) <= 1]
    problem = cp.Problem(objective, constraints)
    problem.solve()
    status = problem.status
    if status == cp.OPTIMAL:
        x.value[x.value < 0] = 0
        return x.value
    else:
        return None
