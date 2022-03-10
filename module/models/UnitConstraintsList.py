from typing import Iterable

import numpy as np

from module.models import UnitConstraint


class UnitConstraintsList:
    def __init__(self, constraints: Iterable[UnitConstraint]):
        self._constraints = list(constraints)

    def check(self, vector: np.ndarray):
        for constraint in self._constraints:
            if not constraint.check(vector):
                return False
        return True

    def check_including_bounds(self, vector: np.ndarray):
        for constraint in self._constraints:
            if not constraint.check_including_bounds(vector):
                return False
        return True
