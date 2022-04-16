import numpy as np


class UnitConstraint:
    def __init__(self, x, y, z, value):
        self._coefficients = np.array([x, y, z])
        self._value = value

    def check(self, vector: np.ndarray):
        return sum(self._coefficients * vector) > self._value

    def check_including_bounds(self, vector: np.ndarray):
        return sum(self._coefficients * vector) >= self._value
