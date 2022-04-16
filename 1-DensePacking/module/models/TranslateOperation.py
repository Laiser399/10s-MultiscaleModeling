import numpy as np


class TranslateOperation:
    _mappings = {
        'x': np.array([1, 0, 0]),
        'y': np.array([0, 1, 0]),
        'z': np.array([0, 0, 1]),
        'X': np.array([-1, 0, 0]),
        'Y': np.array([0, -1, 0]),
        'Z': np.array([0, 0, -1]),
    }

    def __init__(self, permutation: str, shift: np.ndarray = None):
        if shift is None:
            shift = np.zeros(3)
        if len(permutation) != 3 or len(shift) != 3:
            raise ValueError
        self._permutation = permutation
        self._permutation_matrix = TranslateOperation._build_permutation_matrix(permutation)
        self._shift = shift

    def translate(self, vector: np.ndarray) -> np.ndarray:
        return self._permutation_matrix.dot(vector) + self._shift

    @staticmethod
    def _build_permutation_matrix(permutation: str):
        permutation_matrix = np.zeros((3, 3))
        for i in range(3):
            permutation_matrix[i] = TranslateOperation._mappings[permutation[i]]
        return permutation_matrix

    def __repr__(self):
        return f'Translate({self._permutation}, {self._shift})'
