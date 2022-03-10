from typing import Iterable, List

import numpy as np


class PositionDefinition:
    def __init__(self, name, base_point: np.ndarray, freedom_vectors: Iterable[np.ndarray] = None):
        self._name = name
        self._base_point = base_point
        self._freedom_vectors = list(freedom_vectors) if freedom_vectors is not None else []

    @property
    def name(self):
        return self._name

    @property
    def base_point(self) -> np.ndarray:
        return self._base_point

    @property
    def freedom_vectors(self) -> List[np.ndarray]:
        return self._freedom_vectors

    @property
    def is_movable(self):
        return self.freedom_degree > 0

    @property
    def freedom_degree(self) -> int:
        return len(self._freedom_vectors)

    def __repr__(self):
        freedom_part = f', {self._freedom_vectors}' if len(self._freedom_vectors) > 0 else ''
        return f'PositionDefinition({self._name}, {self.base_point}{freedom_part})'
