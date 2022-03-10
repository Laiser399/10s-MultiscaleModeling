import collections
from typing import Iterable, List

import numpy as np

from module.helpers.TranslationBuilderHelper import TranslationBuilderHelper
from module.models.TranslateOperation import TranslateOperation


class TranslationPointsBuilder:
    def __init__(self, translations: Iterable[TranslateOperation]):
        self._translations = list(translations)
        self._confirmed_points = []

    def build_from(self, start_point: np.ndarray):
        points_to_check = collections.deque()
        points_to_check.append(start_point)
        while points_to_check:
            point = points_to_check.popleft()
            if self._is_replica(point):
                continue
            self._confirmed_points.append(point)
            for translation in self._translations:
                new_point = translation.translate(point) % 1
                points_to_check.append(new_point)

    def get_result(self) -> List[np.ndarray]:
        return list(self._confirmed_points)

    def reset(self):
        self._confirmed_points.clear()

    def _is_replica(self, point: np.ndarray) -> bool:
        return TranslationBuilderHelper.is_replica(self._confirmed_points, point)
