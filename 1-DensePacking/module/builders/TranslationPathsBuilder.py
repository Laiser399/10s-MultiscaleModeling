import collections
from typing import Iterable, List, Tuple

import numpy as np

from module.helpers import TranslationBuilderHelper
from module.models import TranslateOperation


class TranslationPathsBuilder:
    """
    Вместо конечного результата (размноженных точек),
    сохраняет последовательности трансформаций,
    которые приведут к конечному результату
    """

    def __init__(self, translations: Iterable[TranslateOperation]):
        self._translations = list(translations)
        self._confirmed_points = []
        self._confirmed_paths = []

    def build_from(self, start_point: np.ndarray):
        queue = collections.deque()
        queue.append((start_point, tuple()))
        while queue:
            point, translations_path = queue.popleft()
            if self._is_replica(point):
                continue
            self._confirmed_points.append(point)
            self._confirmed_paths.append(translations_path)
            for translation in self._translations:
                new_point = translation.translate(point) % 1
                queue.append((new_point, translations_path + (translation,)))

    def get_paths(self) -> List[Tuple[TranslateOperation]]:
        return list(self._confirmed_paths)

    def reset(self):
        self._confirmed_points.clear()
        self._confirmed_paths.clear()

    def _is_replica(self, point: np.ndarray) -> bool:
        return TranslationBuilderHelper.is_replica(self._confirmed_points, point)
