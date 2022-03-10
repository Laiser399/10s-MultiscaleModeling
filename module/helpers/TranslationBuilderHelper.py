from typing import List

import numpy as np

from module.helpers import CellHelper


class TranslationBuilderHelper:
    @staticmethod
    def is_replica(points: List[np.ndarray], point_to_check: np.ndarray):
        for point in points:
            distance = CellHelper.calc_periodic_distance(point, point_to_check)
            if distance < 1e-8:
                return True
        return False
