import numpy as np


class CellHelper:
    @staticmethod
    def calc_periodic_distance(a: np.ndarray, b: np.ndarray) -> np.ndarray:
        '''
        Вычисляет расстояние с учетом периодичности границы.
        Например, если расстояние по x внутри решетки равно 0.9, то фактическое расстояние будет равной 0.1
        '''
        return np.linalg.norm((a - b + 0.5) % 1 - 0.5)
