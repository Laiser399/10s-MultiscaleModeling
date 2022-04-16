from typing import Dict, List, Iterable, Tuple

import numpy as np

from module.models.Atom import Atom
from module.models.AtomDefinition import AtomDefinition
from module.models.PositionDefinition import PositionDefinition
from module.models.TranslateOperation import TranslateOperation


class AtomsBuilderHelper:
    @staticmethod
    def build_atoms_by_paths(atom_definition: AtomDefinition,
                             paths: Iterable[Tuple[TranslateOperation]],
                             freedom_coefficients: Dict[PositionDefinition, List[float]] = None) -> List[Atom]:
        position_definition: PositionDefinition = atom_definition.position_definition
        point = AtomsBuilderHelper.get_actual_point(position_definition, freedom_coefficients)
        translated_points = AtomsBuilderHelper._apply_translation_paths(point, paths)
        return list(
            map(
                lambda p: Atom(atom_definition, p),
                translated_points
            )
        )

    @staticmethod
    def _apply_translation_paths(point: np.ndarray,
                                 paths: Iterable[Tuple[TranslateOperation]]) -> List[np.ndarray]:
        result = []
        for path in paths:
            translated_point = point.copy()
            for t in path:
                translated_point = t.translate(translated_point)
            result.append(translated_point)
        return result

    @staticmethod
    def get_actual_point(position_definition: PositionDefinition,
                         freedom_coefficients: Dict[PositionDefinition, List[float]]) -> np.ndarray:
        """
        Возвращает точку, с применением сдвига с помощью векторов степеней свободы
        """
        point: np.ndarray = position_definition.base_point.copy()
        if position_definition.freedom_degree > 0:
            coefficients = AtomsBuilderHelper._get_freedom_coefficients(position_definition, freedom_coefficients)

            for i in range(position_definition.freedom_degree):
                point += coefficients[i] * position_definition.freedom_vectors[i]
        return point % 1

    @staticmethod
    def _get_freedom_coefficients(position_definition: PositionDefinition,
                                  freedom_coefficients: Dict[PositionDefinition, List[float]]) -> List[float]:
        if position_definition not in freedom_coefficients:
            raise ValueError(f'For position definition {position_definition} not found freedom coefficient')
        coefficients = freedom_coefficients[position_definition]
        if len(coefficients) != position_definition.freedom_degree:
            raise ValueError(
                f'Count of freedom coefficients is not equal to count of freedom vectors in position definition: '
                f'{position_definition}')
        return coefficients
