from typing import Iterable, Dict, List

from module.builders.TranslationPointsBuilder import TranslationPointsBuilder
from module.helpers.AtomsBuilderHelper import AtomsBuilderHelper
from module.models.Atom import Atom
from module.models.AtomDefinition import AtomDefinition
from module.models.PositionDefinition import PositionDefinition
from module.models.TranslateOperation import TranslateOperation


class AtomsBuilder:
    def __init__(self, translations: Iterable[TranslateOperation]):
        self._translations_builder = TranslationPointsBuilder(translations)

    def build_atoms(self,
                    atoms_definitions: Iterable[AtomDefinition],
                    freedom_coefficients: Dict[PositionDefinition, List[float]] = None) -> List[Atom]:
        if freedom_coefficients is None:
            freedom_coefficients = dict()

        atoms = []
        for atom_definition in atoms_definitions:
            position_definition = atom_definition.position_definition
            point = AtomsBuilderHelper.get_actual_point(position_definition, freedom_coefficients)

            self._translations_builder.build_from(point)
            for point in self._translations_builder.get_result():
                atoms.append(Atom(atom_definition, point))
            self._translations_builder.reset()
        return atoms
