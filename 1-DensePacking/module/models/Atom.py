import numpy as np

from module.models import AtomDefinition


class Atom:
    def __init__(self, atom_definition: AtomDefinition, point: np.ndarray):
        self._atom_definition = atom_definition
        self._point = point

    @property
    def atom_definition(self):
        return self._atom_definition

    @property
    def point(self):
        return self._point

    def __repr__(self):
        return f'Atom({self._atom_definition.name}, {self._atom_definition.position_definition.name}, {self._point})'
