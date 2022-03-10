from module.models import PositionDefinition


class AtomDefinition:
    def __init__(self, name, radius, position_definition: PositionDefinition):
        self._name = name
        self._radius = radius
        self._position_definition = position_definition

    @property
    def name(self):
        return self._name

    @property
    def radius(self):
        return self._radius

    @property
    def position_definition(self):
        return self._position_definition

    def __repr__(self):
        return f'AtomDefinition({self._name}, {self._radius}, {self._position_definition})'
