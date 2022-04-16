from abc import ABC, abstractmethod
from typing import Sequence

from pydantic.dataclasses import dataclass


class BaseBlock(ABC):
    @abstractmethod
    def create_block(self) -> str:
        pass


@dataclass
class ControlBlock(BaseBlock):
    calculation: str
    prefix: str
    pseudo_dir: str
    outdir: str

    def create_block(self) -> str:
        return f'&control\n' \
               f'\tcalculation = \'{self.calculation}\'\n' \
               f'\tprefix = \'{self.prefix}\'\n' \
               f'\tpseudo_dir = \'{self.pseudo_dir}\'\n' \
               f'\toutdir = \'{self.outdir}\'\n' \
               f'/\n'


@dataclass
class SystemBlock(BaseBlock):
    ibrav: int
    A: float
    nat: int
    ntyp: int
    ecutwfc: int
    ecutrho: int

    def create_block(self) -> str:
        return f'&system\n' \
               f'\tibrav = {self.ibrav}\n' \
               f'\tA = {self.A}\n' \
               f'\tnat = {self.nat}\n' \
               f'\tntyp = {self.ntyp}\n' \
               f'\tecutwfc = {self.ecutwfc}\n' \
               f'\tecutrho = {self.ecutrho}\n' \
               f'/\n'


@dataclass
class ElectronsBlock(BaseBlock):
    def create_block(self) -> str:
        return f'&electrons\n' \
               f'/\n'


@dataclass
class CustomBlock(BaseBlock):
    block_name: str
    options: Sequence[str] = tuple()
    lines: Sequence[str] = tuple()

    def create_block(self) -> str:
        options_str = ''.join(map(
            lambda x: f' {x}',
            self.options
        ))
        lines_str = ''.join(map(
            lambda x: f'\t{x}\n',
            self.lines
        ))
        return f'{self.block_name}{options_str}\n' \
               f'{lines_str}'


@dataclass
class QEConfiguration(BaseBlock):
    control: ControlBlock
    system: SystemBlock
    electrons: ElectronsBlock
    custom_blocks: Sequence[CustomBlock]

    def create_block(self) -> str:
        custom_blocks_str = ''.join(map(
            lambda x: x.create_block(),
            self.custom_blocks
        ))

        return self.control.create_block() \
               + self.system.create_block() \
               + self.electrons.create_block() \
               + custom_blocks_str


c = QEConfiguration(
    ControlBlock(
        'scf',
        'cubic',
        './SSSP_1.1.2_PBE_precision/',
        './out/',
    ),
    SystemBlock(
        1,
        6.46,
        8,
        1,
        70,
        560
    ),
    ElectronsBlock(),
    [
        CustomBlock(
            'ATOMIC_SPECIES',
            lines=[
                'Sn 118.71 Sn_pbe_v1.uspp.F.UPF'
            ]
        )
    ]
)

print(c.create_block())
