from abc import ABC, abstractmethod
from typing import Sequence

from pydantic import BaseModel


class BaseBlock(ABC, BaseModel):
    @abstractmethod
    def create_block(self) -> str:
        pass


class ControlBlock(BaseBlock):
    calculation: str
    prefix: str
    pseudo_dir: str
    outdir: str

    def create_block(self) -> str:
        return f'&control\n' \
               f'    calculation = \'{self.calculation}\'\n' \
               f'    prefix = \'{self.prefix}\'\n' \
               f'    pseudo_dir = \'{self.pseudo_dir}\'\n' \
               f'    outdir = \'{self.outdir}\'\n' \
               f'/\n'


class SystemBlock(BaseBlock):
    ibrav: int
    A: float
    nat: int
    ntyp: int
    ecutwfc: int
    ecutrho: int

    def create_block(self) -> str:
        return f'&system\n' \
               f'    ibrav = {self.ibrav}\n' \
               f'    A = {self.A}\n' \
               f'    nat = {self.nat}\n' \
               f'    ntyp = {self.ntyp}\n' \
               f'    ecutwfc = {self.ecutwfc}\n' \
               f'    ecutrho = {self.ecutrho}\n' \
               f'/\n'


class ElectronsBlock(BaseBlock):
    def create_block(self) -> str:
        return f'&electrons\n' \
               f'/\n'


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
            lambda x: f'    {x}\n',
            self.lines
        ))
        return f'{self.block_name}{options_str}\n' \
               f'{lines_str}'


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


# todo заменить значения в соответствии с вариантом задачи
base_configuration = QEConfiguration(
    control=ControlBlock(
        calculation='scf',
        prefix='base',
        pseudo_dir='./SSSP_1.1.2_PBE_precision/',
        outdir='./out/test/',
    ),
    system=SystemBlock(
        ibrav=2,
        A=6.46,
        nat=2,
        ntyp=1,
        ecutwfc=70,
        ecutrho=560
    ),
    electrons=ElectronsBlock(),
    custom_blocks=[
        CustomBlock(
            block_name='ATOMIC_SPECIES',
            lines=[
                'Sn 118.71 Sn_pbe_v1.uspp.F.UPF'
            ]
        ),
        CustomBlock(
            block_name='ATOMIC_POSITIONS',
            options=['crystal'],
            lines=[
                'Sn 0.0 0.0 0.0',
                'Sn 0.25 0.25 0.25'
            ]
        ),
        CustomBlock(
            block_name='K_POINTS',
            options=['automatic'],
            lines=['8 8 8 0 0 0']
        )
    ]
)

with open('./in/test.txt', 'w') as output_file:
    output_file.write(base_configuration.create_block())
