from abc import ABC, abstractmethod
from typing import Sequence, Optional

from pydantic import BaseModel


class BaseBlock(ABC, BaseModel):
    @abstractmethod
    def create_block(self) -> str:
        pass

    @classmethod
    def _create_field(cls, name: str, value: Optional, quoted: bool = False):
        if value is None:
            return ''
        if quoted:
            return f'    {name} = \'{value}\'\n'
        return f'    {name} = {value}\n'


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
    ecutwfc: float
    ecutrho: float
    occupations: Optional[str]
    nspin: Optional[int]
    tot_magnetization: Optional[int]

    def create_block(self) -> str:
        return f'&system\n' \
               f'    ibrav = {self.ibrav}\n' \
               f'    A = {self.A}\n' \
               f'    nat = {self.nat}\n' \
               f'    ntyp = {self.ntyp}\n' \
               f'    ecutwfc = {self.ecutwfc}\n' \
               f'    ecutrho = {self.ecutrho}\n' \
               f'{self._create_field("occupations", self.occupations, True)}' \
               f'{self._create_field("nspin", self.nspin)}' \
               f'{self._create_field("tot_magnetization", self.tot_magnetization)}' \
               f'/\n'


class ElectronsBlock(BaseBlock):
    conv_thr: Optional[float]

    def create_block(self) -> str:
        conv_thr_str = f'    conv_thr = {self.conv_thr}\n' if self.conv_thr is not None else ''

        return f'&electrons\n' \
               f'{conv_thr_str}' \
               f'/\n'


class IonsBlock(BaseBlock):
    def create_block(self) -> str:
        return f'&ions\n' \
               f'/\n'


class CellBlock(BaseBlock):
    def create_block(self) -> str:
        return f'&cell\n' \
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
    ions: Optional[IonsBlock]
    cell: Optional[CellBlock]
    additional_blocks: Sequence[BaseBlock]

    def create_block(self) -> str:
        ions_str = self.ions.create_block() if self.ions is not None else ''
        cell_str = self.cell.create_block() if self.cell is not None else ''
        additional_blocks_str = ''.join(map(
            lambda x: x.create_block(),
            self.additional_blocks
        ))

        return self.control.create_block() \
               + self.system.create_block() \
               + self.electrons.create_block() \
               + ions_str \
               + cell_str \
               + additional_blocks_str
