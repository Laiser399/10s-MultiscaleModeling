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
    ecutwfc: float
    ecutrho: float

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
                'Ge 72.63 ge_pbe_v1.4.uspp.F.UPF'
            ]
        ),
        CustomBlock(
            block_name='ATOMIC_POSITIONS',
            options=['crystal'],
            lines=[
                'Ge 0.0 0.0 0.0',
                'Ge 0.25 0.25 0.25'
            ]
        ),
        CustomBlock(
            block_name='K_POINTS',
            options=['automatic'],
            lines=['8 8 8 0 0 0']
        )
    ]
)


def fix_K_BLOCK(k_block: CustomBlock, k: int) -> CustomBlock:
    k_block.lines = [f'{k} {k} {k} 0 0 0']
    return k_block


def fix_k_block_within_custom_blocks(custom_blocks: Sequence[CustomBlock], k: int):
    list(map(
        lambda x: fix_K_BLOCK(x, k),
        filter(
            lambda x: x.block_name == 'K_POINTS',
            custom_blocks
        )
    ))


def get_config_input_file_path(config_name: str):
    return f'./in/{config_name}.txt'


def create_start_bash_script(file_path: str, config_names: Sequence[str]):
    with open(file_path, 'w') as output_file:
        output_file.write('#!/bin/bash\n\n')

        output_file.write('set -e;\n\n')

        for config_name in config_names:
            config_input_file_path = get_config_input_file_path(config_name)

            output_file.write(f'mkdir ./out/{config_name};\n')
            output_file.write(f'pw.x -in {config_input_file_path} > ./out/{config_name}/console_output.txt;\n')
            output_file.write(f'echo "done {config_name}";\n')


def generate_config(configuration: QEConfiguration, config_name, k: int, ecutwfc: float):
    ecutrho_coefficient = 8

    config_input_file_path = get_config_input_file_path(config_name)

    configuration.control.outdir = f'./out/{config_name}'
    configuration.system.ecutwfc = ecutwfc
    configuration.system.ecutrho = ecutrho_coefficient * ecutwfc
    fix_k_block_within_custom_blocks(configuration.custom_blocks, k)

    with open(config_input_file_path, 'w') as output_file:
        output_file.write(configuration.create_block())


def generate_ecutwfc(k: int):
    configuration = base_configuration.copy(deep=True)
    config_names = []
    for ecutwfc in range(10, 100, 10):
        config_name = f'ecutwfc_{k}_{ecutwfc}'

        generate_config(configuration, config_name, k, ecutwfc)

        config_names.append(config_name)

    create_start_bash_script('start_ecutwfc.sh', config_names)


def generate_k(ecutwfc: float):
    configuration = base_configuration.copy(deep=True)
    config_names = []
    for k in range(4, 16, 2):
        config_name = f'k_{k}_{ecutwfc}'

        generate_config(configuration, config_name, k, ecutwfc)

        config_names.append(config_name)

    create_start_bash_script('start_k.sh', config_names)


# dos2unix ./start_ecutwfc.sh; bash ./start_ecutwfc.sh;
generate_ecutwfc(8)
# dos2unix ./start_k.sh; bash ./start_k.sh;
generate_k(60)
