from typing import Sequence

from src.quantum_espresso import QEConfiguration, ControlBlock, SystemBlock, ElectronsBlock, CustomBlock, BaseBlock

base_configuration = QEConfiguration(
    control=ControlBlock(
        calculation='scf',
        prefix='base',
        pseudo_dir='./SSSP_1.1.2_PBE_precision/',
        outdir='./out/test/',
    ),
    system=SystemBlock(
        ibrav=2,
        A=5.66,
        nat=2,
        ntyp=1,
        ecutwfc=70,
        ecutrho=560
    ),
    electrons=ElectronsBlock(),
    additional_blocks=[
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


def fix_k_block_within_custom_blocks(additional_blocks: Sequence[BaseBlock], k: int):
    for block in additional_blocks:
        if isinstance(block, CustomBlock) and block.block_name == 'K_POINTS':
            fix_K_BLOCK(block, k)


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
    fix_k_block_within_custom_blocks(configuration.additional_blocks, k)

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
