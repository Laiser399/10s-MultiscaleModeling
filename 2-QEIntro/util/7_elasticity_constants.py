import os.path
import posixpath
from typing import Callable, List

import numpy as np
from pydantic import BaseModel

from src.quantum_espresso import QEConfiguration, ControlBlock, SystemBlock, ElectronsBlock, CustomBlock, \
    CellParametersBlock, IonsBlock

base_input_dir = './in/elasticity_constants/'
base_output_dir = './out/elasticity_constants/'


def create_B(alpha):
    return np.array([
        [1 + alpha, 0, 0],
        [0, 1 + alpha, 0],
        [0, 0, 1 + alpha],
    ])


def create_C_11(alpha):
    return np.array([
        [1 + alpha, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
    ])


def create_C_22(alpha):
    return np.array([
        [1, 0, 0],
        [0, 1 + alpha, 0],
        [0, 0, 1],
    ])


def create_C_33(alpha):
    return np.array([
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1 + alpha],
    ])


def create_C_44(alpha):
    denominator = (1 - alpha ** 2) ** (1 / 3)
    return np.array([
        [1 / denominator, 0, 0],
        [0, 1 / denominator, alpha / denominator],
        [0, alpha / denominator, 1 / denominator],
    ])


def create_C_55(alpha):
    denominator = (1 - alpha ** 2) ** (1 / 3)
    return np.array([
        [1 / denominator, 0, alpha / denominator],
        [0, 1 / denominator, 0],
        [alpha / denominator, 0, 1 / denominator],
    ])


def create_C_66(alpha):
    denominator = (1 - alpha ** 2) ** (1 / 3)
    return np.array([
        [1 / denominator, alpha / denominator, 0],
        [alpha / denominator, 1 / denominator, 0],
        [0, 0, 1 / denominator],
    ])


def create_C_12(alpha):
    denominator = (1 - alpha ** 2) ** (1 / 3)
    return np.array([
        [(1 + alpha) / denominator, 0, 0],
        [0, (1 - alpha) / denominator, 0],
        [0, 0, 1 / denominator],
    ])


def create_C_13(alpha):
    denominator = (1 - alpha ** 2) ** (1 / 3)
    return np.array([
        [(1 + alpha) / denominator, 0, 0],
        [0, 1 / denominator, 0],
        [0, 0, (1 - alpha) / denominator],
    ])


def create_C_23(alpha):
    denominator = (1 - alpha ** 2) ** (1 / 3)
    return np.array([
        [1 / denominator, 0, 0],
        [0, (1 + alpha) / denominator, 0],
        [0, 0, (1 - alpha) / denominator],
    ])


class DeformationDescriptor(BaseModel):
    name: str
    deformation_matrix_generator: Callable[[float], np.ndarray]
    relax: bool = False


class DeformationParamsDescriptor(BaseModel):
    name: str
    alpha: float


class QEConfigParams(BaseModel):
    name: str
    cell_parameters_matrix: np.ndarray
    relax: bool = False

    class Config:
        arbitrary_types_allowed = True


deformation_descriptors = [
    DeformationDescriptor(name='B', deformation_matrix_generator=create_B),
    DeformationDescriptor(name='C_11', deformation_matrix_generator=create_C_11),
    DeformationDescriptor(name='C_22', deformation_matrix_generator=create_C_22),
    DeformationDescriptor(name='C_33', deformation_matrix_generator=create_C_33),
    DeformationDescriptor(name='C_44', deformation_matrix_generator=create_C_44, relax=True),
    DeformationDescriptor(name='C_55', deformation_matrix_generator=create_C_55, relax=True),
    DeformationDescriptor(name='C_66', deformation_matrix_generator=create_C_66, relax=True),
    DeformationDescriptor(name='C_12', deformation_matrix_generator=create_C_12),
    DeformationDescriptor(name='C_13', deformation_matrix_generator=create_C_13),
    DeformationDescriptor(name='C_23', deformation_matrix_generator=create_C_23),
]

deformation_params_descriptors = [
    DeformationParamsDescriptor(name='low', alpha=-0.01),
    DeformationParamsDescriptor(name='high', alpha=0.01),
]


def get_input_file_path(config_name):
    return posixpath.join(base_input_dir, f'{config_name}.txt')


def get_output_dir_path(config_name):
    return posixpath.join(base_output_dir, config_name)


def get_console_output_file_path(config_name):
    return posixpath.join(base_output_dir, config_name, 'console_output.txt')


def create_configuration(name, cell_parameters_matrix, relax=False):
    calculation = 'relax' if relax else 'scf'
    ions = IonsBlock() if relax else None

    configuration = QEConfiguration(
        control=ControlBlock(
            calculation=calculation,
            prefix=name,
            pseudo_dir='./SSSP_1.1.2_PBE_precision/',
            outdir=get_output_dir_path(name),
        ),
        system=SystemBlock(
            ibrav=0,
            A=5.758879396279999,
            nat=2,
            ntyp=1,
            ecutwfc=60,
            ecutrho=480,
        ),
        electrons=ElectronsBlock(
            conv_thr=1e-8
        ),
        ions=ions,
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
                lines=[
                    '10 10 10 0 0 0'
                ]
            ),
            CellParametersBlock(
                options=['alat'],
                vectors=cell_parameters_matrix
            )
        ]
    )

    return configuration


def create_qe_config_params() -> List[QEConfigParams]:
    base_cell_parameters = np.array([
        [0, 0.5, 0.5],
        [0.5, 0, 0.5],
        [0.5, 0.5, 0],
    ])

    deformations = [
        QEConfigParams(name='normal', cell_parameters_matrix=base_cell_parameters)
    ]

    for deformation_descr in deformation_descriptors:
        for params_descr in deformation_params_descriptors:
            deformation_matrix = deformation_descr.deformation_matrix_generator(params_descr.alpha)
            cell_parameters_matrix = base_cell_parameters.dot(deformation_matrix)

            name_parts = [deformation_descr.name, params_descr.name]
            if deformation_descr.relax:
                name_parts.append('relax')
            name = '_'.join(name_parts)

            deformations.append(QEConfigParams(
                name=name,
                cell_parameters_matrix=cell_parameters_matrix,
                relax=deformation_descr.relax
            ))

    return deformations


def create_input_dir_if_needed():
    if not os.path.exists(base_input_dir):
        os.mkdir(base_input_dir)


def save_configuration(configuration_name, configuration):
    file_path = get_input_file_path(configuration_name)
    with open(file_path, 'w') as output_file:
        output_file.write(configuration.create_block())


def generate_configurations(configs_params: List[QEConfigParams]):
    for params in configs_params:
        configuration = create_configuration(
            params.name,
            params.cell_parameters_matrix,
            params.relax
        )
        save_configuration(params.name, configuration)


def create_start_bash_script(script_path, configuration_names):
    with open(script_path, 'w') as output_file:
        output_file.write('#!/bin/bash\n\n')

        output_file.write('set -e;\n\n')
        for name in configuration_names:
            input_path = get_input_file_path(name)
            output_path = get_output_dir_path(name)
            console_output_path = get_console_output_file_path(name)

            output_file.write(f'mkdir -p {output_path};\n')
            output_file.write(f'pw.x -in {input_path} > {console_output_path};\n')
            output_file.write(f'echo "done {name}";\n')


create_input_dir_if_needed()

configs_params = create_qe_config_params()
generate_configurations(configs_params)
print(f'Generated {len(configs_params)} configurations.')

script_path = './start_elasticity_constant.sh'
create_start_bash_script(script_path, map(
    lambda x: x.name,
    configs_params
))
print(f'Generated start script at path "{script_path}"')

print('Start command:')
start_command = f'dos2unix {script_path}; bash {script_path};'
print(f'    {start_command}')
