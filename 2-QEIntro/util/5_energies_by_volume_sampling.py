import os.path
import posixpath

from src.quantum_espresso import QEConfiguration, ControlBlock, SystemBlock, ElectronsBlock, CustomBlock

base_input_dir = './in/energies_sampling/'
base_output_dir = './out/energies_sampling/'


def get_input_file_path(config_name):
    return posixpath.join(base_input_dir, f'{config_name}.txt')


def get_output_dir_path(config_name):
    return posixpath.join(base_output_dir, config_name)


def get_console_output_file_path(config_name):
    return posixpath.join(base_output_dir, config_name, 'console_output.txt')


def create_configuration(name, lattice_constant):
    return QEConfiguration(
        control=ControlBlock(
            calculation='scf',
            prefix=name,
            pseudo_dir='./SSSP_1.1.2_PBE_precision/',
            outdir=get_output_dir_path(name),
        ),
        system=SystemBlock(
            ibrav=2,
            A=lattice_constant,
            nat=2,
            ntyp=1,
            ecutwfc=60,
            ecutrho=480,
        ),
        electrons=ElectronsBlock(
            conv_thr=1e-8
        ),
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
                lines=[
                    '10 10 10 0 0 0'
                ]
            )
        ]
    )


def save_configuration(configuration_name, configuration):
    file_path = get_input_file_path(configuration_name)
    with open(file_path, 'w') as output_file:
        output_file.write(configuration.create_block())


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


base_lattice_constant = 5.7588793963


def get_config_name(coefficient):
    return f'test_{coefficient}'


def generate_configurations(coefficients):
    configuration_names = [get_config_name(coefficient) for coefficient in coefficients]

    for name, coefficient in zip(configuration_names, coefficients):
        lattice_constant = base_lattice_constant * coefficient
        configuration = create_configuration(name, lattice_constant)
        save_configuration(name, configuration)

    return configuration_names


def create_input_dir_if_needed():
    if not os.path.exists(base_input_dir):
        os.mkdir(base_input_dir)


create_input_dir_if_needed()

coefficients = [1 + i / 100 for i in range(-4, 5)]
configuration_names = generate_configurations(coefficients)
print(f'Generated configurations: {", ".join(configuration_names)}')

script_path = './start_energies_sampling.sh'
create_start_bash_script(script_path, configuration_names)
print(f'Generated start script at path "{script_path}"')

print('Start command:')
start_command = f'dos2unix {script_path}; bash {script_path};'
print(start_command)
