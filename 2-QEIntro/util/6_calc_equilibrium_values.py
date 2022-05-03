import os.path

import matplotlib.pyplot as plt
import numpy as np

from src.common import get_energy

base_output_dir = './out/energies_sampling/'
base_lattice_constant = 5.7588793963
ry_to_ev_coefficient = 13.6056980659
# Коэффициент перевода из эВ/А^3 в ГПа
hard_to_name_coefficient = 160.2


def get_config_file_path(config_name):
    return f'./in/{config_name}.txt'


def get_console_output_file_path(config_name):
    return os.path.join(base_output_dir, config_name, 'console_output.txt')


def get_config_name(coefficient):
    return f'test_{coefficient}'


def extract_volumes_and_energies():
    lattice_coefficients = [1 + i / 100 for i in range(-4, 5)]
    configuration_names = [get_config_name(coefficient) for coefficient in lattice_coefficients]

    volume_list = []
    energy_list = []
    for name, coefficient in zip(configuration_names, lattice_coefficients):
        lattice_constant = base_lattice_constant * coefficient
        output_file_path = get_console_output_file_path(name)
        energy = ry_to_ev_coefficient * get_energy(output_file_path)

        volume_list.append(lattice_constant ** 3)
        energy_list.append(energy)

    return volume_list, energy_list


volume_list, energy_list = extract_volumes_and_energies()

poly_coefficients = np.polyfit(volume_list, energy_list, 2)


def show_graph():
    plt.plot(volume_list, energy_list, '.', label='Actual energies')
    arg = np.linspace(170, 215, 100)
    values = np.polyval(poly_coefficients, arg)
    plt.plot(arg, values, label='Fitted curve')
    plt.xlabel('Volume [$\\AA^3$]')
    plt.ylabel('Energy [eV]')
    plt.title('Fit result')
    plt.tight_layout()
    plt.legend()
    plt.show()


show_graph()

equilibrium_volume = -poly_coefficients[1] / poly_coefficients[0] / 2 * 4
equilibrium_lattice_constant = equilibrium_volume ** (1 / 3)
bulk_modulus = -poly_coefficients[1] * hard_to_name_coefficient

print(f'{equilibrium_volume = } A^3')
print(f'{equilibrium_lattice_constant = } A')
print(f'{bulk_modulus = } GPa')
