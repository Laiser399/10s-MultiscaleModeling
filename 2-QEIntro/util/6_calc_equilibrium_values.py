import os.path

import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt

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

        volume_list.append(lattice_constant ** 3 / 4)
        energy_list.append(energy)

    return volume_list, energy_list


# Birch–Murnaghan equation of state
def state_equation(V, V_0, B_0, B_0_derivative, E_0):
    temp0 = (V_0 / V) ** (2 / 3)
    temp1 = (temp0 - 1) ** 3 * B_0_derivative
    temp2 = (temp0 - 1) ** 2 * (6 - 4 * temp0)
    return E_0 + 9 * V_0 * B_0 * (temp1 + temp2) / 16


def fit_state_equation(volume_list, energy_list):
    volume_lb, volume_rb = min(volume_list), max(volume_list)

    min_energy = min(energy_list)
    energy_lb, energy_rb = min_energy * 0.9, min_energy * 1.1
    if energy_lb > energy_rb:
        energy_lb, energy_rb = energy_rb, energy_lb

    optimized_values, _ = opt.curve_fit(
        state_equation,
        volume_list,
        energy_list,
        bounds=(
            [volume_lb, -np.inf, -np.inf, energy_lb],
            [volume_rb, np.inf, np.inf, energy_rb]
        ),
        maxfev=4000
    )

    return optimized_values


volume_list, energy_list = extract_volumes_and_energies()
V_0, B_0, B_0_derivative, E_0 = fit_state_equation(volume_list, energy_list)

print(f'{V_0 = } A^3')
print(f'B_0 = {B_0 * hard_to_name_coefficient} GPa')
print(f'B_0\' = {B_0_derivative}')
print(f'{E_0 = } eV')


def show_graph():
    plt.plot(volume_list, energy_list, '.', label='Actual energies')

    arg = np.linspace(42, 54, 100)
    values = state_equation(arg, V_0, B_0, B_0_derivative, E_0)
    plt.plot(arg, values, label='Fitted curve')

    plt.xlabel('Volume [$\\AA^3$]')
    plt.ylabel('Energy [eV]')
    plt.title('Fit result')
    plt.tight_layout()
    plt.legend()
    plt.show()


show_graph()
