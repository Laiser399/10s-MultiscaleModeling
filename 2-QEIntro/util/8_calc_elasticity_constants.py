import posixpath

from src.common import get_energy

base_output_dir = './out/elasticity_constants/'

ry_to_ev_coefficient = 13.6056980659
# Коэффициент перевода из эВ/А^3 в ГПа
hard_to_name_coefficient = 160.2

alpha = 0.01
lattice_constant = 6.64970915332
system_volume = lattice_constant ** 3 / 4


def get_console_output_file_path(config_name):
    return posixpath.join(base_output_dir, config_name, 'console_output.txt')


transformation_names = [
    'B',
    'C_11',
    'C_22',
    'C_33',
    'C_44',
    'C_55',
    'C_66',
    'C_12',
    'C_13',
    'C_23',
]


def load_energies():
    energies = {
        'normal': ry_to_ev_coefficient * get_energy(get_console_output_file_path('normal'))
    }

    for transformation_name in transformation_names:
        low_config_name = f'{transformation_name}_low'
        high_config_name = f'{transformation_name}_high'
        energies[low_config_name] = ry_to_ev_coefficient * get_energy(get_console_output_file_path(low_config_name))
        energies[high_config_name] = ry_to_ev_coefficient * get_energy(get_console_output_file_path(high_config_name))

    return energies


energies = load_energies()


def calculate_diff(transformation_name):
    low_energy = energies[f'{transformation_name}_low']
    normal_energy = energies['normal']
    high_energy = energies[f'{transformation_name}_high']

    res = (low_energy - 2 * normal_energy + high_energy) / alpha ** 2
    res *= hard_to_name_coefficient
    return res


def calculate_B():
    return calculate_diff('B') / 9 / system_volume


def calculate_C_11():
    return calculate_diff('C_11') / system_volume


def calculate_C_22():
    return calculate_diff('C_22') / system_volume


def calculate_C_33():
    return calculate_diff('C_33') / system_volume


def calculate_C_44():
    return calculate_diff('C_44') / 4 / system_volume


def calculate_C_55():
    return calculate_diff('C_55') / 4 / system_volume


def calculate_C_66():
    return calculate_diff('C_66') / 4 / system_volume


def calculate_C_12():
    tm = calculate_diff('C_11') + calculate_diff('C_22') - calculate_diff('C_12')
    return tm / 2 / system_volume


def calculate_C_13():
    tm = calculate_diff('C_11') + calculate_diff('C_33') - calculate_diff('C_13')
    return tm / 2 / system_volume


def calculate_C_23():
    tm = calculate_diff('C_22') + calculate_diff('C_33') - calculate_diff('C_23')
    return tm / 2 / system_volume


print(f'B = {calculate_B()} GPa')
print(f'C_11 = {calculate_C_11()} GPa')
print(f'C_22 = {calculate_C_22()} GPa')
print(f'C_33 = {calculate_C_33()} GPa')
print(f'C_44 = {calculate_C_44()} GPa')
print(f'C_55 = {calculate_C_55()} GPa')
print(f'C_66 = {calculate_C_66()} GPa')
print(f'C_12 = {calculate_C_12()} GPa')
print(f'C_13 = {calculate_C_13()} GPa')
print(f'C_23 = {calculate_C_23()} GPa')
