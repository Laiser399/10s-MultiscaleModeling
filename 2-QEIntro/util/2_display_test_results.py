import matplotlib.pyplot as plt

from src.common import get_energy


def get_energy_from_k_test(k: int, ecutwfc: float):
    config_name = f'k_{k}_{ecutwfc}'
    file_path = f'./out/{config_name}/console_output.txt'
    return get_energy(file_path)


def get_energy_from_ecutwfc_test(k: int, ecutwfc: float):
    config_name = f'ecutwfc_{k}_{ecutwfc}'
    file_path = f'./out/{config_name}/console_output.txt'
    return get_energy(file_path)


def show_ecutwfc_test():
    fixed_k = 8

    ecutwfc_list = range(10, 100, 10)
    energies = [
        get_energy_from_ecutwfc_test(fixed_k, ecutwfc)
        for ecutwfc in ecutwfc_list
    ]

    plt.plot(ecutwfc_list, energies, '-o')
    plt.xlabel('ecutwfc')
    plt.ylabel('Total energy (Ry)')
    plt.title(f'Fixed k = {fixed_k}')
    plt.legend(['E(ecutwfc), k = const'])
    plt.show()


def show_k_test():
    fixed_ecutwfc = 60

    k_list = [k for k in range(4, 16, 2)]
    energies = [
        get_energy_from_k_test(k, fixed_ecutwfc)
        for k in k_list
    ]

    plt.plot(k_list, energies, '-o')
    plt.xlabel('k')
    plt.ylabel('Total energy (Ry)')
    plt.legend(['E(k), ecutwfc = const'])
    plt.title(f'Fixed ecutwfc = {fixed_ecutwfc}')
    plt.show()


show_ecutwfc_test()
show_k_test()
