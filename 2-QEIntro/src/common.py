import re

total_energy_regex = re.compile(r'!\s*total energy\s*=\s*(.*) Ry\n')


def get_energy(file_path: str) -> float:
    with open(file_path, 'r') as input_file:
        content = input_file.read()

    matches = list(total_energy_regex.finditer(content))
    return float(matches[-1].group(1))
