import re

total_energy_regex = re.compile(r'!\s*total energy\s*=\s*(.*) Ry\n')


def get_energy(file_path: str) -> float:
    with open(file_path, 'r') as input_file:
        content = input_file.read()

    match = total_energy_regex.search(content)

    return float(match.group(1))
