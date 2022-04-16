from src.quantum_espresso import QEConfiguration, ControlBlock, SystemBlock, ElectronsBlock, CustomBlock

configuration = QEConfiguration(
    control=ControlBlock(
        calculation='scf',
        prefix='single',
        pseudo_dir='./SSSP_1.1.2_PBE_precision/',
        outdir='./out/single/',
    ),
    system=SystemBlock(
        ibrav=1,
        A=20,
        nat=1,
        ntyp=1,
        ecutwfc=60,
        ecutrho=480,
        occupations='fixed',
        nspin=2,
        tot_magnetization=2
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
                'Ge 0.0 0.0 0.0'
            ]
        ),
        CustomBlock(
            block_name='K_POINTS',
            options=['gamma']
        )
    ]
)


def create_configuration():
    with open('./in/single.txt', 'w') as output_file:
        output_file.write(configuration.create_block())


# set -e; mkdir ./out/single; pw.x -in ./in/single.txt | tee ./out/single/console_output.txt;
create_configuration()
