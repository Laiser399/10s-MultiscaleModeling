from src.quantum_espresso import QEConfiguration, ControlBlock, SystemBlock, ElectronsBlock, IonsBlock, CellBlock, \
    CustomBlock

configuration = QEConfiguration(
    control=ControlBlock(
        calculation='vc-relax',
        prefix='base',
        pseudo_dir='./SSSP_1.1.2_PBE_precision/',
        outdir='./out/relax/',
    ),
    system=SystemBlock(
        ibrav=2,
        A=5.66,
        nat=2,
        ntyp=1,
        ecutwfc=60,
        ecutrho=480
    ),
    electrons=ElectronsBlock(
        conv_thr=1e-8
    ),
    ions=IonsBlock(),
    cell=CellBlock(),
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
            lines=['10 10 10 0 0 0']
        )
    ]
)


def create_configuration():
    with open('./in/relax.txt', 'w') as output_file:
        output_file.write(configuration.create_block())


# set -e; mkdir ./out/relax; pw.x -in ./in/relax.txt | tee ./out/relax/console_output.txt;
create_configuration()
