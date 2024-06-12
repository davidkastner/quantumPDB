import os
from glob import glob
from qp.job_manager.create import get_charge


ATOM_ELECTRON_MAP = {
    "H": 1,
    "C": 6,
    "N": 7,
    "O": 8,
    "S": 16,
    "Cl": 17,
    "Br": 35,
    "P": 15,
    "Se": 34,
    "I": 53,
    "F": 9
}


def count_electron(xyzfile):
    with open(xyzfile) as f:
        lines = f.readlines()[2:]
    count = 0
    for line in lines:
        if line:
            atom_type = line.split()[0].strip()
            count += ATOM_ELECTRON_MAP.get(atom_type, 0)
    return count


def check_charge(cluster_path):
    charge, extraspin = get_charge(cluster_path)
    spinmult = 1 + extraspin
    xyzfile = glob(os.path.join(cluster_path, "*.xyz"))[0]
    num_electron = count_electron(xyzfile)
    if (num_electron - charge) % 2 == spinmult % 2:
        print(f"charge {charge}, spin multiplicity {spinmult}, electron number {num_electron} in {cluster_path} are not consistent!")
