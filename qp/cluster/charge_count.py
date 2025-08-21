import os
from glob import glob
from qp.manager.create import get_charge


ATOM_ELECTRON_MAP = {
    "H": 1,
    "B": 5,
    "C": 6,
    "N": 7,
    "O": 8,
    "F": 9,
    "Na": 11,
    "Mg": 12,
    "Al": 13,
    "P": 15,
    "S": 16,
    "Cl": 17,
    "K": 19,
    "Ca": 20,
    "Zn": 30,
    "Se": 34,
    "Br": 35,
    "Rb": 37,
    "Cd": 48,
    "I": 53,
    "Cs": 55,
    "Pb":82,
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
    print(f"charge: {charge}, extraspin: {extraspin} for {os.path.basename(cluster_path)}")
    spinmult = 1 + extraspin
    xyzfile = glob(os.path.join(cluster_path, "*.xyz"))[0]
    num_electron = count_electron(xyzfile)
    if (num_electron - charge) % 2 == spinmult % 2:
        print(f"> ERROR: charge {charge}, spin multiplicity {spinmult}, electron number {num_electron} in {os.path.basename(cluster_path)} are not consistent!")
