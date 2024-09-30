"""Calculate QM properties of the ligands

**Usage**

>>> protonate.compute_charge("path/to/OUT.sdf")
{"AAG_A331": 0, "SO4_A901": -2, "SO4_A325": -2, 
"SO4_A903": -2, "AKG_A330": -2,  "GOL_A328": 0, 
"GOL_A329": 0, "GOL_A900": 0, "GOL_A904": 0}

"""

def compute_charge(path_ligand, path_pdb):
    """
    Computes the total charge of each ligand

    Parameters
    ----------
    path: str
        Path to ligand SDF file

    Returns
    -------
    charge: dict
        Keyed by ligand ID
    """
    with open(path_ligand, "r") as f:
        sdf = f.read()
    with open(path_pdb, "r") as f:
        pdb_lines = f.readlines()
    ligands = [
        [t for t in s.splitlines() if t != ""]
        for s in sdf.split("$$$$")
        if s != "\n" and s != ""
    ]

    charge = {}
    for l in ligands:
        title = l[0].split("_")
        name_id = [
            (res_name, chain_id, res_id) 
            for res_name, chain_id, res_id 
            in zip(title[::3], title[1::3], title[2::3])
        ]
        name = " ".join([f"{res_name}_{chain_id}{res_id}" for res_name, chain_id, res_id in name_id])

        c = 0
        n_atom = 0
        for line in l:
            if "V2000" in line:
                n_atom = int(line[:3]) # https://discover.3ds.com/sites/default/files/2020-08/biovia_ctfileformats_2020.pdf
            if line.startswith("M  RGP"):
                # R# atom is found in the addition between CSO and TAN
                # It's not a real atom and recorded in the RGP entry
                n_atom -= int(line.split()[2])
            if line.startswith("M  CHG"):
                c += sum([int(x) for x in line.split()[4::2]])
                break
        charge[name] = c
        cnt = 0
        for res_name, chain_id, res_id in name_id:
            if res_name == "NO":
                cnt = n_atom
                break
            # to detect removed atoms
            # NO is corrected from NH2OH, thus the deleted 3 hydrogen shouldn't affect the charge
            for line in pdb_lines:
                if line[17:20].strip() == res_name and line[21] == chain_id and line[22:26].strip() == res_id and line[0:3] != "TER":
                    cnt += 1
        charge[name] -= (n_atom - cnt)

    # check Na+, K+, Mg2+, Ca2+
    for line in pdb_lines:
        if len(line) > 26 and line[0:3] != "TER":
            res_name = line[17:20].strip()
            chain_id = line[21]
            res_id = line[22:26].strip()
            if res_name in {"NA", "K"}:
                print(f"Find {res_name} in protein, charge + 1")
                charge[f"{res_name}_{chain_id}{res_id}"] = 1
            elif res_name in {"MG", "CA"}:
                print(f"Find {res_name} in protein, charge + 2")
                charge[f"{res_name}_{chain_id}{res_id}"] = 2
    return charge


def compute_spin(path_ligand):
    with open(path_ligand, "r") as f:
        sdf = f.read()

    ligands = [
        [t for t in s.splitlines() if t != ""]
        for s in sdf.split("$$$$")
        if s != "\n" and s != ""
    ]
    spin = {}
    for l in ligands:
        title = l[0].split("_")
        name_id = [
            (res_name, chain_id, res_id) 
            for res_name, chain_id, res_id 
            in zip(title[::3], title[1::3], title[2::3])
        ]
        name = " ".join([f"{res_name}_{chain_id}{res_id}" for res_name, chain_id, res_id in name_id])
        for res_name, chain_id, res_id in name_id:
            if res_name == "NO":
                spin[name] = 1
            elif res_name == "OXY":
                spin[name] = 2
    return spin
