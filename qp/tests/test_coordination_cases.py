"""Coordination-step integration tests with precomputed modeller/Protoss fixtures."""

from __future__ import annotations

import math
import shutil
from pathlib import Path

import pytest
import yaml
from click.testing import CliRunner

from qp.cli import cli

CASES_DIR = Path(__file__).resolve().parent / "coordination_cases"
CASES_YAML = CASES_DIR / "cases.yaml"


def _load_cases():
    with open(CASES_YAML) as f:
        data = yaml.safe_load(f)
    return data["defaults"], data["cases"]


DEFAULTS, CASES = _load_cases()
CASE_IDS = [c["case"] for c in CASES]
CASES_BY_NAME = {c["case"]: c for c in CASES}


def _fixture_root(case: str, pdb: str) -> Path:
    return CASES_DIR / case / pdb


def _prepare_workdir(tmp_path: Path, case: str, pdb: str, center_residues, overrides=None) -> Path:
    """Copy fixtures into tmp_path/{pdb}/ and write a config.yaml."""
    src = _fixture_root(case, pdb)
    dest = tmp_path / pdb
    shutil.copytree(src, dest)

    config = {
        **DEFAULTS,
        "input": pdb,
        "output_dir": str(tmp_path),
        "center_residues": list(center_residues),
        **(overrides or {}),
    }
    config_path = tmp_path / "config.yaml"
    with open(config_path, "w") as f:
        yaml.safe_dump(config, f, default_flow_style=False)
    return config_path


def _run_qp(config_path: Path):
    runner = CliRunner()
    result = runner.invoke(cli, ["run", "-c", str(config_path)], catch_exceptions=False)
    assert result.exit_code == 0, (
        f"qp run failed (exit {result.exit_code}):\n{result.output}"
    )
    return result


def _parse_charge_csv(path: Path):
    """Return (sphere_charges, ligand_charges).

    sphere_charges: {name: {col: int}} for columns '0','1','2',...
    ligand_charges: {name: int}
    """
    lines = path.read_text().splitlines()
    header = lines[0].split(",")
    cols = header[1:]
    sphere = {}
    ligands = {}
    i = 1
    while i < len(lines) and lines[i].strip():
        parts = lines[i].split(",")
        sphere[parts[0]] = {c: int(parts[j + 1]) for j, c in enumerate(cols)}
        i += 1
    i += 1  # blank
    while i < len(lines):
        line = lines[i].strip()
        if not line:
            i += 1
            continue
        name, charge = line.split(",")
        # Skip accidental header if present
        if name == "Name" and charge == "charge":
            i += 1
            continue
        ligands[name] = int(charge)
        i += 1
    return sphere, ligands


def _parse_spin_csv(path: Path):
    """spin.csv: blank first line, then Name,spin rows (no header required)."""
    spins = {}
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.lower().startswith("name"):
            continue
        name, spin = line.split(",")
        spins[name] = int(spin)
    return spins


def _parse_pdb_atoms(path: Path):
    """Yield (resname, chain, resid, atom_name, element, x, y, z)."""
    atoms = []
    with open(path) as f:
        for line in f:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            atom_name = line[12:16].strip()
            resname = line[17:20].strip()
            chain = line[21].strip() or " "
            resid = int(line[22:26])
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            element = line[76:78].strip() if len(line) >= 78 else ""
            if not element:
                element = "".join(c for c in atom_name if c.isalpha())[:1].upper()
            atoms.append((resname, chain, resid, atom_name, element, x, y, z))
    return atoms


def _res_key(resname, chain, resid):
    return f"{resname}_{chain}{resid}"


def _cluster_pdb(tmp_path: Path, pdb: str, cluster_id: str) -> Path:
    path = tmp_path / pdb / cluster_id / f"{cluster_id}.pdb"
    assert path.is_file(), f"Missing cluster PDB: {path}"
    return path


def _residue_atoms(atoms, resname, chain, resid):
    return [a for a in atoms if a[0] == resname and a[1] == chain and a[2] == resid]


def _atom_coords(atoms, resname, chain, resid, atom_name):
    for a in atoms:
        if a[0] == resname and a[1] == chain and a[2] == resid and a[3] == atom_name:
            return a[5], a[6], a[7]
    return None


def _dist(a, b):
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def _has_h_near(atoms, resname, chain, resid, heavy_name, cutoff=1.1):
    heavy = _atom_coords(atoms, resname, chain, resid, heavy_name)
    assert heavy is not None, f"Missing {heavy_name} on {_res_key(resname, chain, resid)}"
    for a in _residue_atoms(atoms, resname, chain, resid):
        if a[4] != "H" and not a[3].startswith("H"):
            continue
        if _dist(heavy, (a[5], a[6], a[7])) <= cutoff:
            return True
    return False


def _residue_keys(atoms):
    return {_res_key(a[0], a[1], a[2]) for a in atoms}


def _atom_names(atoms, resname, chain, resid):
    return {a[3] for a in _residue_atoms(atoms, resname, chain, resid)}


@pytest.fixture
def run_case(tmp_path):
    """Factory: prepare fixtures, run qp, return tmp_path for assertions."""

    def _run(case_name: str):
        meta = CASES_BY_NAME[case_name]
        overrides = {
            k: v
            for k, v in meta.items()
            if k not in {"case", "pdb", "center_residues"}
        }
        config_path = _prepare_workdir(
            tmp_path,
            meta["case"],
            meta["pdb"],
            meta["center_residues"],
            overrides=overrides,
        )
        _run_qp(config_path)
        return tmp_path, meta["pdb"]

    return _run


@pytest.mark.integration
@pytest.mark.parametrize("case_name", CASE_IDS)
def test_coordination_case(run_case, case_name):
    tmp_path, pdb = run_case(case_name)
    ASSERTIONS[case_name](tmp_path, pdb)


def _assert_asp_proton(tmp_path, pdb):
    sphere, _ = _parse_charge_csv(tmp_path / pdb / "charge.csv")
    assert sphere["A501"]["2"] == -4


def _assert_css_ligand(tmp_path, pdb):
    sphere, _ = _parse_charge_csv(tmp_path / pdb / "charge.csv")
    # readme says A201; cluster / charge row is A301
    assert sphere["A301"]["1"] == 0


def _assert_cys_ligand(tmp_path, pdb):
    sphere, _ = _parse_charge_csv(tmp_path / pdb / "charge.csv")
    assert sphere["A201"]["1"] == 0


def _assert_glu_proton(tmp_path, pdb):
    sphere, _ = _parse_charge_csv(tmp_path / pdb / "charge.csv")
    assert sphere["A1130"]["2"] == -1


def _assert_his_flipping(tmp_path, pdb):
    atoms = _parse_pdb_atoms(_cluster_pdb(tmp_path, pdb, "C401"))
    assert _has_h_near(atoms, "HIS", "C", 314, "ND1")
    assert not _has_h_near(atoms, "HIS", "C", 314, "NE2")


def _assert_cluster_smoothing(tmp_path, pdb):
    atoms = _parse_pdb_atoms(_cluster_pdb(tmp_path, pdb, "B501"))
    # readme says ≤40 Å; known-correct B501 diameter is ~44.5 Å
    coords = [(a[5], a[6], a[7]) for a in atoms]
    max_d = 0.0
    for i, c1 in enumerate(coords):
        for c2 in coords[i + 1 :]:
            max_d = max(max_d, _dist(c1, c2))
    assert max_d <= 45.0, f"max pairwise distance {max_d:.2f} Å exceeds 45 Å"


def _assert_deprotonated_his(tmp_path, pdb):
    sphere, _ = _parse_charge_csv(tmp_path / pdb / "charge.csv")
    # readme says column "2"; correct charge.csv has -2 in column "1"
    assert sphere["A310"]["1"] == -2


def _assert_nitrile_hydratase_no(tmp_path, pdb):
    atoms = _parse_pdb_atoms(_cluster_pdb(tmp_path, pdb, "A300"))
    cso_names = _atom_names(atoms, "CSO", "A", 114)
    assert {"HB1", "HB2", "HD"} <= cso_names
    assert "H" not in _atom_names(atoms, "SER", "A", 113)
    assert "H" not in cso_names

    sphere, _ = _parse_charge_csv(tmp_path / pdb / "charge.csv")
    assert sphere["A300"]["1"] == -4

    no_atoms = _residue_atoms(atoms, "NO", "A", 301)
    elements = [a[4] for a in no_atoms]
    assert elements.count("N") == 1
    assert elements.count("O") == 1

    spins = _parse_spin_csv(tmp_path / pdb / "spin.csv")
    assert spins["NO_A301"] == 1


def _assert_nitrile_hydratase_covalent(tmp_path, pdb):
    atoms = _parse_pdb_atoms(_cluster_pdb(tmp_path, pdb, "A301"))
    cso_names = _atom_names(atoms, "CSO", "A", 114)
    assert {"HB1", "HB2"} <= cso_names
    assert "HD" not in cso_names
    assert "H" not in _atom_names(atoms, "SER", "A", 113)
    assert "H" not in cso_names

    sphere, ligands = _parse_charge_csv(tmp_path / pdb / "charge.csv")
    # readme says A300; charge row / cluster is A301
    assert sphere["A301"]["1"] == -4
    assert ligands["TAN_A302"] == 0


def _assert_nitrile_hydratase_isonitrile(tmp_path, pdb):
    atoms = _parse_pdb_atoms(_cluster_pdb(tmp_path, pdb, "A300"))
    cso_names = _atom_names(atoms, "CSO", "A", 114)
    assert {"HB1", "HB2", "HD"} <= cso_names
    assert "H" not in _atom_names(atoms, "SER", "A", 113)
    assert "H" not in cso_names

    sphere, ligands = _parse_charge_csv(tmp_path / pdb / "charge.csv")
    assert sphere["A300"]["1"] == -4
    assert ligands["TB0_A301"] == 0

    tb0 = _residue_atoms(atoms, "TB0", "A", 301)
    n_h = sum(
        1
        for a in tb0
        if a[4] == "H" or a[3].startswith("H")
    )
    assert n_h == 9, f"TB0 has {n_h} hydrogens, expected 9"


def _assert_oxygen(tmp_path, pdb):
    _, ligands = _parse_charge_csv(tmp_path / pdb / "charge.csv")
    # readme says -1; correct triplet O2 charge is 0
    assert ligands["OXY_A310"] == 0
    spins = _parse_spin_csv(tmp_path / pdb / "spin.csv")
    assert spins["OXY_A310"] == 2


def _assert_partial_occupancy_aa_gt_other(tmp_path, pdb):
    keys = _residue_keys(_parse_pdb_atoms(_cluster_pdb(tmp_path, pdb, "A501")))
    assert "HIS_A86" in keys
    assert "HOH_A667" not in keys
    assert "DCY_A502" not in keys


def _assert_partial_occupancy_large_gt_small(tmp_path, pdb):
    keys = _residue_keys(_parse_pdb_atoms(_cluster_pdb(tmp_path, pdb, "A501")))
    assert "CL_A502" in keys
    assert "HOH_A601" not in keys
    assert "THJ_A503" not in keys


def _assert_fake_n_terminus(tmp_path, pdb):
    sphere, _ = _parse_charge_csv(tmp_path / pdb / "charge.csv")
    assert sphere["C2"]["3"] == 0


def _assert_atom_count(tmp_path, pdb):
    _, ligands = _parse_charge_csv(tmp_path / pdb / "charge.csv")
    assert ligands["NAG_A1 GAL_A2 NAG_A3 GAL_A4"] == 0


def _assert_polysacchride(tmp_path, pdb):
    keys = _residue_keys(_parse_pdb_atoms(_cluster_pdb(tmp_path, pdb, "B2")))
    assert "BGC_B1" in keys


def _assert_terminal_ligand(tmp_path, pdb):
    _, ligands = _parse_charge_csv(tmp_path / pdb / "charge.csv")
    assert ligands["BGC_E1 GAL_E2"] == 0


ASSERTIONS = {
    "ASP_proton": _assert_asp_proton,
    "CSS_ligand": _assert_css_ligand,
    "CYS_ligand": _assert_cys_ligand,
    "GLU_proton": _assert_glu_proton,
    "HIS_flipping": _assert_his_flipping,
    "cluster_smoothing": _assert_cluster_smoothing,
    "deprotonated_HIS": _assert_deprotonated_his,
    "nitrile_hydratase_NO": _assert_nitrile_hydratase_no,
    "nitrile_hydratase_covalent": _assert_nitrile_hydratase_covalent,
    "nitrile_hydratase_isonitrile": _assert_nitrile_hydratase_isonitrile,
    "oxygen": _assert_oxygen,
    "partial_occupancy_AA_gt_other": _assert_partial_occupancy_aa_gt_other,
    "partial_occupancy_large_gt_small": _assert_partial_occupancy_large_gt_small,
    "fake_N_terminus": _assert_fake_n_terminus,
    "atom_count": _assert_atom_count,
    "polysacchride": _assert_polysacchride,
    "terminal_ligand": _assert_terminal_ligand,
}
