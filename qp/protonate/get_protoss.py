"""Add hydrogens using Protoss

**Usage**

#. Submitting existing or custom PDB file::

    >>> from qp.protonate import get_protoss
    >>> pid = get_protoss.upload("path/to/PDB.pdb")
    >>> job = get_protoss.submit(pid)
    >>> get_protoss.download(job, "path/to/OUT.pdb")

#. Submitting PDB code::

    >>> from qp.protonate import get_protoss
    >>> pdb = "1dry"
    >>> job = get_protoss.submit(pdb)
    >>> get_protoss.download(job, "path/to/OUT.pdb")
    >>> get_protoss.download(job, "path/to/OUT.sdf", "ligands")

Protoss automatically removes alternative conformations and overlapping entries.
Under ProteinsPlus API v1, download the log file (``key="log"`` in
``get_protoss.download``) to see affected atoms. API v2 no longer exposes
that clash log; ``key="log"`` writes an empty placeholder.

Some metal-coordinating residues may be incorrectly protonated. Use
``get_protoss.adjust_activesites(path, metals)`` with the metal IDs to deprotonate
these residues. 
"""

import os
import json
import time
import requests
from Bio.PDB import PDBParser, PDBIO

API_BASE = "https://proteins.plus/api/v2/"
UPLOAD_URL = API_BASE + "molecule_handler/upload/"
UPLOAD_JOBS_URL = API_BASE + "molecule_handler/upload/jobs/"
PROTEINS_URL = API_BASE + "molecule_handler/proteins/"
LIGANDS_URL = API_BASE + "molecule_handler/ligands/"
PROTOSS_URL = API_BASE + "protoss/"
PROTOSS_JOBS_URL = API_BASE + "protoss/jobs/"


def _retry_settings(default_retries=5, default_delay=60):
    """Return (retries, delay) from env, falling back to defaults.

    ``PROTOSS_UPLOAD_RETRIES`` / ``PROTOSS_UPLOAD_RETRY_DELAY`` let CI use a
    short retry budget so flaky ProteinsPlus outages become fast skips instead
    of multi-minute failures.
    """
    retries = int(os.environ.get("PROTOSS_UPLOAD_RETRIES", default_retries))
    delay = int(os.environ.get("PROTOSS_UPLOAD_RETRY_DELAY", default_delay))
    return max(1, retries), max(0, delay)



def _poll_job(job_id, jobs_url, timeout=600, poll_interval=2):
    """Poll a ProteinsPlus v2 job until it leaves pending/running."""
    deadline = time.time() + timeout
    url = f"{jobs_url}{job_id}/"
    while True:
        r = requests.get(url, timeout=60)
        r.raise_for_status()
        job = r.json()
        status = job.get("status")
        if status not in ("pending", "running"):
            return job
        if time.time() > deadline:
            raise TimeoutError(f"Job {job_id} did not finish within {timeout}s (last status={status})")
        time.sleep(poll_interval)


def _is_pdb_code(pid):
    """Return True if pid looks like a 4-character PDB code (not a UUID)."""
    return isinstance(pid, str) and len(pid) == 4 and "-" not in pid


class _UploadedProtein(str):
    """Protein UUID that retains the source file for ligand-preserving Protoss."""

    def __new__(cls, protein_id, source_path):
        value = super().__new__(cls, protein_id)
        value.source_path = source_path
        return value


def upload(path):
    """Upload a PDB file to the ProteinsPlus web server (API v2).

    Parameters
    ----------
    path : str
        Path to the PDB file to upload.

    Returns
    -------
    str
        ProteinsPlus protein UUID for the uploaded structure.

    Raises
    ------
    ValueError
        If the upload is rejected or preprocessing fails.
    KeyError
        If the upload fails after all retry attempts.
    """
    retries, delay = _retry_settings()

    for attempt in range(retries):
        try:
            
            with open(path, "rb") as fh:
                pp = requests.post(
                    UPLOAD_URL,
                    files={"protein_file": (os.path.basename(path), fh)},
                    timeout=180,
                )
            
            if pp.status_code == 400:
                raise ValueError("Bad request")
            if pp.status_code not in (200, 202):
                raise ValueError(f"Upload failed with HTTP {pp.status_code}")

            job_id = pp.json()["job_id"]
            job = _poll_job(job_id, UPLOAD_JOBS_URL)
            
            if job.get("status") != "success" or not job.get("output_protein"):
                raise ValueError(job.get("error") or "Upload preprocessing failed")
            return _UploadedProtein(job["output_protein"], path)
        except (KeyError, json.JSONDecodeError, requests.RequestException, TimeoutError) as e:
            print(f"> Upload error ({type(e).__name__}). Retrying in {delay} seconds...")
            
            time.sleep(delay)

    raise KeyError(f"> Failed to upload the file and retrieve protein id after {retries} attempts.")


def submit(pid):
    """Submit a PDB code or ProteinsPlus protein UUID to the Protoss web API (v2).

    Parameters
    ----------
    pid : str
        Four-character PDB code or ProteinsPlus protein UUID from :func:`upload`.

    Returns
    -------
    str
        URL of the Protoss job location for status polling.

    Raises
    ------
    ValueError
        If the PDB code is invalid (server returns 400) or submission fails.
    KeyError
        If the submission fails after all retry attempts.
    """
    retries, delay = _retry_settings()

    for attempt in range(retries):
        try:
            protein_id = pid
            if _is_pdb_code(pid):
                up = requests.post(UPLOAD_URL, data={"pdb_code": pid.lower()}, timeout=60)
                if up.status_code == 400:
                    raise ValueError("Invalid PDB code")
                if up.status_code not in (200, 202):
                    raise ValueError(f"PDB code upload failed with HTTP {up.status_code}")
                up_job = _poll_job(up.json()["job_id"], UPLOAD_JOBS_URL)
                if up_job.get("status") != "success" or not up_job.get("output_protein"):
                    raise ValueError(up_job.get("error") or "Invalid PDB code")
                protein_id = up_job["output_protein"]

            ligand_protoss = requests.post(
                PROTOSS_URL,
                data={"protein_id": protein_id},
                headers={"Accept": "application/json"},
                timeout=60,
            )
            structure_protoss = ligand_protoss
            if hasattr(pid, "source_path"):
                with open(pid.source_path, "rb") as fh:
                    structure_protoss = requests.post(
                        PROTOSS_URL,
                        files={"protein_file": (os.path.basename(pid.source_path), fh)},
                        headers={"Accept": "application/json"},
                        timeout=180,
                    )
            
            if structure_protoss.status_code == 400 or ligand_protoss.status_code == 400:
                raise ValueError("Invalid protein id / PDB code")
            if structure_protoss.status_code not in (200, 202):
                raise ValueError(f"Protoss structure submit failed with HTTP {structure_protoss.status_code}")
            if ligand_protoss.status_code not in (200, 202):
                raise ValueError(f"Protoss ligand submit failed with HTTP {ligand_protoss.status_code}")
            structure_job_id = structure_protoss.json()["job_id"]
            ligand_job_id = ligand_protoss.json()["job_id"]
            return {
                "protein": f"{PROTOSS_JOBS_URL}{structure_job_id}/",
                "ligands": f"{PROTOSS_JOBS_URL}{ligand_job_id}/",
            }
        except (KeyError, json.JSONDecodeError, requests.RequestException, TimeoutError) as e:
            print(f"> Submit error ({type(e).__name__}). Retrying in {delay} seconds...")
            time.sleep(delay)

    raise KeyError(f"> Failed to submit and retrieve Protoss job after {retries} attempts.")



def download(job, out, key="protein"):
    """Download a Protoss output file (API v2).

    Polls the Protoss job URL until completion, then downloads the
    requested output file.

    Parameters
    ----------
    job : str
        URL of the Protoss job location from :func:`submit`.
    out : str
        Path to the output file (directory created if needed).
    key : str, optional
        File type to download: ``'protein'`` (protonated PDB),
        ``'ligand'`` / ``'ligands'`` (ligand SDF), or ``'log'`` (processing log).
        Default is ``'protein'``.

    Raises
    ------
    KeyError
        If the download fails after 5 retry attempts.
    """
    retries = 5
    delay = 60  # seconds
    if isinstance(job, dict):
        job_url = job["ligands"] if key in ("ligand", "ligands") else job["protein"]
    else:
        job_url = job
    job_id = job_url.rstrip("/").split("/")[-1]

    for attempt in range(retries):
        try:
            job_data = _poll_job(job_id, PROTOSS_JOBS_URL)
            if job_data.get("status") != "success":
                raise ValueError(job_data.get("error") or f"Protoss job status={job_data.get('status')}")

            output_protein = job_data.get("output_protein")
            if not output_protein:
                raise KeyError("output_protein")

            os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)

            if key == "protein":
                prot = requests.get(f"{PROTEINS_URL}{output_protein}/", timeout=120)
                prot.raise_for_status()
                content = prot.json()["file_string"]
            elif key in ("ligand", "ligands"):
                prot = requests.get(f"{PROTEINS_URL}{output_protein}/", timeout=120)
                prot.raise_for_status()
                parts = []
                for lig_id in prot.json().get("ligand_set") or []:
                    lig = requests.get(f"{LIGANDS_URL}{lig_id}/", timeout=60)
                    lig.raise_for_status()
                    # API v2 returns each ligand ending in "$$$$" without a
                    # trailing newline; concatenate as standard multi-mol SDF.
                    mol = lig.json()["file_string"].rstrip("\n")
                    if mol.endswith("$$$$"):
                        mol = mol[:-4].rstrip("\n")
                    parts.append(mol + "\n$$$$\n")
                content = "".join(parts)
            elif key == "log":
                # ProteinsPlus API v2 no longer exposes a Protoss clash log.
                content = ""
            else:
                raise KeyError(key)

            
            with open(out, "w") as f:
                f.write(content)
            return
        except (KeyError, json.JSONDecodeError, requests.RequestException, ValueError, TimeoutError) as e:
            
            time.sleep(delay)
            continue

    raise KeyError(f"> Failed to download the file with key '{key}' after {retries} attempts.")


def repair_ligands(path, orig):
    """Repair ligands that Protoss renamed to ``MOL``.

    Protoss sometimes replaces unrecognized ligand residues with a generic
    ``MOL`` label. This function restores the original residue names and
    structures by matching them back from the pre-Protoss PDB, then
    reassigns hydrogen atoms to the closest heavy atoms.

    Parameters
    ----------
    path : str
        Path to the Protoss output PDB file (modified in place).
    orig : str
        Path to the original (pre-Protoss) PDB file.
    """
    parser = PDBParser(QUIET=True)
    prot_structure = parser.get_structure("Prot", path)
    orig_structure = parser.get_structure("Orig", orig)

    for res in prot_structure[0].get_residues():
        if res.get_resname() == "MOL":
            resid = res.get_id()
            chain = res.get_parent()
            chain.detach_child(resid)

            missing = []
            found = False
            for r in orig_structure[0][chain.get_id()].get_residues():
                if r.get_id()[1] == resid[1]:
                    found = True
                if found:
                    if r.get_id() not in chain:
                        chain.add(r)
                        missing.append(r)
                    else:
                        break

            for r in missing:
                for a in r.get_unpacked_list():
                    if a.element == "H":
                        r.detach_child(atom.get_id())
            
            for atom in res.get_unpacked_list():
                if atom.element != "H":
                    continue
                closest = None
                for r in missing:
                    for a in r.get_unpacked_list():
                        if closest is None or atom - a < atom - closest:
                            closest = a
                closest.get_parent().add(atom)

    io = PDBIO()
    io.set_structure(prot_structure)
    io.save(path)
