"""Command-line interface (CLI) entry point.

**Usage** ``qpdb [OPTIONS]``

Options:

==================== =====
-i                   Input PDB code, PDB file, or batch file [required]
-o                   Output directory  [required]
-m, -\-modeller      Use MODELLER to add missing loops 
-p, -\-protoss       Use Protoss to add hydrogens
-c, -\-coordination  Select the first, second, etc. coordination spheres
-s, -\-skip          Skips rerunning MODELLER or Protoss if an
                     existing output file is found
==================== =====

"""

# Print first to welcome the user while it waits to load the modules
print("\n.-------------------------------.")
print("| WELCOME TO THE QUANTUMPDB CLI |")
print(".-------------------------------.")
print("Default programmed actions for the quantumPDB package.")
print("GitHub: https://github.com/davidkastner/quantumpdb")
print("Documenation: https://quantumpdb.readthedocs.io\n")

import os
import click
from qp.checks import fetch_pdb
from Bio.PDB.PDBExceptions import PDBIOException

@click.command()
@click.option("-i", required=True, multiple=True, help="Input PDB code, PDB file, or batch file")
@click.option("-o", required=True, type=click.Path(file_okay=False), help="Output directory")
@click.option("--modeller", "-m", is_flag=True, help="Use MODELLER to add missing loops")
@click.option("--protoss", "-p", is_flag=True, help="Use Protoss to add hydrogens")
@click.option("--coordination", "-c", is_flag=True, help="Select the first, second, etc. coordination spheres")
@click.option("--skip", "-s", type=click.Choice(["modeller", "protoss", "all"]), is_flag=False, flag_value="all",
                              help="Skips rerunning MODELLER or Protoss if an existing output file is found")
def cli(
    i,
    o,
    modeller,
    protoss,
    coordination,
    skip
    ):
    """
    The overall command-line interface (CLI) entry point.
    The CLI interacts with the rest of the package.

    A complete reference of quantumPDB functionality.
    This is advantagous because it quickly introduces the package.
    Specifically, to the complete scope of available functionality.
    It also improves long-term maintainability and readability.

    """
    o = os.path.abspath(o)

    # Store input PDBs as a tuple of 
    #   parsed ID (PDB code or filename)
    #   path to PDB file (existing or to download)
    pdb_all = []
    for p in i:
        if os.path.isfile(p):
            pdb, ext = os.path.splitext(os.path.basename(p))
            pdb = pdb.replace(".", "_")
            if ext == ".pdb":
                pdb_all.append((pdb, p))
            else:
                with open(p, "r") as f:
                    pdb_all.extend([(pdb, f"{o}/{pdb}/{pdb}.pdb") 
                                    for pdb in f.read().splitlines()])
        else:
            pdb_all.append((p, f"{o}/{p}/{p}.pdb"))

    if modeller:
        from qp.structure import missing_loops
        click.echo("MODELLER parameters:")
        optimize = int(click.prompt("> Optimize\n   0: None\n   1: [Missing]\n   2: All\n ", 
                                    type=click.Choice(["0", "1", "2"]), default="1", show_default=False))
        click.echo("")

    if coordination:
        from qp.cluster import coordination_spheres
        click.echo("Coordination sphere parameters:")
        metals = click.prompt("> Active site metals", default="FE FE2").split(" ")
        limit = click.prompt("> Number of spheres", default=2)
        ligands = click.prompt("> Additional ligands [AAs and waters]", default=[], show_default=False)
        capping = int(click.prompt("> Capping (requires Protoss)\n   0: [None]\n   1: H\n   2: ACE/NME\n ", 
                                   type=click.Choice(["0", "1", "2"]), default="0", show_default=False))
        charge = click.confirm("> Compute charges (requires Protoss)")
        if capping or charge:
            protoss = True
        click.echo("")

    if protoss:
        from qp.structure import add_hydrogens

    err = {
        "PDB": [],
        "Protoss": [],
        "Coordination sphere": [],
        "Other": []
    }
    for pdb, path in pdb_all:
        click.secho(pdb, bold=True)

        # Skips fetching if PDB file exists
        if not os.path.isfile(path):
            click.echo(f"> Fetching PDB file")
            try:
                fetch_pdb(pdb, path)
            except ValueError:
                click.secho("Error fetching PDB file\n", italic=True, fg="red")
                err["PDB"].append(pdb)
                continue

        mod_path = f"{o}/{pdb}/{pdb}_modeller.pdb"
        if modeller:
            if skip in ["modeller", "all"] and os.path.isfile(mod_path):
                click.echo("> MODELLER file found")
            else:
                click.echo("> Building model")
                residues = missing_loops.get_residues(path)
                ali_path = f"{o}/{pdb}/{pdb}.ali"
                missing_loops.write_alignment(residues, pdb, path, ali_path)
                missing_loops.build_model(residues, pdb, ali_path, mod_path, optimize)

        prot_path = f"{o}/{pdb}/Protoss"
        if protoss: 
            if skip in ["protoss", "all"] and os.path.isfile(f"{prot_path}/{pdb}_protoss.pdb"):
                click.echo("> Protoss file found")
            else:
                click.echo("> Running Protoss")
                if modeller:
                    path = mod_path
                try:
                    pid = add_hydrogens.upload(path)
                except ValueError:
                    click.secho("Error uploading PDB file\n", italic=True, fg="red")
                    # Occurs when uploading a PDB file > 4MB
                    # TODO retry with parsed PDB code? Will raise Bio.PDB error if  
                    #      number of atoms in output exceeds 99999
                    err["Protoss"].append(pdb)
                    continue
                job = add_hydrogens.submit(pid)
                add_hydrogens.download(job, f"{prot_path}/{pdb}_protoss.pdb", "protein")
                add_hydrogens.download(job, f"{prot_path}/{pdb}_ligands.sdf", "ligands")
                add_hydrogens.download(job, f"{prot_path}/{pdb}_log.txt", "log")
        
        if coordination:
            click.echo("> Extracting clusters")
            if modeller:
                path = mod_path
            if protoss:
                path = f"{prot_path}/{pdb}_protoss.pdb"
            charge_path = f"{o}/{pdb}/{pdb}_charge.csv" if charge else None

            try:
                if protoss:
                    add_hydrogens.adjust_active_sites(path, metals)
                clusters = coordination_spheres.extract_clusters(path, f"{o}/{pdb}", metals, limit, ligands, capping, charge_path)
            except (ValueError, PDBIOException):
                click.secho("Residue or atom limit exceeded\n", italic=True, fg="red")
                err["Other"].append(pdb)
                continue
            except KeyError:
                click.secho("Missing template atoms for capping\n", italic=True, fg="red")
                err["Coordination sphere"].append(pdb)
                continue

            if charge:
                ligand_charge = add_hydrogens.compute_charge(f"{prot_path}/{pdb}_ligands.sdf")
                with open(charge_path, "a") as f:
                    f.write("\n")
                    for k, v in ligand_charge.items():
                        f.write(f"{k},{v}\n")
        
        click.echo("")

    for k, v in err.items():
        if v:
            click.echo(click.style(k + " errors: ", bold=True, fg="red") + ", ".join(v))

if __name__ == "__main__":
    # Run the command-line interface when this script is executed
    cli()
