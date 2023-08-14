"""Command-line interface (CLI) entry point."""

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

@click.command()
@click.option("-i", required=True, multiple=True, help="Input PDB code, PDB file, or batch file")
@click.option("-o", required=True, type=click.Path(file_okay=False), help="Output directory")
@click.option("--modeller", "-m", is_flag=True, help="Use MODELLER to add missing loops")
@click.option("--protoss", "-p", is_flag=True, help="Use Protoss to add hydrogens")
@click.option("--coordination", "-c", is_flag=True, help="Select the first, second, etc. coordination spheres")
@click.option("--skip", "-s", type=click.Choice(["modeller", "protoss", "all"]), help="Skips rerunning MODELLER or Protoss if an existing output file is found")
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
            if ext == "pdb":
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
        if capping:
            protoss = True
        click.echo("")

    if protoss:
        from qp.structure import protoss

    for pdb, path in pdb_all:
        click.secho(pdb, bold=True)

        # Skips fetching if PDB file exists
        if not os.path.isfile(path):
            click.echo(f"> Fetching PDB file")
            try:
                fetch_pdb(pdb, path)
            except ValueError:
                click.secho("Error fetching PDB file", italic=True, fg="red")

        mod_path = f"{o}/{pdb}/{pdb}_modeller.pdb"
        if modeller:
            if skip in ["modeller", "all"] and os.path.isfile(mod_path):
                click.echo("> MODELLER file found")
            else:
                click.echo("> Building model")
                residues = missing_loops.get_residues(path)
                ali_path = f"{o}/{pdb}/{pdb}.ali"
                missing_loops.write_alignment(residues, pdb, path, ali_path)
                missing_loops.build_model(residues, ali_path, pdb, mod_path, optimize)

        prot_path = f"{o}/{pdb}/{pdb}_protoss.pdb"
        if protoss: 
            if skip in ["protoss", "all"] and os.path.isfile(prot_path):
                click.echo("> Protoss file found")
            else:
                click.echo("> Running Protoss")
                if modeller:
                    path = mod_path
                try:
                    pid = protoss.upload(path)
                except ValueError:
                    click.secho("Error uploading PDB file", italic=True, fg="red")
                    # Occurs when uploading a PDB file > 4MB
                    # TODO retry with parsed PDB code? Will raise Bio.PDB error if  
                    #      number of atoms in output exceeds 99999
                    continue
                job = protoss.submit(pid)
                protoss.download(job, prot_path)
        
        if coordination:
            click.echo("> Extracting clusters")
            if modeller:
                path = mod_path
            if protoss:
                path = prot_path
                protoss.adjust_active_sites(path, metals)
            try:
                coordination_spheres.extract_clusters(path, f"{o}/{pdb}", metals, limit, ligands, capping)
            except KeyError:
                click.secho("Missing template atoms for capping", italic=True, fg="red")


if __name__ == "__main__":
    # Run the command-line interface when this script is executed
    cli()
