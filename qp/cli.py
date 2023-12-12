"""Command-line interface (CLI) entry point.

**Usage** ``qp [OPTIONS]``

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

import click

@click.group()
def cli():
    # Print first to welcome the user while it waits to load the modules
    print("\n.-------------------------------.")
    print("| WELCOME TO THE QUANTUMPDB CLI |")
    print(".-------------------------------.")
    print("Default programmed actions for the quantumPDB package.")
    print("GitHub: https://github.com/davidkastner/quantumpdb")
    print("Documenation: https://quantumpdb.readthedocs.io\n")
    print("Run the workflow: qp run -i pdbs.txt -o datasets/ -m -p -c")
    print("Submit QM calculations: qp submit -j\n")


@cli.command()
@click.option("-i", required=True, multiple=True, help="Input PDB code, PDB file, or batch file")
@click.option("-o", required=True, type=click.Path(file_okay=False), help="Output directory")
@click.option("--modeller", "-m", is_flag=True, help="Use MODELLER to add missing loops")
@click.option("--protoss", "-p", is_flag=True, help="Use Protoss to add hydrogens")
@click.option("--coordination", "-c", is_flag=True, help="Select the first, second, etc. coordination spheres")
@click.option("--skip", "-s", type=click.Choice(["modeller", "protoss", "all"]), is_flag=False, flag_value="all",
                              help="Skips rerunning MODELLER or Protoss if an existing output file is found")
def run(i,
        o,
        modeller,
        protoss,
        coordination,
        skip,):
    """Generates quantumPDB structures and files."""
    
    import os
    from Bio.PDB.PDBExceptions import PDBIOException
    from qp.checks import fetch_pdb

    o = os.path.abspath(o)
    pdb_all = fetch_pdb.parse_input(i, o)

    if modeller:
        from qp.structure import missing_loops
        click.echo("MODELLER parameters:")
        optimize = int(
            click.prompt(
                "> Optimize\n   0: None\n   1: [Missing]\n   2: All\n ",
                type=click.Choice(["0", "1", "2"]),
                default="1",
                show_default=False,
            )
        )
        click.echo("")

    if coordination:
        from qp.cluster import coordination_spheres
        click.echo("Coordination sphere parameters:")
        metals = click.prompt("> Active site metals", default="FE FE2").split(" ")
        limit = click.prompt("> Number of spheres", default=2)
        ligands = click.prompt(
            "> Additional ligands [AAs]", default=[], show_default=False
        )
        capping = int(
            click.prompt(
                "> Capping (requires Protoss)\n   0: [None]\n   1: H\n   2: ACE/NME\n ",
                type=click.Choice(["0", "1", "2"]),
                default="0",
                show_default=False,
            )
        )
        charge = click.confirm("> Compute charges (requires Protoss)", default=True)
        count = click.confirm("> Count residues", default=True)
        xyz = click.confirm("> Write XYZ files", default=True)
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
                fetch_pdb.fetch_pdb(pdb, path)
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
                AA = missing_loops.define_residues()
                residues = missing_loops.get_residues(path, AA)
                # residues = missing_loops.clean_termini(residues, AA)
                ali_path = f"{o}/{pdb}/{pdb}.ali"
                missing_loops.write_alignment(residues, pdb, path, ali_path)
                # missing_loops.strip_ends_from_ali(ali_path)
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

            try:
                if protoss:
                    add_hydrogens.adjust_activesites(path, metals)
                    add_hydrogens.rename_nterminal(path)
                clusters = coordination_spheres.extract_clusters(
                    path, f"{o}/{pdb}", metals,
                    limit, ligands, capping, charge, count, xyz
                )
            except (ValueError, PDBIOException):  # TODO add custom exceptions
                click.secho("Residue or atom limit exceeded\n", italic=True, fg="red")
                err["Other"].append(pdb)
                continue
            except KeyError as e:
                click.secho(f"Missing template atoms for capping {e}\n", italic=True, fg="red")
                err["Coordination sphere"].append(pdb)
                continue

            if charge:
                ligand_charge = add_hydrogens.compute_charge(f"{prot_path}/{pdb}_ligands.sdf")
                with open(f"{o}/{pdb}/charge.csv", "a") as f:
                    f.write("\n")
                    for k, v in sorted(ligand_charge.items()):
                        f.write(f"{k},{v}\n")

        click.echo("")

    for k, v in err.items():
        if v:
            click.echo(click.style(k + " errors: ", bold=True, fg="red") + ", ".join(v))


@cli.command()
@click.option("--job_manager", "-j", is_flag=True, help="Submit qp generated QM jobs")
@click.option("--find_incomplete", "-f", is_flag=True, help="Find failed structures")
def submit(job_manager,
           find_incomplete,):
    """Handles the submission of jobs for the quantumPDB."""

    if job_manager:
        from qp.manager import job_manager
        
        job_count = click.prompt("> Jobs count to be submitted in this batch", default=80)
        
        # Get the users perferred functional
        choice_map = {"0": "uwpbeh", "1": "ugfn2xtb", "2": "gfn2xtb"}
        method_number = click.prompt(
            "> Requested functional:\n   0: [uwpbeh]\n   1: ugfn2xtb\n   2: gfn2xtb\n ",
            type=click.Choice(["0", "1", "2"]),
            default="0",
            show_default=False)
        method = choice_map[method_number]
        click.echo("")

        if method == "uwpbeh":
            basis = "lacvps_ecp"
            guess = "generate"
            constraint_freeze = ""
            gpus = 1
            memory = "8G"
        elif method == "ugfn2xtb":
            basis = "gfn2xtb"
            guess = "hcore"
            constraint_freeze = ""
            gpus = 1
            memory = "8G"
        elif method == "gfn2xtb":
            basis = "gfn2xtb"
            guess = "hcore"
            constraint_freeze = True
            gpus = 1
            memory = "8G"

        job_manager.submit_jobs(job_count, basis, method, guess,constraint_freeze, gpus, memory)

        if find_incomplete:
            from qp.manager import find_incomplete
            find_incomplete.find()



if __name__ == "__main__":
    # Run the command-line interface when this script is executed
    cli()
