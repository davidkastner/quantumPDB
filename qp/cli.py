"""Command-line interface (CLI) entry point.

**Usage** ``qp [OPTIONS]``

Options:

==================== =====
-i                   Input PDB code, PDB file, or batch file [required]
-o                   Output directory  [required]
-m, --modeller       Use MODELLER to add missing loops 
-p, --protoss        Use Protoss to add hydrogens
-c, --coordination   Select the first, second, etc. coordination spheres
-s, --skip           Skips rerunning MODELLER or Protoss if an
                     existing output file is found
==================== =====

"""

import click

def welcome():
    """Print first to welcome the user while it waits to load the modules"""
    
    print("\n        ╔═══════════════════════════════╗")
    print("        ║ .---------------------------. ║")
    print("        ║ |         __________        | ║")
    print("        ║ |       / ____/\____ \      | ║")
    print("        ║ |      < <_|  ||  |_> >     | ║")
    print("        ║ |       \__   ||   __/      | ║")
    print("        ║ |          |__||__|         | ║")
    print("        ║ |                           | ║")
    print("        ║ |   WELCOME TO QUANTUMPDB   | ║")
    print("        ║ '---------------------------' ║")
    print("        ╚═══════════════════════════════╝\n")

    print("GitHub: https://github.com/davidkastner/quantumpdb")
    print("Documenation: https://quantumpdb.readthedocs.io")
    print("• Flags: -m (modeller) -p (protoss) -c (spheres) -s (use previous modeller/protoss)")
    print("• Generate structures: qp run -i pdbs.txt -o datasets/ -m -p -c")
    print("• Submit QM: qp submit -j\n")

# Welcome even if no flags
welcome()

@click.group()
def cli():
    """CLI entry point"""
    pass

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
                "> Optimize select residues\n   0: None\n   1: [Missing]\n   2: All\n ",
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
        limit = click.prompt("> Number of spheres", default=3)
        first_sphere_radius = click.prompt("> Radius of first spheres", type=float, default=4.0)
        ligands = click.prompt("> Additional ligands [unnatural AAs] in outer coordination spheres", default=[], show_default=False)
        capping = int(
            click.prompt(
                "> Capping (requires Protoss)\n   0: None\n   1: [H]\n   2: ACE/NME\n ",
                type=click.Choice(["0", "1", "2"]),
                default="1",
                show_default=False,
            )
        )

        # Prompt user for their preferred cluster smoothing method
        choice = click.prompt(
            "> Smoothing method\n   0: Box plot\n   1: DBSCAN\n   2: [Dummy Atom]\n   3: None\n",
            type=click.Choice(["0", "1", "2", "3"]),
            default="2",
            show_default=False,)
        smooth_options = {"0": {}, "1": {"eps": 6, "min_samples": 3}, "2": {"mean_distance": 3}, "3": {}}
        smooth_method_options = {"0": "box_plot", "1": "dbscan", "2": "dummy_atom", "3": False}
        smooth_params = smooth_options[choice]
        smooth_method = smooth_method_options[choice]

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
                # Remove trailing missing residues from the ends of all chains
                residues = missing_loops.clean_termini(residues)

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

            # try:
            if protoss:
                add_hydrogens.adjust_activesites(path, metals)

            clusters = coordination_spheres.extract_clusters(
                path, f"{o}/{pdb}", metals,
                limit, ligands, capping, charge, count, xyz, first_sphere_radius,
                smooth_method=smooth_method,
                **smooth_params
            )
            # except (ValueError, PDBIOException):  # TODO add custom exceptions
            #     click.secho("Residue or atom limit exceeded\n", italic=True, fg="red")
            #     err["Other"].append(pdb)
            #     continue
            # except KeyError as e:
            #     raise e
            #     click.secho(f"Missing template atoms for capping {e}\n", italic=True, fg="red")
            #     err["Coordination sphere"].append(pdb)
            #     continue

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
@click.option("--failure_checkup", "-f", is_flag=True, help="Find failed structures")
def submit(job_manager,
           failure_checkup,):
    """Handles the submission of jobs for the quantumPDB."""

    if job_manager:
        from qp.manager import job_manager
        
        job_count = click.prompt("> Jobs count to be submitted in this batch", default=80)

        # Single point or geometry optimization
        choice_map = {"0": False, "1": True}
        calculation_type = click.prompt(
            "> Calculation type:\n   0: [single point]\n   1: geometry optimization\n",
            type=click.Choice(["0", "1"]),
            default="0",
            show_default=False)
        minimization = choice_map[calculation_type]
        click.echo("")

        # Get the users perferred functional
        choice_map = {"0": "uwpbeh", "1": "ub3lyp", "2": "ugfn2xtb", "3": "gfn2xtb"}
        method_number = click.prompt(
            "> Requested functional:\n   0: [uwPBEh]\n   1: uB3LYP\n   2: uGFN2xTB\n   3: GFN2xTB\n ",
            type=click.Choice(["0", "1", "2", "3"]),
            default="0",
            show_default=False)
        method = choice_map[method_number]

        click.echo("")

        if method == "uwpbeh":
            basis = "lacvps_ecp"
            guess = "generate"
            gpus = 1
            memory = "8G"
        elif method == "ub3lyp":
            basis = "lacvps_ecp"
            guess = "generate"
            gpus = 1
            memory = "8G"
        elif method == "ugfn2xtb":
            basis = "gfn2xtb"
            guess = "hcore"
            gpus = 1
            memory = "8G"
        elif method == "gfn2xtb":
            basis = "gfn2xtb"
            guess = "hcore"
            gpus = 1
            memory = "8G"
        
        url = "https://docs.google.com/spreadsheets/d/1St_4YEKcWzrs7yS1GTehfAKabtGCJeuM0sC0u1JG8ZE/gviz/tq?tqx=out:csv&sheet=Sheet1"
        master_list_path = job_manager.get_master_list(url)
        job_manager.submit_jobs(job_count, master_list_path, minimization, basis, method, guess, gpus, memory)

    if failure_checkup:
        from qp.manager import failure_checkup
        qm_job_dir = input("What is the name of your QM job directory? ")
        failure_counts = failure_checkup.check_all_jobs(qm_job_dir)
        failure_checkup.plot_failures(failure_counts)


if __name__ == "__main__":
    # Run the command-line interface when this script is executed
    welcome()
    cli()
