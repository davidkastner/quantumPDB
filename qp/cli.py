"""Command-line interface (CLI) entry point.

**Usage** ``qp run --config path/to/config.yaml``
"""

import click
import yaml

def welcome():
    """Print first to welcome the user while it waits to load the modules"""
    print("\n")
    print("             ╔════════════════════════╗             ")
    print("             ║       __________       ║             ")
    print("             ║     / ____/\____ \     ║             ")
    print("             ║    < <_|  ||  |_> >    ║             ")
    print("             ║     \__   ||   __/     ║             ")
    print("             ║        |__||__|        ║             ")
    print("             ║                        ║             ")
    print("             ║       QUANTUMPDB       ║             ")
    print("             ║  [quantumpdb.rtfd.io]  ║             ")
    print("             ╚═══════════╗╔═══════════╝             ")
    print("                 ╔═══════╝╚═══════╗                 ")
    print("                 ║ THE KULIK LAB  ║                 ")
    print("                 ╚═══════╗╔═══════╝                 ")
    print("  ╔══════════════════════╝╚══════════════════════╗  ")
    print("  ║   Code: github.com/davidkastner/quantumpdb   ║  ")
    print("  ║   Docs: quantumpdb.readthedocs.io            ║  ")
    print("  ║      - Clusters: qp run -c config.yaml       ║  ")
    print("  ║      - QM calcs: qp submit -j                ║  ")
    print("  ╚══════════════════════════════════════════════╝  \n")

# Welcome even if no flags
welcome()

# Read in the configuration yaml file
def read_config(config_file):
    with open(config_file, 'r') as file:
        return yaml.safe_load(file)

@click.group()
def cli():
    """CLI entry point"""
    pass

@cli.command()
@click.option("--config", "-c", required=True, type=click.Path(exists=True), help="Path to the configuration YAML file")
def run(config):
    """Generates quantumPDB structures and files."""

    config_data = read_config(config)

    # Parse configuration parameters
    modeller = config_data.get('modeller', False)
    protoss = config_data.get('protoss', False)
    coordination = config_data.get('coordination', False)
    skip = config_data.get('skip', False)

    import os
    from Bio.PDB.PDBExceptions import PDBIOException
    from qp.checks import fetch_pdb

    if modeller:
        from qp.structure import missing_loops
        click.echo("MODELLER parameters:")
        optimize = config_data.get('optimize_select_residues', 1)

    if coordination:
        from qp.cluster import coordination_spheres
        limit = config_data.get('number_of_spheres', 2)
        first_sphere_radius = config_data.get('radius_of_first_sphere', 4.0)
        ligands = config_data.get('additional_ligands', [])
        capping = config_data.get('capping_method', 1)

        # Prompt user for their preferred cluster smoothing method
        smooth_choice = config_data.get('smoothing_method', 2)
        smooth_options = {0: {}, 1: {"eps": 6, "min_samples": 3}, 2: {"mean_distance": 3}, 3: {}}
        smooth_method_options = {0: "box_plot", 1: "dbscan", 2: "dummy_atom", 3: False}
        smooth_params = smooth_options[smooth_choice]
        smooth_method = smooth_method_options[smooth_choice]

        center_residues = config_data.get('center_residues', [])
        charge = config_data.get('compute_charges', True)
        count = config_data.get('count_residues', True)
        xyz = config_data.get('write_xyz', True)

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

    i = config_data.get('input', [])
    # If the input was a path or a single PDB, convert it to a list
    if isinstance(i, str):
        i = [i]
    o = config_data.get('output_dir', '')
    o = os.path.abspath(o)
    pdb_all = fetch_pdb.parse_input(i, o)

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
                import qp
                pdbl = pdb.lower()
                prepared_flag = False
                for path in qp.__path__:
                    prepared_prot_path = os.path.join(path, f"prepared/{pdbl}/Protoss")
                    if os.path.exists(prepared_prot_path):
                        prepared_flag = True
                if prepared_flag:
                    os.makedirs(prot_path, exist_ok=True)
                    from shutil import copy
                    copy(os.path.join(prepared_prot_path, f"{pdbl}_protoss.pdb"), f"{prot_path}/{pdb}_protoss.pdb")
                    copy(os.path.join(prepared_prot_path, f"{pdbl}_ligands.sdf"), f"{prot_path}/{pdb}_ligands.sdf")
                    click.echo("> Using pre-prepared protoss file")
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
                old_path = f"{prot_path}/{pdb}_protoss_old.pdb"
                if not os.path.exists(old_path):
                    from shutil import copy
                    copy(path, old_path)
                add_hydrogens.adjust_activesites(path, center_residues)

            if charge:
                ligand_charge = add_hydrogens.compute_charge(f"{prot_path}/{pdb}_ligands.sdf", path)
                ligand_spin = add_hydrogens.compute_spin(f"{prot_path}/{pdb}_ligands.sdf")
            else:
                ligand_charge = dict()

            cluster_paths = coordination_spheres.extract_clusters(
                path, f"{o}/{pdb}", center_residues,
                limit, ligands, capping, charge, count, xyz, first_sphere_radius, 
                ligand_charge=ligand_charge,
                smooth_method=smooth_method,
                **smooth_params
            )

            if charge:
                with open(f"{o}/{pdb}/charge.csv", "a") as f:
                    f.write("\n")
                    for k, v in sorted(ligand_charge.items()):
                        f.write(f"{k},{v}\n")
                with open(f"{o}/{pdb}/spin.csv", "a") as f:
                    f.write("\n")
                    for k, v in sorted(ligand_spin.items()):
                        f.write(f"{k},{v}\n")
                
                from qp.checks.charge_count import check_charge
                for cluster_path in cluster_paths:
                    check_charge(cluster_path)
                        

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
        job_manager.manage_jobs(job_count, master_list_path, minimization, basis, method, guess, gpus, memory)


    if failure_checkup:
        from qp.manager import failure_checkup
        qm_job_dir = input("What is the name of your QM job directory? ")
        failure_counts = failure_checkup.check_all_jobs(qm_job_dir)
        failure_checkup.plot_failures(failure_counts)


if __name__ == "__main__":
    # Run the command-line interface when this script is executed
    welcome()
    cli()
