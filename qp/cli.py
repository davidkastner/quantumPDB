"""Command-line interface (CLI) entry point.

**Usage** ``qp run --config path/to/config.yaml``
"""

import os
import yaml
import click

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
    print("  ║      - QM calcs: qp submit -c config.yaml    ║  ")
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
    skip = config_data.get('skip', 'all')

    import os
    from Bio.PDB.PDBExceptions import PDBIOException
    from qp.checks import fetch_pdb

    if modeller:
        from qp.structure import missing_loops
        optimize = config_data.get('optimize_select_residues', 1)
        convert_to_oxo = config_data.get('convert_to_oxo', False)

    if coordination:
        from qp.cluster import coordination_spheres
        center_residues = config_data.get('center_residues', [])
        sphere_count = config_data.get('number_of_spheres', 2)
        first_sphere_radius = config_data.get('radius_of_first_sphere', 4.0)
        max_atom_count = config_data.get('max_atom_count', None)
        merge_cutoff = config_data.get('merge_distance_cutoff', 0.0)

        include_ligands = config_data.get('include_ligands', 2)
        ligands = config_data.get('additional_ligands', [])
        capping = config_data.get('capping_method', 1)
        charge = config_data.get('compute_charges', True)
        count = config_data.get('count_residues', True)
        xyz = config_data.get('write_xyz', True)
        hetero_pdb = config_data.get('write_hetero_pdb', False)

        # Prompt user for their preferred cluster smoothing method
        smooth_choice = config_data.get('smoothing_method', 2)
        smooth_options = {0: {}, 1: {"eps": 6, "min_samples": 3}, 2: {"mean_distance": 3}, 3: {}}
        smooth_method_options = {0: "box_plot", 1: "dbscan", 2: "dummy_atom", 3: False}
        smooth_params = smooth_options[smooth_choice]
        smooth_method = smooth_method_options[smooth_choice]

        if capping or charge:
            protoss = True

    if protoss:
        from qp.structure import add_hydrogens

    err = {
        "PDB": [],
        "Protoss": [],
        "Coordination sphere": [],
        "Other": []
    }

    input = config_data.get('input', [])
    # If the input was a path or a single PDB, convert it to a list
    if isinstance(input, str):
        input = [input]
    output = config_data.get('output_dir', '')
    output = os.path.abspath(output)
    pdb_all = fetch_pdb.parse_input(input, output)

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

        mod_path = f"{output}/{pdb}/{pdb}_modeller.pdb"
        if modeller:
            if skip in ["modeller", "all"] and os.path.isfile(mod_path):
                click.echo("> MODELLER file found")
            else:
                click.echo("> Building model")
                AA = missing_loops.define_residues()
                residues = missing_loops.get_residues(path, AA)
                # Remove trailing missing residues from the ends of all chains
                residues = missing_loops.clean_termini(residues)

                ali_path = f"{output}/{pdb}/{pdb}.ali"
                missing_loops.write_alignment(residues, pdb, path, ali_path)
                missing_loops.build_model(residues, pdb, path, ali_path, mod_path, optimize)

            if convert_to_oxo:
                click.echo("> Requested to convert AKG to reactive OXO and SUC")
                from qp.structure.convert_to_oxo import add_oxo_and_suc
                add_oxo_and_suc(mod_path)

        prot_path = f"{output}/{pdb}/Protoss"
        if protoss:
            if skip in ["protoss", "all"] and os.path.isfile(f"{prot_path}/{pdb}_protoss.pdb"):
                click.echo("> Protoss file found")
            else:
                import qp
                pdbl = pdb.lower()
                prepared_flag = False
                for qp_path in qp.__path__:
                    prepared_prot_path = os.path.join(qp_path, f"prepared/{pdbl}/Protoss")
                    if os.path.exists(prepared_prot_path):
                        prepared_flag = True
                if prepared_flag:
                    os.makedirs(prot_path, exist_ok=True)
                    from shutil import copy
                    copy(os.path.join(prepared_prot_path, f"{pdbl}_protoss.pdb"), f"{prot_path}/{pdb}_protoss.pdb")
                    copy(os.path.join(prepared_prot_path, f"{pdbl}_ligands.sdf"), f"{prot_path}/{pdb}_ligands.sdf")
                    click.echo("> Using pre-prepared protoss file")
                else:
                    from qp.structure.add_hydrogens import clean_partial_occupancy
                    click.echo("> Running Protoss")
                    if modeller:
                        path = mod_path
                    clean_partial_occupancy(path, center_residues)

                    try:
                        pid = add_hydrogens.upload(path)
                    except ValueError:
                        click.secho("Error uploading PDB file\n", italic=True, fg="red")
                        # Occurs when uploading a PDB file > 4MB
                        # TODO retry with parsed PDB code? Will raise Bio.PDB error if
                        #      number of atoms in output exceeds 99999
                        err["Protoss"].append(pdb)
                        continue
                    
                    protoss_pdb = f"{prot_path}/{pdb}_protoss.pdb"
                    job = add_hydrogens.submit(pid)
                    add_hydrogens.download(job, protoss_pdb, "protein")
                    add_hydrogens.download(job, f"{prot_path}/{pdb}_ligands.sdf", "ligands")
                    add_hydrogens.download(job, f"{prot_path}/{pdb}_log.txt", "log")
                    add_hydrogens.repair_ligands(protoss_pdb, path)

            if convert_to_oxo:
                from qp.structure.convert_to_oxo import remove_oxo_hydrogens
                remove_oxo_hydrogens(protoss_pdb)


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
                path, f"{output}/{pdb}", center_residues, sphere_count, 
                first_sphere_radius, max_atom_count, merge_cutoff, smooth_method,
                ligands, capping, charge, ligand_charge, count, xyz, hetero_pdb, include_ligands,
                **smooth_params
            )

            if charge:
                with open(f"{output}/{pdb}/charge.csv", "a") as f:
                    f.write("\n")
                    for k, v in sorted(ligand_charge.items()):
                        f.write(f"{k},{v}\n")
                if ligand_spin:
                    with open(f"{output}/{pdb}/spin.csv", "w") as f:
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
@click.option("--config", "-c", required=True, type=click.Path(exists=True), help="Path to the configuration YAML file")
@click.option("--failure_checkup", "-f", is_flag=True, help="Find failed structures")
def submit(config, failure_checkup):
    """Handles the submission of jobs for the quantumPDB."""

    from qp.job_manager import create
    from qp.job_manager import submit
    
    if config:
        config_data = read_config(config)
        optimization = config_data.get('optimization', False)
        method = config_data.get('method', 'wpbeh')
        basis = config_data.get('basis', 'lacvps_ecp')
        guess = config_data.get('guess', 'generate')
        gpus = config_data.get('gpus', 1)
        memory = config_data.get('memory', '8G')
        scheduler = config_data.get('scheduler', 'slurm')
        pcm_radii_file = config_data.get('pcm_radii_file', 'pcm_radii')
        job_count = config_data.get('job_count', 80)
        charge_embedding = config_data.get('charge_embedding', False)
        charge_embedding_cutoff = config_data.get('charge_embedding_cutoff', 20)
        dielectric = config_data.get('dielectric', 10)
        submit_jobs = config_data.get('submit_jobs', False)

        # Check if a config file and end if it was a PDB
        output = config_data.get('output_dir', '') # Ensure execution from the correct directory
        input = config_data.get('input', [])
        
        if not os.path.exists(input):
            raise FileNotFoundError(f"Could not find input file named {input}.")
        input = os.path.abspath(input)

        click.echo("   > Creating job files for QM calculations")
        create.create_jobs(input, output, optimization, basis, method, guess, charge_embedding, charge_embedding_cutoff, gpus, memory, scheduler, pcm_radii_file, dielectric)
        if submit_jobs:
            click.echo("\n   > Submitting QM calculations")
            submit.manage_jobs(output, job_count, method, scheduler)


    if failure_checkup:
        from qp.job_manager import failure_checkup
        qm_job_dir = input("What is the name of your QM job directory? ")
        failure_counts = failure_checkup.check_all_jobs(qm_job_dir)
        failure_checkup.plot_failures(failure_counts)


if __name__ == "__main__":
    # Run the command-line interface when this script is executed
    welcome()
    cli()
