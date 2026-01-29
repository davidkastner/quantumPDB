"""Command-line interface (CLI) entry point.

**Usage** ``qp run --config path/to/config.yaml``
"""

import os
import sys
import click
import traceback

def welcome():
    """Print the QuantumPDB welcome banner with ASCII art and usage hints.

    Displays the package logo, links to documentation and GitHub, and
    basic usage examples. Called automatically when the CLI is invoked.
    """
    click.secho("\n")
    click.secho("             ╔════════════════════════╗             ", bold=True)
    click.secho("             ║       __________       ║             ", bold=True)
    click.secho("             ║     / ____/\____ \     ║             ", bold=True)
    click.secho("             ║    < <_|  ||  |_> >    ║             ", bold=True)
    click.secho("             ║     \__   ||   __/     ║             ", bold=True)
    click.secho("             ║        |__||__|        ║             ", bold=True)
    click.secho("             ║                        ║             ", bold=True)
    click.secho("             ║       QUANTUMPDB       ║             ", bold=True)
    click.secho("             ║  [quantumpdb.rtfd.io]  ║             ", bold=True)
    click.secho("             ╚═══════════╗╔═══════════╝             ", bold=True)
    click.secho("                 ╔═══════╝╚═══════╗                 ", bold=True)
    click.secho("                 ║ THE KULIK LAB  ║                 ", bold=True)
    click.secho("                 ╚═══════╗╔═══════╝                 ", bold=True)
    click.secho("  ╔══════════════════════╝╚══════════════════════╗  ", bold=True)
    click.secho("  ║   Code: github.com/davidkastner/quantumpdb   ║  ", bold=True)
    click.secho("  ║   Docs: quantumpdb.readthedocs.io            ║  ", bold=True)
    click.secho("  ║      - Clusters: qp run -c config.yaml       ║  ", bold=True)
    click.secho("  ║      - QM calcs: qp submit -c config.yaml    ║  ", bold=True)
    click.secho("  ╚══════════════════════════════════════════════╝\n", bold=True)

# Welcome even if no flags
welcome()

@click.group()
def cli():
    """QuantumPDB command-line interface.

    Use ``qp run``, ``qp submit``, or ``qp analyze`` with a configuration
    YAML file to execute different stages of the pipeline.
    """
    pass

@cli.command()
@click.option("--config", "-c", required=True, type=click.Path(exists=True), help="Path to the configuration YAML file")
def run(config):
    """Execute the structure preparation, protonation, and cluster extraction pipeline.

    This command runs stages 1-3 of the QuantumPDB workflow: fetching PDB
    structures, modeling missing residues with MODELLER, assigning protonation
    states with Protoss, and extracting QM cluster models using Voronoi
    tessellation.
    """
    import os
    from shutil import copy
    from Bio.PDB.PDBExceptions import PDBIOException
    from qp.structure import setup

    config_data = setup.read_config(config)
    err = {"PDB": [], "Protoss": [], "Coordination sphere": [], "Other": []}

    # Parse configuration parameters
    modeller = config_data.get('modeller', False)
    protoss = config_data.get('protoss', False)
    coordination = config_data.get('coordination', False)
    skip = config_data.get('skip', 'all')
    max_clash_refinement_iter = config_data.get('max_clash_refinement_iter', 5)
    
    input = config_data.get('input', [])
    output = config_data.get('output_dir', '')
    center_yaml_residues = config_data.get('center_residues', [])
    pdb_all, center_residues = setup.parse_input(input, output, center_yaml_residues)

    if modeller:
        from qp.structure import missing
        optimize = config_data.get('optimize_select_residues', 1)
        convert_to_nhie_oxo = config_data.get('convert_to_nhie_oxo', False)
    if protoss:
        from qp.protonate import get_protoss, fix
    if coordination:
        from qp.cluster import spheres
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
        smooth_choice = config_data.get('smoothing_method', 2)
        smooth_options = {0: {}, 1: {"eps": 6, "min_samples": 3}, 2: {"mean_distance": 3}, 3: {}}
        smooth_method_options = {0: "box_plot", 1: "dbscan", 2: "dummy_atom", 3: False}
        smooth_params = smooth_options[smooth_choice]
        smooth_method = smooth_method_options[smooth_choice]

        if capping or charge:
            protoss = True

    for pdb, path in pdb_all:
        try:
            click.secho("╔══════╗", bold=True)
            click.secho(f"║ {pdb.upper()} ║", bold=True)
            click.secho("╚══════╝", bold=True)

            # Skips fetching if PDB file exists
            if not os.path.isfile(path):
                click.echo(f"> Fetching PDB file")
                try:
                    setup.fetch_pdb(pdb, path)
                except ValueError as e:
                    # Catches invalid PDB ID
                    click.secho(f"> Error: {e}\n", italic=True, fg="red")
                    err["PDB"].append(pdb)
                    continue
                except IOError as e:
                    # Catches network and server issues
                    click.secho(f"> Error: A server or network issue occurred. {e}\n", italic=True, fg="red")
                    err["PDB"].append(pdb)
                    continue
            
            # Extract the current center residue from the list of all residues
            from qp.cluster.spheres import CenterResidue
            center_residue = CenterResidue(center_residues.pop(0))
            click.echo(f"> Using center residue: {center_residue}")
            
            residues_with_clashes = [] # Start by assuming no protoss clashes
            for i in range(max_clash_refinement_iter):
                mod_path = f"{output}/{pdb}/{pdb}_modeller.pdb"
                if modeller:
                    if skip in ["modeller", "all"] and os.path.isfile(mod_path) and not residues_with_clashes:
                        click.echo("> MODELLER file found")
                    else:
                        click.echo("> Building model")
                        AA = missing.define_residues()
                        residues = missing.get_residues(path, AA)

                        # Remove trailing missing residues from the ends of all chains
                        residues = missing.clean_termini(residues)

                        # Update residues with clashes if any
                        if residues_with_clashes:
                            missing.delete_clashes(path, residues_with_clashes)
                            for residue_with_clashes in residues_with_clashes:
                                res_id = residue_with_clashes[0]
                                one_letter_code = residue_with_clashes[1]
                                chain_index = residue_with_clashes[2]
                                for j, residue in enumerate(residues[chain_index]):
                                        if residue[0][0] == res_id:
                                            if residue[1] == one_letter_code:
                                                residues[chain_index][j] = ((res_id, ' '), one_letter_code, 'R')
                                            else:
                                                print("> Residue index matched but residue name did not.")
                                                exit()

                        ali_path = f"{output}/{pdb}/{pdb}.ali"
                        missing.write_alignment(residues, pdb, path, ali_path)
                        print("> Generated alignment file and starting Modeller run:\n")
                        if not os.path.isfile(f"{output}/{pdb}/{pdb}.pdb"):
                            copy(path, f"{output}/{pdb}/{pdb}.pdb")
                        missing.build_model(residues, pdb, path, ali_path, mod_path, optimize)

                prot_path = f"{output}/{pdb}/Protoss"
                protoss_log_file = f"{prot_path}/{pdb}_log.txt"
                protoss_pdb = f"{prot_path}/{pdb}_protoss.pdb"
                ligands_sdf = f"{prot_path}/{pdb}_ligands.sdf"

                if protoss:
                    if skip in ["protoss", "all"] and os.path.isfile(protoss_pdb) and not residues_with_clashes:
                        click.echo("> Protoss file found")
                    else:
                        import qp
                        pdbl = pdb.lower()
                        prepared_flag = False
                        for qp_path in qp.__path__:
                            prepared_prot_path = os.path.join(qp_path, f"resources/prepared/{pdbl}/Protoss")
                            if os.path.exists(prepared_prot_path):
                                prepared_flag = True
                        if prepared_flag:
                            os.makedirs(prot_path, exist_ok=True)
                            copy(os.path.join(prepared_prot_path, f"{pdbl}_protoss.pdb"), protoss_pdb)
                            copy(os.path.join(prepared_prot_path, f"{pdbl}_ligands.sdf"), ligands_sdf)
                            click.echo("> Using pre-prepared protoss file")
                        else:
                            from qp.protonate.fix import clean_occupancy, fix_OXT
                            click.echo("> Running Protoss")
                            if modeller:
                                pdb_path = mod_path
                            else:
                                pdb_path = path
                              
                            fix_OXT(pdb_path)
                            clean_occupancy(pdb_path, center_residue)

                            try:
                                pid = get_protoss.upload(pdb_path)
                            except ValueError:
                                click.secho("> Error: Could not upload PDB file\n", italic=True, fg="red")
                                # Occurs when uploading a PDB file > 4MB
                                # TODO retry with parsed PDB code? Will raise Bio.PDB error if
                                #      number of atoms in output exceeds 99999
                                err["Protoss"].append(pdb)
                                continue
                            
                            job = get_protoss.submit(pid)
                            get_protoss.download(job, protoss_pdb, "protein")
                            get_protoss.download(job, ligands_sdf, "ligands")
                            get_protoss.download(job, protoss_log_file, "log")
                            get_protoss.repair_ligands(protoss_pdb, pdb_path)

                # Get any residues identified by Protoss as problematic
                from qp.protonate.parse_output import parse_log
                AA = missing.define_residues()
                residues_with_clashes = parse_log(protoss_log_file, protoss_pdb, AA)
                if residues_with_clashes:
                    print(f"> WARNING: Protoss removed {len(residues_with_clashes)} residues due to clashes")
                    clash_details = ", ".join([f"{res_name}{res_id} in chain {chain}" for res_id, _, _, chain, res_name in sorted(list(residues_with_clashes))])
                    print(f"> WARNING: Remodeling {clash_details} with Modeller.")
                else:
                    break # Don't loop again as no clashes were detected

            if protoss and convert_to_nhie_oxo:
                click.echo("> Converting AKG to reactive OXO and SIN state")
                from qp.structure.convert_nhie_oxo import add_oxo_sin
                add_oxo_sin(protoss_pdb)
                click.echo("> Removing hydrogens from OXO's")
                
                from qp.structure.convert_nhie_oxo import remove_oxo_hydrogens, update_oxo_sdf, update_sin_sdf
                remove_oxo_hydrogens(protoss_pdb)
                
                click.echo("> Adding OXO to ligands.sdf")
                update_oxo_sdf(protoss_pdb, ligands_sdf)
                
                click.echo("> Updating SIN in ligands.sdf")
                sin_ligands_sdf = f"{prot_path}/{pdb}_ligands_sin.sdf"
                click.echo("> Uploading new succinate structure to Protoss")
                pid = get_protoss.upload(protoss_pdb)
                job = get_protoss.submit(pid) # Send Protoss the new SIN ligands
                get_protoss.download(job, sin_ligands_sdf, "ligands") # Get the .sdf file for the SIN ligands
                update_sin_sdf(ligands_sdf, sin_ligands_sdf)

            # Set the path to the current structure, there may be a better way to do this
            if modeller:
                path = mod_path
            if protoss:
                # Check if any atoms changed from HETATM to ATOM or vice versa for future troubleshooting purposes
                from qp.protonate.parse_output import record_type
                changed_residues = record_type(mod_path, protoss_pdb) # I want to do something with this eventually

                path = f"{prot_path}/{pdb}_protoss.pdb"
                old_path = f"{prot_path}/{pdb}_protoss_orig.pdb"
                if not os.path.exists(old_path):
                    from shutil import copy
                    copy(path, old_path)
                fix.adjust_activesites(path, center_residue)

            if coordination:
                from qp.protonate.ligand_prop import compute_charge, compute_spin
                click.echo("> Extracting clusters")
                if charge:
                    ligand_charge = compute_charge(f"{prot_path}/{pdb}_ligands.sdf", path)
                    ligand_spin = compute_spin(f"{prot_path}/{pdb}_ligands.sdf")
                else:
                    ligand_charge = dict()
                cluster_paths = spheres.extract_clusters(
                    path, f"{output}/{pdb}", center_residue, sphere_count, 
                    first_sphere_radius, max_atom_count, merge_cutoff, smooth_method,
                    ligands, capping, charge, ligand_charge, count, xyz, hetero_pdb, include_ligands,
                    **smooth_params
                )

                if charge:
                    charge_csv_path = f"{output}/{pdb}/charge.csv"
                    with open(charge_csv_path, "a") as f:
                        f.write("\n")
                        for k, v in sorted(ligand_charge.items()):
                            f.write(f"{k},{v}\n")
                    if ligand_spin:
                        with open(f"{output}/{pdb}/spin.csv", "w") as f:
                            f.write("\n")
                            for k, v in sorted(ligand_spin.items()):
                                f.write(f"{k},{v}\n")
                    
                    from qp.cluster.charge_count import check_charge
                    for cluster_path in cluster_paths:
                        check_charge(cluster_path)               
            click.echo("")

            for k, v in err.items():
                if v:
                    click.echo(click.style(k + " errors: ", bold=True, fg="red") + ", ".join(v))
    
        except:
            # Log the exception details to stderr, which is already redirected to log.out
            click.echo(f"> CRITICAL FAILURE: Error processing {pdb.upper()}", err=True)
            traceback.print_exc(file=sys.stderr)


@cli.command()
@click.option("--config", "-c", required=True, type=click.Path(exists=True), help="Path to the configuration YAML file")
def submit(config):
    """Create QM input files and submit calculations to the job scheduler.

    This command runs stage 4 of the QuantumPDB workflow: generating
    TeraChem or ORCA input files and submitting them to SLURM or SGE
    job schedulers.
    """

    from qp.structure import setup
    from qp.manager import create
    from qp.manager import submit
    
    # Parse configuration parameters
    config_data = setup.read_config(config)
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
    create_jobs = config_data.get('create_jobs', False)
    submit_jobs = config_data.get('submit_jobs', False)
    input = config_data.get('input', [])
    output = config_data.get('output_dir', '')
    
    if not os.path.exists(input):
        raise FileNotFoundError(f"Could not find input file named {input}.")
    input = os.path.abspath(input)

    if create_jobs:
        click.echo("> Creating job files for QM calculations")
        create.create_jobs(input, output, optimization, basis, method, guess, charge_embedding, charge_embedding_cutoff, gpus, memory, scheduler, pcm_radii_file, dielectric)
    if submit_jobs:
        click.echo("\n> Submitting QM calculations")
        submit.manage_jobs(output, job_count, method, scheduler)


@cli.command()
@click.option("--config", "-c", required=True, type=click.Path(exists=True), help="Path to the configuration YAML file")
def analyze(config):
    """Check job status and run post-processing analysis with Multiwfn.

    This command runs stage 5 of the QuantumPDB workflow: classifying
    job outcomes (done, failed, queued), and optionally computing
    partial charges and dipole moments using Multiwfn.
    """

    from qp.structure import setup

    config_data = setup.read_config(config)
    method = config_data.get('method', 'wpbeh')
    job_checkup = config_data.get('job_checkup', False)
    delete_queued = config_data.get('delete_queued', False)
    calc_charge_schemes = config_data.get('calc_charge_schemes', False)
    calc_dipole = config_data.get('calc_dipole', False)
    multiwfn_path = config_data.get('multiwfn_path', 'Multiwfn')
    charge_scheme = config_data.get('charge_scheme', 'Hirshfeld')
    input = config_data.get('input', [])
    output = config_data.get('output_dir', '')
    center_yaml_residues = config_data.get('center_residues', [])
    pdb_all, center_residues = setup.parse_input(input, output, center_yaml_residues)
    
    if job_checkup:
        from qp.analyze import checkup
        click.echo("> Checking to see the status of the jobs...")
        failure_counts, author_counts = checkup.check_all_jobs(method, output, delete_queued)
        checkup.plot_failures(failure_counts)
        checkup.plot_authors(author_counts)
        checkup.write_author_credit_csv(author_counts)
        checkup.plot_failure_modes_from_csv(os.path.join("checkup", "failure_modes.csv"))

    if calc_charge_schemes or calc_dipole:
        from qp.analyze import multiwfn
        click.echo("> Activate Multiwfn module (e.g., module load multiwfn/noGUI_3.7)")
        settings_ini_path = multiwfn.get_settings_ini_path()
        atmrad_path = multiwfn.get_atmrad_path()
        if calc_charge_schemes:
            click.echo("> Calculating charge schemes...")
            multiwfn.iterate_qm_output(pdb_all, method, output, multiwfn_path, settings_ini_path, atmrad_path, charge_scheme, multiwfn.charge_scheme)
        if calc_dipole:
            click.echo("> Calculating dipoles using center of mass as the reference...")
            multiwfn.iterate_qm_output(pdb_all, method, output, multiwfn_path, settings_ini_path, atmrad_path, charge_scheme, multiwfn.calc_dipole)


if __name__ == "__main__":
    # Run the command-line interface when this script is executed
    welcome()
    cli()
