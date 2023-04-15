"""Command-line interface (CLI) entry point."""

# Print first to welcome the user while it waits to load the modules
print("\n.-------------------------------.")
print("| WELCOME TO THE QUANTUMPDB CLI |")
print(".-------------------------------.")
print("Default programmed actions for the quantumPDB package.")
print("GitHub: https://github.com/davidkastner/quantumpdb")
print("Documenation: https://quantumpdb.readthedocs.io\n")

import click

@click.command()
@click.option("--first_task", "-a", is_flag=True, help="First task.")
def cli(
    first_task,
    ):
    """
    The overall command-line interface (CLI) entry point.
    The CLI interacts with the rest of the package.

    A complete reference of quantumPDB functionality.
    This is advantagous because it quickly introduces so molecuLearn.
    Specificaly, to the complete scope of available functionality.
    It also improves long-term maintainability and readability.

    """
    
    if first_task
        click.echo("> First processing task:")
        click.echo("> Loading...")
        import qp.process
        

    else:
        click.echo("No functionality was requested.\nTry --help.")


if __name__ == "__main__":
    # Run the command-line interface when this script is executed
    cli()
