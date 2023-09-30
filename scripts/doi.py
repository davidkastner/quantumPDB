import pandas as pd
import requests

def get_doi_url(pdb_code):
    """Fetch the DOI for a given PDB code using the PDB REST API."""
    url = f'https://data.rcsb.org/rest/v1/core/entry/{pdb_code}'
    response = requests.get(url)
    
    if response.status_code != 200:
        print(f"Error: Failed to fetch data for PDB: {pdb_code}. Status code: {response.status_code}")
        return None

    response_data = response.json()

    # Extract the DOI from the response data
    doi = response_data.get('pdbx_database_id_doi')
    if doi:
        return doi
    else:
        print(f"No DOI found for PDB: {pdb_code}")
        return None

def main():
    # Read the CSV file
    df = pd.read_csv('protein_master_list_charge_multiplicity.csv')

    # Total number of PDBs
    total_pdbs = len(df)

    # Counter for API calls since last save
    api_calls_since_last_save = 0

    # For each PDB code, if there isn't already a DOI URL, fetch the DOI and set the DOI URL in Column 6
    for index, row in df.iterrows():
        pdb_code = row.iloc[0]
        existing_doi_url = row.iloc[5]
        
        # Check if the DOI URL is already present
        if pd.isnull(existing_doi_url):
            doi = get_doi_url(pdb_code)
            if doi:
                doi_url = f'https://doi.org/{doi}'
                df.at[index, df.columns[5]] = doi_url

            # Increment the API call counter
            api_calls_since_last_save += 1

            # Print the progress update
            print(f"> Finished PDB: {pdb_code} ({index + 1}/{total_pdbs})")

            # If 5 API calls have been made since the last save, save the dataframe and reset the counter
            if api_calls_since_last_save == 5:
                df.to_csv('doi_charge_multiplicity.csv', index=False)
                api_calls_since_last_save = 0

    # Save one final time after the loop is done to catch any remaining changes
    df.to_csv('doi_charge_multiplicity.csv', index=False)

if __name__ == "__main__":
    main()
