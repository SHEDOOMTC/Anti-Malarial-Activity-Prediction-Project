import pandas as pd
import requests
from io import StringIO

def download_pubchem_csv(aid):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/{aid}/CSV"
    response = requests.get(url)
    if response.status_code == 200:
        df = pd.read_csv(StringIO(response.text))
        df['AID'] = aid
        return df
    else:
        print(f"Failed to fetch AID {aid}")
        return None

# Download and merge
aids = [449703, 1619, 1828,  1815,  489011,  489016, 524796, 504696, 2195, 2196, 1822, 606570, 1159566]
dfs = [download_pubchem_csv(aid) for aid in aids if download_pubchem_csv(aid) is not None]
merged = pd.concat(dfs, ignore_index=True)

# Save merged dataset
merged.to_csv("merged_assays.csv", index=False)
print("Merged dataset saved as 'merged_assays.csv'")
