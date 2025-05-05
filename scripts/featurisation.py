from rdkit import Chem
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator
import pandas as pd
import numpy as np

# Load dataset
df = pd.read_csv("data/merged_dataset.csv")

# Initialize Morgan fingerprint generator
generator = GetMorganGenerator(radius=2, fpSize=2048)

# Function to generate fingerprint from SMILES
def smiles_to_fingerprint(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    fp = generator.GetFingerprint(mol)
    return np.array(fp)

# Apply fingerprint function to each SMILES
fingerprints = df["canonical_smiles"].apply(smiles_to_fingerprint)

# Drop rows with invalid SMILES
valid_indices = fingerprints.notnull()
df = df[valid_indices].reset_index(drop=True)
fingerprints = fingerprints[valid_indices].reset_index(drop=True)

# Stack all fingerprints into a 2D array
X = np.stack(fingerprints)

# Save features and labels
np.save("data/features.npy", X)
df[["activity_label"]].to_csv("data/labels.csv", index=False)

print(f"✅ Featurization complete. Shape: {X.shape}")



