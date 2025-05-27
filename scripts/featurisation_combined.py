import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from tqdm import tqdm

# Load the saved Morgan fingerprint array
X_fp = np.load('data/features.npy')  # Replace with your actual path

# Load your CSV with SMILES
df = pd.read_csv('data/merged_dataset.csv')  # Make sure it has a 'smiles' column

# Define descriptors to compute
descriptor_names = [
    'MolWt', 'MolLogP', 'NumRotatableBonds',
    'NumHAcceptors', 'NumHDonors', 'TPSA', 'RingCount'
]

def calc_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return [np.nan] * len(descriptor_names)
    return [
        Descriptors.MolWt(mol),
        Descriptors.MolLogP(mol),
        Descriptors.NumRotatableBonds(mol),
        Descriptors.NumHAcceptors(mol),
        Descriptors.NumHDonors(mol),
        Descriptors.TPSA(mol),
        Descriptors.RingCount(mol)
    ]

# Generate descriptor features
descriptor_features = [calc_descriptors(smi) for smi in tqdm(df['canonical_smiles'], desc="Calculating descriptors")]
X_desc = np.array(descriptor_features)

# Combine Morgan fingerprints with descriptors
X_combined = np.concatenate([X_fp, X_desc], axis=1)

# Generate feature names
fp_feature_names = [f'bit_{i}' for i in range(X_fp.shape[1])]
combined_feature_names = fp_feature_names + descriptor_names

# Save combined features with headers
df_combined = pd.DataFrame(X_combined, columns=combined_feature_names)
df_combined.to_csv('data/features_combined.csv', index=False)

# Save the combined features
np.save('data/features_combined.npy', X_combined)

print("Combined feature shape:", X_combined.shape)

