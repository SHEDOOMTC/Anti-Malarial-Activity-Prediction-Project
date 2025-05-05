import pandas as pd

# Load the cleaned datasets
chembl_df = pd.read_csv("data/cleaned_chembl_data.csv")
pubchem_df = pd.read_csv("data/cleaned_pubchem_data.csv")

# Rename identifier columns for uniformity
chembl_df = chembl_df.rename(columns={"molecule_chembl_id": "compound_id"})
pubchem_df = pubchem_df.rename(columns={"PUBCHEM_CID": "compound_id"})
pubchem_df = pubchem_df.rename(columns={"activity_class": "activity_label"})
pubchem_df = pubchem_df.rename(columns={"Standard Type": "standard_type"})

# Concatenate
merged_df = pd.concat([chembl_df, pubchem_df], ignore_index=True)

# Drop missing SMILES or standard_value_nM
merged_df.dropna(subset=["canonical_smiles", "standard_value_nM"], inplace=True)

# Ensure standard_value_nM is numeric
merged_df["standard_value_nM"] = pd.to_numeric(merged_df["standard_value_nM"], errors="coerce")


# Retain only the row with the lowest standard_value_nM for each canonical_smiles
merged_df = merged_df.loc[merged_df.groupby("canonical_smiles")["standard_value_nM"].idxmin()].reset_index(drop=True)

# Final shape and class distribution
print(f"Merged dataset shape: {merged_df.shape}")
print("Class distribution:")
print(merged_df["activity_label"].value_counts())

# Save
merged_df.to_csv("data/merged_dataset.csv", index=False)
