import pandas as pd

# Load dataset
df = pd.read_csv("data/antimalarial_activities_selected.csv")

# Convert standard_value to numeric, coercing errors
df['standard_value'] = pd.to_numeric(df['standard_value'], errors='coerce')

# Drop ligand_efficiency column if it exists
if 'ligand_efficiency' in df.columns:
    df = df.drop(columns=['ligand_efficiency'])

# Filter for standard_types of interest
valid_types = ['IC50', 'EC50', 'ED50']
df = df[df['standard_type'].isin(valid_types)]

# Drop rows with missing data in important fields
df = df.dropna(subset=['standard_value', 'standard_units', 'canonical_smiles'])

# Convert all values to nM
def convert_to_nM(row):
    if row['standard_units'].lower() == 'um':
        return row['standard_value'] * 1000  # µM to nM
    elif row['standard_units'].lower() == 'nm':
        return row['standard_value']
    else:
        return None  # discard unknown units

df['standard_value_nM'] = df.apply(convert_to_nM, axis=1)

# Drop rows where conversion failed
df = df.dropna(subset=['standard_value_nM'])

# Assign activity labels
df['activity_label'] = df['standard_value_nM'].apply(lambda x: 1 if x <= 1000 else 0)

# Keep only necessary columns
columns_to_keep = ['molecule_chembl_id', 'canonical_smiles', 'standard_type',
                   'standard_value_nM', 'activity_label']
df = df[columns_to_keep]

# Save cleaned and labeled dataset
df.to_csv("data/antimalarial_cleaned_labeled.csv", index=False)

print("✅ Cleaned and labeled data saved to 'data/antimalarial_cleaned.csv'")
print("\n📊 Activity Label Distribution:")
print(df['activity_label'].value_counts())