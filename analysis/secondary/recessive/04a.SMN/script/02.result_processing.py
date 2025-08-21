#!/usr/bin/env python3
import pandas as pd

# Load data
smn_output = "/home/_shared/jscliu/project/2025/Flagship/analysis/secondary/recessive/04a.SMN/summary/smn_output.txt"
filtered_output = "/home/_shared/jscliu/project/2025/Flagship/analysis/secondary/recessive/04a.SMN/summary/filtered_smn_output.tsv"
df = pd.read_csv(smn_output, sep='\t')

# Filter for carrier + silent carrier
filtered_df = df[
    (df['isCarrier'] == True) |
    (
        (df['isCarrier'] == False) &
        (df['SMN1_CN'] == 2) &
        (df['g.27134T>G_CN'] == 1)
    )
].copy()

filtered_df.to_csv(filtered_output, index=False)
print(filtered_df)
