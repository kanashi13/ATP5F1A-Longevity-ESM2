from Bio import SeqIO
import pandas as pd

fasta_file = "FASTA/ATP5F1A_combined_sequence.txt"
records = list(SeqIO.parse(fasta_file, "fasta"))

sequences = {rec.description.split("[")[-1].split("]")[0].strip(): str(rec.seq) for rec in records}

df = pd.DataFrame.from_dict(sequences, orient='index')
df = df.transpose()

variable_sites = []
for pos in range(len(df)):
    unique_aa = df.iloc[pos].unique()
    if len(unique_aa) > 1:
        variable_sites.append({
            'Position': pos+1,
            'Amino Acids': {sp: aa for sp, aa in df.iloc[pos].items()},
            'Num Variants': len(unique_aa)
        })

var_df = pd.DataFrame(variable_sites)
print(f"Number of variable positions: {len(var_df)} out of 553")
print(var_df)

var_df.to_excel("variable_positions.xlsx", index=False)
print("File variable_positions.xlsx created")