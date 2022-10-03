import pandas as pd
from Bio import SeqIO

# read in the header of the mito scaffold
blast_output = pd.read_csv(input[1], sep="\t", header = None)
for index, value in blast_output.iterrows():
    if index == 0:
        mito_tig = value[1
# code to tag the mito contig header
to_add = "mitochondrial_"
with open(output[0], "w") as outputs:
    for r in SeqIO.parse(input[0], "fasta"):
        if r.id = str(mito_tig)
        r.id = (to_add + r.description).replace(" ", "_")
        r.description = r.id
        SeqIO.write(r, output[0], "fasta"
# code to remove the mito contig from the assembly
with open(out_folder + assembly.replace(".fasta", ".no_mito.fasta")) as mito_purge: # etc
    for seq in Seq.io(input[0]):
        if seq.header != mito_tig:
        # write to new file
        SeqIO.write(r, output[1], "fasta")