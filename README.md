# Biopython_Project
In-silico Analysis and Functional Prediction of an Unknown DNA-Encoded Protein from Escherichia col

Original file is located at
    https://colab.research.google.com/drive/1kIXfwEsFQyau6kobXhuXz6lLKQ3-_1xa

##Install Biopython
"""

pip install biopython

"""##Import the necesarry packages"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction

"""##Read the Fasta file"""

record = SeqIO.read("sequence.fasta","fasta")
dna_seq = record.seq.upper()
print(f"Sequence_ID : {record.id}")
print(f"Sequence_length : {len(dna_seq)}bp")

"""##Quality analysis of Unknown Sequence

GC Content
"""

gc_cont = gc_fraction(dna_seq)*100
print(f"GC_Content : {gc_cont:.2f}%")

"""Checking DNA validation"""

valid_bases = set("ATGCN")
invalid_bases = set(dna_seq) - valid_bases
print(f"Total invalid bases : {invalid_bases}")

"""Scientific interpretation:

“The DNA sequence is very large (~4.6 Mb), GC-balanced, and contains no invalid nucleotide characters, indicating a high-quality genomic DNA sequence suitable for downstream gene prediction analysis.”

##ORF
"""

from Bio.Seq import Seq

seq = dna_seq
cleaned_seq_str = "".join([base if base in "ATGCN" else "N" for base in str(seq)])
cleaned_seq = Seq(cleaned_seq_str)

translated_frames = []

for strand, nuc in [("+", cleaned_seq), ("-", cleaned_seq.reverse_complement())]:
    for frame in range(3):
        protein = nuc[frame:].translate()
        translated_frames.append((strand, frame, protein))

"""##selecting the ORF with minimum aa length 300"""

orfs = []
min_aa_length = 300

for strand, frame, protein in translated_frames:
    for orf in protein.split("*"):
        if len(orf) >= min_aa_length:
            orfs.append((strand, frame, len(orf), orf))

print("Number of ORFs after filtering:", len(orfs))

orfs_sorted = sorted(orfs, key=lambda x: x[2], reverse=True)

for i, orf in enumerate(orfs_sorted[:5], start=1):
    print(f"ORF {i}: Length = {orf[2]} aa | Strand = {orf[0]} | Frame = {orf[1]}")

"""##Filtering the Longest ORF"""

from Bio.Seq import Seq


orfs = []
min_aa_length = 300


for strand, frame, protein in translated_frames:
    for orf in protein.split("*"):
        if len(orf) >= min_aa_length:
            orfs.append((strand, frame, len(orf), orf))


if orfs:
    longest_orf = max(orfs, key=lambda x: x[2])

    strand, frame, aa_length, protein_seq = longest_orf

    print("Longest ORF details:")
    print("Protein length (aa):".ljust(20), aa_length)
    print("Strand:".ljust(20), strand)
    print("Frame:".ljust(20), frame)
else:
    print(f"No ORFs found with a minimum amino acid length of {min_aa_length} in any frame.")
    protein_seq = Seq("")

print(protein_seq)

"""“Since the input sequence represents a complete bacterial genome, sequence filtering was performed using ORF prediction. The genome was translated in all six reading frames using Biopython, and ORFs shorter than 300 amino acids were excluded. Among the remaining candidates, the longest ORF (2387 amino acids) encoded on the forward strand was selected for downstream analysis.”

##Converting the longest orf seq into a fastafile
"""

from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

orf_record = SeqRecord(
    protein_seq,
    id="Longest_ORF_2387aa",
    description="Predicted protein from ORF filtering (Strand +, Frame 0)"
)

SeqIO.write(orf_record, "orf_protein.fasta", "fasta")

print("ORF protein saved as orf_protein.fasta")

"""“Following ORF-based sequence filtering, the longest predicted ORF (2387 amino acids) was extracted and saved as a protein FASTA file for downstream homology-based comparison.”

#SEQUENCE COMPARISON (HOMOLOGY SEARCH)

PURPOSE OF THIS STEP (VERY IMPORTANT)

To identify whether the selected ORF is similar to any known genes/proteins and to infer its possible function based on evolutionary conservation.

This step answers:

“What is this gene similar to?”

“Does it already exist in other organisms?”
"""

from Bio.Blast import NCBIWWW
from Bio import SeqIO

record = SeqIO.read("orf_protein.fasta", "fasta")

result_handle = NCBIWWW.qblast(
    program="blastp",
    database="nr",
    sequence=record.seq
)

with open("blast_result.xml", "w") as out_handle:
    out_handle.write(result_handle.read())

print("BLAST results saved as blast_result.xml")

from Bio.Blast import NCBIXML

with open("blast_result.xml") as result_handle:
    blast_record = NCBIXML.read(result_handle)

top_alignment = blast_record.alignments[0]
top_hsp = top_alignment.hsps[0]

print("Top hit:", top_alignment.hit_def)
print("E-value:", top_hsp.expect)
print("Identity:", top_hsp.identities)

"""“Homology-based similarity analysis was performed using BLASTp against the NCBI non-redundant protein database. The selected ORF showed a highly significant match (E-value = 0.0) to a putative adhesin protein from Escherichia coli, with near-complete sequence identity, indicating strong evolutionary conservation.”

##Converting interpro_results.tsv file into dataframe
"""

import pandas as pd

interpro_df = pd.read_csv("interpro_results.tsv", sep="\t")

interpro_df.head()

"""“InterProScan analysis identified multiple conserved immunoglobulin-like domains (IPR013783) distributed throughout the protein sequence. These domains are characteristic of bacterial adhesins and are involved in protein–protein and cell–surface interactions. The presence of repeated Ig-like folds strongly supports the functional annotation of the protein as an adhesin.”

“The detection of multiple Ig-like domains suggests a modular architecture, which is typical of large surface-exposed bacterial adhesion proteins.”
"""
