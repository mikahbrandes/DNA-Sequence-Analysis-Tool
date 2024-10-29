# Mikah's DNA sequence analysis tool 

# Analyzes nucleuotide composition and identifies coding regions (ORFs)

from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load DNA sequence data
def load_fasta(filepath):
    with open(filepath, "r") as file:
        sequences = list(SeqIO.parse(file, "fasta"))
    return sequences

sequences = load_fasta("data/nrf2_seq.fasta")

# Function to calculate nucleotide composition
def nucleotide_composition(sequence):
    composition = {
        "A": sequence.count("A"),
        "T": sequence.count("T"),
        "G": sequence.count("G"),
        "C": sequence.count("C")
    }
    return composition

# Test on the first sequence
composition = nucleotide_composition(str(sequences[0].seq))
print(composition)
{'A': 1773, 'T': 1825, 'G': 1291, 'C': 919}
# Plot nucleotide composition
def plot_composition(composition):
    df = pd.DataFrame(list(composition.items()), columns=["Nucleotide", "Count"])
    sns.barplot(x="Nucleotide", y="Count", data=df)
    plt.title("Nucleotide Composition")
    plt.show()

plot_composition(composition)
![image](https://github.com/user-attachments/assets/0eb679a2-f549-4dc0-bd19-85a4d8f68d02)
# Function to find ORFs
def find_orfs(sequence, min_length=100):
    orfs = []
    for frame in range(3):
        for i in range(frame, len(sequence) - min_length, 3):
            if sequence[i:i+3] == "ATG":
                for j in range(i, len(sequence) - min_length, 3):
                    if sequence[j:j+3] in ["TAA", "TAG", "TGA"]:
                        orfs.append(sequence[i:j+3])
                        break
    return orfs

orfs = find_orfs(str(sequences[0].seq))
print(f"Number of ORFs found: {len(orfs)}")
Number of ORFs found: 110
