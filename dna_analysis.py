
#Define Initial DNA Sequence
dna_sequence = """ATTGGGGAGGAGGCGAGTTGAGCGGCGGCAGTTCGCCTGCGTGCGCTGCGCGG
CGTCGACATCTGATCCGCACCATGGAAATCCCCGCTCAATCTTTGGAGCAGGGAT
GCGGGGCGATCAAGATGGGGATGCGGGATGGGGGCGACGGTGTATTTCCGCCAG
AAGATTTCGCCGCGGGAGCTCGCGGTGCGTACGTGCATGTTCAAACGCACGGTG
CGCGCATGGCAGTGGCAGACTGATCAACGCAGCTGGAAGCATCCGAAGCGCGCG
GGCACGCGTGTCCTCGACGCGTGGCCTCACATGCTGTCGGGTCGGTTCAAGACC
GAAAGCCACCGACCGACGCGCGAGCAATGCGCTACGCGGATCGCGTTCGACACG
AGCCGCGCGCGAGGCAAGGCCGACGTATTCGATCTTCCAGAGGAAGCCTATTGG
CTCGAGTCGTAGTGCTCGATATGGTAGAGCAACATGAATCCCGGGCTAAGTACAA
GAAGTAACCCGGCAACGAGTGAGATTGCGACGAATAAACGCTTCACCATGATCGC
GCTCCTGAGTTGGTTGAGGTGAATTGGAAAGTCGATTCCTGGGGGATCATTCCCG
GCAAGGCGCGCAATCCCCGCATTGTTCTCAAGATCGCAACGCGATTCGTCAGGCC
GATCTTCATGGGGTGTCTCGCTGGTAGTGATTCCGTCGTGGCCCGCGCATGTGCA
TGACGGCATCCGGGGAG"""

# Remove newlines 
dna_sequence = dna_sequence.replace('\n', '')  
#Remove whitespace
dna_sequence = dna_sequence.replace(' ', '')   

# 1. Calculate sequence length
sequence_length = len(dna_sequence)
print("=" * 60)
print("DNA Sequence Analysis")
print("=" * 60)
print(f"\n1. Sequence Length: {sequence_length} nucleotides\n")

# 2. Calculate triplet frequencies
print("=" * 60)
print("2. Triplet Frequencies")
print("=" * 60)

triplet_counts = {}
total_triplets = 0

# Count triplets
for i in range(0, len(dna_sequence) - 2, 3):
    triplet = dna_sequence[i:i+3]
    if len(triplet) == 3:
        triplet_counts[triplet] = triplet_counts.get(triplet, 0) + 1
        total_triplets += 1

# Sort triplets alphabetically
sorted_triplets = sorted(triplet_counts.items())

# Print table header
print(f"\n{'Triplet':<10} {'Count':<10} {'Frequency':<15}")
print("-" * 35)

# Print triplet frequencies
for triplet, count in sorted_triplets:
    frequency = count / total_triplets
    print(f"{triplet:<10} {count:<10} {frequency:.6f}")

print(f"\nTotal triplets analyzed: {total_triplets}")

# 3. Find genetic complement
print("\n" + "=" * 60)
print("3. Genetic Complement")
print("=" * 60)

complement_map = {
    'A': 'T',
    'T': 'A',
    'G': 'C',
    'C': 'G'
}

#Convert DNA Sequence to Complement
complement_list = []
for base in dna_sequence:
    complement_list.append(complement_map[base])
complement_sequence = ''.join(complement_list)

print(f"\nOriginal sequence (first 60 bases):")
print(dna_sequence[:60])
print(f"\nComplement sequence (first 60 bases):")
print(complement_sequence[:60])
print(f"\nFull complement sequence:")
print(complement_sequence)

# Save results to file
f = open('dna_analysis_results.txt', 'w')
f.write("DNA SEQUENCE ANALYSIS RESULTS\n")
f.write("=" * 60 + "\n\n")
f.write(f"1. Sequence Length: {sequence_length} nucleotides\n\n")

f.write("2. Triplet Frequencies:\n")
f.write(f"{'Triplet':<10} {'Count':<10} {'Frequency':<15}\n")
f.write("-" * 35 + "\n")
for triplet, count in sorted_triplets:
    frequency = count / total_triplets
    f.write(f"{triplet:<10} {count:<10} {frequency:.6f}\n")
f.write(f"\nTotal triplets: {total_triplets}\n\n")

f.write("3. Genetic Complement:\n")
f.write(complement_sequence + "\n")
f.close()  # Must remember to close!

print("\n" + "=" * 60)
print("Results saved to 'dna_analysis_results.txt'")
print("=" * 60)
