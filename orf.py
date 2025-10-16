
class DNASequence:
#Represent a DNA sequence
    
    def __init__(self, header, gene_id, sequence):
        self.header = header
        self.gene_id = gene_id
        self.sequence = sequence.replace('\n', '').replace(' ', '').upper()
    
    def get_length(self):
        return len(self.sequence)
    
    def find_orfs_in_frame(self, frame_number):
      #Return ORFS using start and stop codon
        START_CODON = 'ATG'
        STOP_CODONS = {'TAA', 'TAG', 'TGA'}
        
        orfs = []
        i = frame_number
        
        # Search for all start codons in this frame
        while i <= len(self.sequence) - 3:
            codon = self.sequence[i:i+3]
            
            # Look for start codon
            if codon == START_CODON:
                start_pos = i
                orf_sequence = START_CODON
                j = i + 3
                
                # Extend until stop codon or end of sequence
                found_stop = False
                while j <= len(self.sequence) - 3:
                    next_codon = self.sequence[j:j+3]
                    orf_sequence += next_codon
                    
                    if next_codon in STOP_CODONS:
                        found_stop = True
                        orfs.append((orf_sequence, start_pos, len(orf_sequence), frame_number + 1))
                        break
                    
                    j += 3
                
                i += 3
            else:
                i += 3
        
        return orfs
    
    def find_all_orfs(self):
#Find all ORFS
        all_orfs = []
        
        for frame in range(3):
            orfs = self.find_orfs_in_frame(frame)
            all_orfs.extend(orfs)
        
        return all_orfs
    
    def find_longest_orf(self):
        #Find the longest ORF
        all_orfs = self.find_all_orfs()
        
        if not all_orfs:
            return None
        
        longest = max(all_orfs, key=lambda x: x[2])
        return longest
    
    def __repr__(self):
        return f"DNASequence({self.gene_id}, length={self.get_length()})"


class FastaDNAParser:
    #parse DNA FASTA files
    
    def __init__(self, fasta_file):
        self.fasta_file = fasta_file
        self.sequences = []
    
    def parse(self):
        #Parse the FASTA file and create DNASequence
        with open(self.fasta_file, 'r') as f:
            current_header = None
            current_gene_id = None
            current_sequence = []
            
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    # Save previous sequence if exists
                    if current_header:
                        seq = ''.join(current_sequence)
                        self.sequences.append(
                            DNASequence(current_header, current_gene_id, seq)
                        )
                    
                    # Parse new header
                    current_header = line[1:].strip()
                    
                    # Extract gene_id from 2nd column
                    parts = current_header.split('|')
                    if len(parts) >= 2:
                        current_gene_id = parts[1]
                    else:
                        current_gene_id = parts[0]
                    
                    current_sequence = []
                else:
                    # Add sequence line
                    if line:
                        current_sequence.append(line)
            
            # Don't forget the last sequence
            if current_header:
                seq = ''.join(current_sequence)
                self.sequences.append(
                    DNASequence(current_header, current_gene_id, seq)
                )
        
        return self.sequences


def main():
    #main function to analyze sequences
    
    print("=" * 70)
    print("ORF FINDER - DNA Sequence Analysis (Problem 3)")
    print("=" * 70)
    
    # Parse FASTA file
    print("\nParsing dna.fasta...")
    parser = FastaDNAParser('dna.fasta')
    sequences = parser.parse()
    print(f"Found {len(sequences)} sequences\n")
    
    # Find longest and shortest sequences
    print("=" * 70)
    print("SEQUENCE LENGTH ANALYSIS")
    print("=" * 70)
    
    if sequences:
        longest_seq = max(sequences, key=lambda x: x.get_length())
        shortest_seq = min(sequences, key=lambda x: x.get_length())
        
        print(f"\nLongest sequence:")
        print(f"  Gene ID: {longest_seq.gene_id}")
        print(f"  Length: {longest_seq.get_length()} bp")
        print(f"  Header: {longest_seq.header[:60]}...")
        
        print(f"\nShortest sequence:")
        print(f"  Gene ID: {shortest_seq.gene_id}")
        print(f"  Length: {shortest_seq.get_length()} bp")
        print(f"  Header: {shortest_seq.header[:60]}...")
    
    # Find longest ORF in each sequence
    print("\n" + "=" * 70)
    print("OPEN READING FRAME (ORF) ANALYSIS")
    print("=" * 70)
    print("\nFinding longest ORF in each sequence...\n")
    
    results = []
    for seq in sequences:
        orf_data = seq.find_longest_orf()
        results.append((seq, orf_data))
        
        print(f"Gene ID: {seq.gene_id}")
        print(f"  Sequence length: {seq.get_length()} bp")
        
        if orf_data:
            orf_seq, start_pos, orf_length, frame = orf_data
            print(f"  Longest ORF length: {orf_length} bp")
            print(f"  Found in frame: {frame}")
            print(f"  Start position: {start_pos}")
            print(f"  ORF sequence (first 60 bp): {orf_seq[:60]}{'...' if len(orf_seq) > 60 else ''}")
        else:
            print(f"  No ORF found")
        print()
    
    # Summary statistics
    print("=" * 70)
    print("SUMMARY STATISTICS")
    print("=" * 70)
    
    orfs_found = [r[1] for r in results if r[1] is not None]
    
    print(f"\nTotal sequences analyzed: {len(sequences)}")
    print(f"Sequences with ORFs: {len(orfs_found)}")
    print(f"Sequences without ORFs: {len(sequences) - len(orfs_found)}")
    
    if orfs_found:
        avg_orf_length = sum(r[2] for r in orfs_found) / len(orfs_found)
        max_orf = max(orfs_found, key=lambda x: x[2])
        min_orf = min(orfs_found, key=lambda x: x[2])
        
        print(f"\nAverage ORF length: {avg_orf_length:.2f} bp")
        print(f"Longest ORF: {max_orf[2]} bp")
        print(f"Shortest ORF: {min_orf[2]} bp")
        
        # Frame distribution
        frame_counts = {1: 0, 2: 0, 3: 0}
        for orf_data in orfs_found:
            frame = orf_data[3]
            frame_counts[frame] += 1
        
        print(f"\nORF distribution by reading frame:")
        for frame in sorted(frame_counts.keys()):
            print(f"  Frame {frame}: {frame_counts[frame]} ORFs")
    
    # Save results to file
    print("\n" + "=" * 70)
    print("Saving results to problem3_answers.txt...")
    
    with open('problem3_answers.txt', 'w') as f:
        f.write("PROBLEM 3 - ORF ANALYSIS RESULTS\n")
        f.write("=" * 70 + "\n\n")
        
        f.write("SEQUENCE LENGTH ANALYSIS:\n")
        f.write("-" * 70 + "\n")
        f.write(f"Longest sequence: {longest_seq.gene_id} ({longest_seq.get_length()} bp)\n")
        f.write(f"Shortest sequence: {shortest_seq.gene_id} ({shortest_seq.get_length()} bp)\n\n")
        
        f.write("LONGEST ORF IN EACH SEQUENCE:\n")
        f.write("-" * 70 + "\n")
        for seq, orf_data in results:
            f.write(f"\nGene ID: {seq.gene_id}\n")
            f.write(f"  Sequence length: {seq.get_length()} bp\n")
            if orf_data:
                orf_seq, start_pos, orf_length, frame = orf_data
                f.write(f"  Longest ORF: {orf_length} bp in frame {frame}\n")
                f.write(f"  Start position: {start_pos}\n")
            else:
                f.write(f"  No ORF found\n")
    
    print("Results saved to problem3_answers.txt")
    print("=" * 70)


if __name__ == "__main__":
    main()
