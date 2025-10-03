class GeneInfo:
    """Class to represent gene information from protein-coding_gene.txt"""
    
    def __init__(self, gene_name):
        self.gene_name = gene_name
        self.synonyms = []
        self.location = ""
    
    def add_synonyms(self, synonyms_str):
        # Add synonyms from comma-separated string
        if synonyms_str and synonyms_str.strip():
            # Split by comma and clean up each synonym
            syn_list = [s.strip() for s in synonyms_str.split(',')]
            self.synonyms = [s for s in syn_list if s]
    
    def set_location(self, location):
        # Set the gene location
        self.location = location.strip() if location else ""
    
    def get_synonyms_string(self):
        # Return synonyms as comma-separated string
        return ','.join(self.synonyms) if self.synonyms else "N/A"
    
    def __repr__(self):
        return f"GeneInfo({self.gene_name}, synonyms={self.synonyms}, location={self.location})"


class ProteinSequence:
    # Class to represent a protein sequence from FASTA file
    
    def __init__(self, gene_name, sequence):
        self.gene_name = gene_name
        self.sequence = sequence.replace('\n', '').replace(' ', '')
        self.gene_info = None
    
    def attach_gene_info(self, gene_info):
        # Attach gene information to this sequence
        self.gene_info = gene_info
    
    def get_formatted_header(self):
        # Return formatted FASTA header with gene info
        if self.gene_info:
            synonyms = self.gene_info.get_synonyms_string()
            location = self.gene_info.location if self.gene_info.location else "Unknown"
            return f">{self.gene_name}|{synonyms}|{location}"
        else:
            return f">{self.gene_name}|N/A|Unknown"
    
    def get_formatted_sequence(self, line_length=60):
        # Return sequence formatted in blocks of specified length
        formatted_lines = []
        for i in range(0, len(self.sequence), line_length):
            formatted_lines.append(self.sequence[i:i+line_length])
        return '\n'.join(formatted_lines)
    
    def __repr__(self):
        return f"ProteinSequence({self.gene_name}, length={len(self.sequence)})"


class FastaParser:
    # Class to parse FASTA files
    
    def __init__(self, fasta_file):
        self.fasta_file = fasta_file
        self.sequences = []
    
    def parse(self):
        # Parse the FASTA file and create ProteinSequence objects
        with open(self.fasta_file, 'r') as f:
            current_gene = None
            current_sequence = []
            
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    # Save previous sequence if exists
                    if current_gene:
                        seq = ''.join(current_sequence)
                        self.sequences.append(ProteinSequence(current_gene, seq))
                    
                    # Start new sequence - gene name is everything after '>'
                    current_gene = line[1:].strip()
                    current_sequence = []
                else:
                    # Add sequence line
                    if line:  # Skip empty lines
                        current_sequence.append(line)
            
            # Don't forget the last sequence
            if current_gene:
                seq = ''.join(current_sequence)
                self.sequences.append(ProteinSequence(current_gene, seq))
        
        return self.sequences


class GeneInfoParser:
    # Class to parse gene information file
    
    def __init__(self, gene_info_file):
        self.gene_info_file = gene_info_file
        self.gene_dict = {}
    
    def parse(self):
        # Parse the gene information file and create GeneInfo objects
        with open(self.gene_info_file, 'r') as f:
            # Skip header line
            header = f.readline()
            
            for line in f:
                line = line.strip()
                if not line:
                    continue
                
                # Split by tab
                parts = line.split('\t')
                
                # Need at least 9 columns based on your format
                if len(parts) >= 8:
                    # Column indices (0-based):
                    # 1 = Approved Symbol (gene name)
                    # 6 = Synonyms
                    # 7 = Chromosome (location)
                    
                    gene_name = parts[1].strip()
                    synonyms_str = parts[6].strip() if len(parts) > 6 else ""
                    location = parts[7].strip() if len(parts) > 7 else ""
                    
                    gene_info = GeneInfo(gene_name)
                    gene_info.add_synonyms(synonyms_str)
                    gene_info.set_location(location)
                    
                    self.gene_dict[gene_name] = gene_info
        
        return self.gene_dict


class FastaWriter:
    # Class to write formatted FASTA files
    
    def __init__(self, output_file):
        self.output_file = output_file
    
    def write(self, sequences):
        # Write sequences to output file
        with open(self.output_file, 'w') as f:
            for seq in sequences:
                f.write(seq.get_formatted_header() + '\n')
                f.write(seq.get_formatted_sequence() + '\n')


def main():
    # Main function to run the parsing pipeline
    
    print("=" * 70)
    print("FASTA PARSER - Problem 2")
    print("=" * 70)
    
    # Parse FASTA file
    print("\nParsing human.fa...")
    fasta_parser = FastaParser('human.fa')
    sequences = fasta_parser.parse()
    print(f"Found {len(sequences)} sequences")
    
    # Parse gene information
    print("\nParsing protein-coding_gene.txt...")
    gene_parser = GeneInfoParser('protein-coding_gene.txt')
    gene_dict = gene_parser.parse()
    print(f"Found information for {len(gene_dict)} genes")
    
    # Attach gene information to sequences
    print("\nMatching sequences with gene information...")
    matched = 0
    unmatched = 0
    for seq in sequences:
        if seq.gene_name in gene_dict:
            seq.attach_gene_info(gene_dict[seq.gene_name])
            print(f"  ✓ Matched: {seq.gene_name}")
            matched += 1
        else:
            print(f"  ✗ Warning: No info found for {seq.gene_name}")
            unmatched += 1
    
    print(f"\nMatching summary: {matched} matched, {unmatched} unmatched")
    
    # Write output
    print("\nWriting output to human_annotated.fa...")
    writer = FastaWriter('human_annotated.fa')
    writer.write(sequences)
    
    print("\n" + "=" * 70)
    print("COMPLETE! Output written to human_annotated.fa")
    print("=" * 70)
    
    # Print detailed summary
    print("\nDETAILED SUMMARY:")
    print("-" * 70)
    for seq in sequences:
        print(f"\nGene: {seq.gene_name}")
        print(f"  Sequence length: {len(seq.sequence)} amino acids")
        if seq.gene_info:
            print(f"  Synonyms: {seq.gene_info.get_synonyms_string()}")
            print(f"  Location: {seq.gene_info.location}")
        else:
            print(f"  No gene information found")


if __name__ == "__main__":
    main()



