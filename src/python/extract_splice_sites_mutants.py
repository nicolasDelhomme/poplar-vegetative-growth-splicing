import pandas as pd

# Read reference genome
def read_genome_file(genome_file):
    genome_sequences = {}
    current_chr = None
    with open(genome_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                current_chr = line.strip()[1:]
                genome_sequences[current_chr] = []
            else:
                if current_chr:
                    genome_sequences[current_chr].append(line.strip())
    return {chr_name: ''.join(seq_lines) for chr_name, seq_lines in genome_sequences.items()}

# Extract nucleotide sequences around 5' donor and 3' acceptor splice sites
def extract_sites(excel_file, genome_file, donor_output, acceptor_output):
    donor_records = []
    acceptor_records = []
    
    df = pd.read_excel(excel_file)
    genome_sequences = read_genome_file(genome_file)

    for index, row in df.iterrows():
        event = row['Event']
        if event.startswith('IR'):
            region = row['Region']
            chr_name, coords = region.split(':')
            donor, acceptor = map(int, coords.split('-'))
            #print(region)
            #print(chr_name)
            #print(coords)
            #print(donor)
            #print(acceptor)
            # Extracting 4 bp upstream and 8 bp downstream of the donor coordinate
            donor_region_start = max(1, donor - 4)
            donor_region_end = donor + 8

            # Extracting 12 bp upstream and 2 bp downstream of the acceptor coordinate
            acceptor_region_start = max(1, acceptor - 12)
            acceptor_region_end = acceptor + 2

            
            chromosome_seq = genome_sequences.get(chr_name, None)
            
            if chromosome_seq:
                # Extract sequences
                donor_sequence = chromosome_seq[donor_region_start-1:donor_region_end]
                acceptor_sequence = chromosome_seq[acceptor_region_start-1:acceptor_region_end]
            if donor_sequence:
                donor_records.append(f">{chr_name}:{donor}:donor\n{donor_sequence}\n")
            if acceptor_sequence:
                acceptor_records.append(f">{chr_name}:{acceptor}:acceptor\n{acceptor_sequence}\n")

            # Write donor and acceptor sequences to FASTA files
    with open(donor_output, 'w') as donor_file:
        donor_file.writelines(donor_records)
    with open(acceptor_output, 'w') as acceptor_file:
        acceptor_file.writelines(acceptor_records)

    donor_file.close()
    acceptor_file.close()


# Read and save to files
excel_file_line3 = '../../doc/line3_control_stranded.xlsx'
excel_file_line26 = '../../doc/line26_control_stranded.xlsx'
genome_file = "../../reference/fasta/Potra02_genome_hardmasked.fasta"
line3_donor_output = "Line3_IR_donor_sequences.fa"
line3_acceptor_output = "Line3_IR_acceptor_sequences.fa"
line26_donor_output = "Line26_IR_donor_sequences.fa"
line26_acceptor_output = "Line26_IR_acceptor_sequences.fa"

extract_sites(excel_file_line3, genome_file, line3_donor_output, line3_acceptor_output)

extract_sites(excel_file_line26, genome_file, line26_donor_output, line26_acceptor_output)
