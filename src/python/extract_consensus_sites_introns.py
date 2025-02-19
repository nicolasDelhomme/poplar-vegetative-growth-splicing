import gzip

# Read genome file
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


def extract_sequences(intron_file, genome_file, donor_output, acceptor_output):
  donor_records = []
  acceptor_records = []
  
  genome_sequences = read_genome_file(genome_file)

  with gzip.open(intron_file, 'rb') as f:
    for line in f:
      line = line.decode()
      if not line.startswith('#'):
        parts = line.strip().split('\t')
        chr_name = parts[0]
        start = int(parts[3])
        end = int(parts[4])
        strand = parts[6]
        parent = parts[8].split('=')[1]
        #print(parts)
        #print(chr_name)
        #print(start)
        #print(end)
        #print(strand)
        #print(parent)
        
        donor_site = start if strand == "+" else end
        acceptor_site = end if strand == "+" else start
        
        donor_start = donor_site - 4 if strand == "+" else donor_site - 8 
        donor_end = donor_site + 8 if strand == "+" else donor_site + 4 
        acceptor_start = acceptor_site - 12 if strand == "+" else acceptor_site - 2 
        acceptor_end = acceptor_site + 2 if strand == "+" else acceptor_site + 12

        chromosome_seq = genome_sequences.get(chr_name, None)

        if chromosome_seq:
            donor_seq = chromosome_seq[donor_start - 1:donor_end]
            acceptor_seq = chromosome_seq[acceptor_start - 1:acceptor_end]
            #print("5': "+donor_seq)
            #print("3': "+acceptor_seq)
            if donor_seq:
                donor_records.append(f">{chr_name}:{start}:donor:{strand}:{parent}\n{donor_seq}\n")
            if acceptor_seq:
                acceptor_records.append(f">{chr_name}:{end}:acceptor:{strand}:{parent}\n{acceptor_seq}\n")

            # Write donor and acceptor sequences to FASTA files
    with open(donor_output, 'w') as donor_file:
        donor_file.writelines(donor_records)
    with open(acceptor_output, 'w') as acceptor_file:
        acceptor_file.writelines(acceptor_records)

    donor_file.close()
    acceptor_file.close()

    print(f"Sequences saved to {donor_output} and {acceptor_output}.")

# Define input and output file paths
intron_file = "../../reference/gff/Potra02_genes_intron-only.gff.gz"
genome_file = "../../reference/fasta/Potra02_genome_hardmasked.fasta"
donor_output = "intron_consensus_donor_sequences.fa"
acceptor_output = "intron_consensus_acceptor_sequences.fa"

# Extract sequences and save to FASTA files
extract_sequences(intron_file, genome_file, donor_output, acceptor_output)
