import gzip
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

subject_fasta = "/mnt/picea/storage/reference/Populus-tremula/v2.2/fasta/Potra02_genome.fasta.gz"
bed_file = "../results/SM_LSM_genes.bed"
output_fasta_file = "../results/SM_LSM_genes.fasta"

# Load the subject genome
with gzip.open(subject_fasta, "rt") as handle:
    subject_genome = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))

# Function to extract sequence based on coordinates and strand
def extract_sequence(contig, start, end, strand):
    if strand == '+':
        return subject_genome[contig].seq[start - 1:end]
    else:
        return subject_genome[contig].seq[start - 1:end].reverse_complement()

# List to store SeqRecord objects
seq_records = []

# Read bed file and extract sequences for each gene
with open(bed_file, "r") as bed_handle:
    for line in bed_handle:
        fields = line.strip().split("\t")
        print(fields)
        contig, start, end, strand, gene = fields[0], int(fields[1]), int(fields[2]), fields[3], fields[4]

        extracted_sequence = extract_sequence(contig, start, end, strand)

        seq_record = SeqRecord(extracted_sequence, id=gene, description=f"{gene}")
        seq_records.append(seq_record)

# Save the combined sequences to a single file
with open(output_fasta_file, "w") as file:
    SeqIO.write(seq_records, file, "fasta")
