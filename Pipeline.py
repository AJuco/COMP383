"""
Aaron Juco
Professor Heather Wheeler
COMP 383-001 Spring 2025
02/26/2025
Python Pipeline Project
"""
import argparse
import sys
import os
import subprocess
from Bio import Entrez, SeqIO


# Copy-pasted check_arg statement
def check_arg(args=None):
	parser = argparse.ArgumentParser(description="Outputs number of reads that meet quality score")
	parser.add_argument('-i', "--input",
		help="path to input file",
		required="True"
    )
	parser.add_argument('-o', "--output",
        help="output file name",
		required="True"
    )
	return parser.parse_args(args)

#retrieve command line arguments and assign to variables
args = check_arg(sys.argv[1:])
infile = args.input
PipelineProject = args.output


"""
Step #1- Making the directory and retrieving the transcriptomes and conversion to fastq files using wget
"""
# Specifies the directory name
directory_name = "PipelineProject_<Aaron>_<Juco>"

# Creates the directory
try:
    os.mkdir(directory_name)
    print(f"Directory '{directory_name}' created successfully.")
except FileExistsError:
    print(f"Directory '{directory_name}' already exists.")
except PermissionError:
    print(f"Permission denied: Unable to create '{directory_name}'.")
except Exception as e:
    print(f"An error occurred: {e}")

# List of SRR numbers given from one patient donor 2 and 6 dpi
srr_numbers = ["SRR5660030", "SRR5660033"]

# Create a directory for the FASTQ files
os.makedirs("fastq_data", exist_ok=True)
os.chdir("fastq_data")

def download_fastq(srr_id):
    """Download paired-end FASTQ files using wget."""
    # Construct URLs for paired-end FASTQ files
    url_1 = f"s3://sra-pub-src-6/{srr_id}/V1_S4_L001_R2_001.fastq"
    url_2 = f"s3://sra-pub-src-6/{srr_id}/V1_S4_L001_R2_001.fastq"
    
    subprocess.run(["wget", url_1], check=True)
    subprocess.run(["wget", url_2], check=True)

# Loops through each SRR number and download
for srr in srr_numbers:
    download_fastq(srr)


"""
Step # 2: Outputs the number of reads in each transcriptome before and after bowtie2 mapping 
into the log file called PipelineProject.log found in the directory to be named "PipelineProject__", 
which is preemptively created for you in the script already in lines 37-48
"""

# Initialization of Variables
genome_accession = "NC_006273.2" # HCMV
genome_fasta = f"{genome_accession}.fasta"
bowtie_index = "HCMV_index"
fastq_dir = "fastq_data"
output_dir = "bowtie2_output"
log_file = "mapping_log.txt"
srr_numbers = ["SRR5577751", "SRR5577754"]

# Downloads HCMV Genome
def download_genome():
    """Download the HCMV reference genome using wget."""
    if not os.path.exists(genome_fasta):
        url = f"https://www.ncbi.nlm.nih.gov/search/api/sequence/{genome_accession}/?report=fasta"
        print("Downloading HCMV genome...")
        subprocess.run(["wget", "-O", genome_fasta, url], check=True)

# Creates Bowtie2 Index
def create_bowtie2_index():
    """Create Bowtie2 index for the HCMV genome."""
    print("Creating Bowtie2 index...")
    subprocess.run([
        "bowtie2-build", genome_fasta, bowtie_index
    ], check=True)

# Map Reads to HCMV Index and Filter Mapped Reads
def map_and_filter_reads(srr_id):
    """Map reads to HCMV index and filter mapped reads."""
    # Input FASTQ files
    fastq1 = f"{fastq_dir}/{srr_id}_1.fastq.gz"
    fastq2 = f"{fastq_dir}/{srr_id}_2.fastq.gz"

    # Output files
    sam_file = f"{output_dir}/{srr_id}.sam"
    bam_file = f"{output_dir}/{srr_id}.bam"
    filtered_fastq1 = f"{output_dir}/{srr_id}_mapped_1.fastq"
    filtered_fastq2 = f"{output_dir}/{srr_id}_mapped_2.fastq"

    # Create output directory if not exists
    os.makedirs(output_dir, exist_ok=True)

    # Count reads before mapping
    read_count_before = count_reads(fastq1, fastq2)

    # Map reads using Bowtie2
    print(f"Mapping {srr_id} reads to HCMV index...")
    subprocess.run([
        "bowtie2",
        "-x", bowtie_index,
        "-1", fastq1,
        "-2", fastq2,
        "-S", sam_file
    ], check=True)

    # Convert SAM to BAM and filter mapped reads
    print(f"Filtering mapped reads for {srr_id}...")
    subprocess.run(["samtools", "view", "-Sb", "-F", "4", "-o", bam_file, sam_file], check=True)
    subprocess.run(["samtools", "fastq", "-1", filtered_fastq1, "-2", filtered_fastq2, bam_file], check=True)

    # Count reads after mapping
    read_count_after = count_reads(filtered_fastq1, filtered_fastq2)

    # Log the read counts
    with open(log_file, "a") as log:
        log.write(f"{srr_id} had {read_count_before} read pairs before Bowtie2 filtering and {read_count_after} read pairs after.\n")

# Helper Function: Count Reads
def count_reads(fastq1, fastq2):
    """Count the number of read pairs in the FASTQ files."""
    count1 = int(subprocess.check_output(["zcat", fastq1, "|", "wc", "-l"]).decode().split()[0]) // 4
    count2 = int(subprocess.check_output(["zcat", fastq2, "|", "wc", "-l"]).decode().split()[0]) // 4
    return min(count1, count2)

# Calls method to execute code
download_genome()
create_bowtie2_index()

# Loops through each SRR and process
for srr in srr_numbers:
    map_and_filter_reads(srr)

print("Bowtie2 mapping and filtering completed. Check the log file for details.")

"""
Step 3: Assemble the transcriptomes together to produce 1 assembly via SPAdes and a k-mer size of 99
"""

# Directories and files
output_dir = "spades_output"
bowtie_output_dir = PipelineProject


def assemble_spades():
    # Construct SPAdes command
    spades_command = [
        "spades.py",
        "-k", "99",
        "--only-assembler",
        "-o", output_dir
    ]

    # Add input files to the command
    for r1, r2 in zip(srr_numbers, srr_numbers):
        spades_command.extend(["-1", r1, "-2", r2])

    # Create output directory if not exists
    os.makedirs(output_dir, exist_ok=True)

    # Run SPAdes
    print("Running SPAdes assembly...")
    subprocess.run(spades_command, check=True)

    # Log the SPAdes command
    spades_command_str = " ".join(spades_command)
    with open(PipelineProject, "a") as log:
        PipelineProject.write(f"\nSPAdes command used: {spades_command_str}\n")

# Runs the SPAdes assembly
assemble_spades()

print("SPAdes assembly complete. Check the log file for the command used.")

"""
Step 4: Output the number of contigs with a length greater than 100
into the PipelineProject.log file as well as the length of the assembly
"""
from Bio import SeqIO

# Input and output files
contigs_file = "spades_output/contigs.fasta"
log_file = PipelineProject

# This method calculates the number and total length of contigs greater than 1000 bp
def analyze_contigs():
    
    #   Initializes contig and lengthcount
    contig_count = 0
    total_length = 0

    # Parse the contigs FASTA file
    for record in SeqIO.parse(contigs_file, "fasta"):
        contig_length = len(record.seq)
        # Check if the contig is longer than 1000 bp
        if contig_length > 1000:
            contig_count += 1
            total_length += contig_length

    # Write results to the log file
    with open(log_file, "a") as log:
        log.write(f"\nThere are {contig_count} contigs > 1000 bp in the assembly.\n")
        log.write(f"There are {total_length} bp in the assembly.\n")

    print("Contig analysis complete. Results written to log file.")

# Run the contig analysis
analyze_contigs()

"""
Step 5: Retrieve the longest contig from the spades assembly using blastn and write the top 10 hits': 
subject accession (sacc), percent identity (pident), alignment length (length), start of alignment in query (qstart), 
end of alignment in query (qend), start of alignment in subject (sstart), end of alignment in subject (send), bitscore,
e-value,and subject title (stitle)
"""

# NCBI Entrez settings
Entrez.email = "ajuco@luc.edu" # Change depending on your email
betaherpes_db = "Betaherpesvirinae_db.fasta"
local_db_name = "Betaherpesvirinae_db"

# This method downloads sequences from the Betaherpesvirinae subfamily 
def download_betaherpes_sequences():
    print("Searching for Betaherpesvirinae sequences...")
    query = "Betaherpesvirinae[Organism] AND srcdb_refseq[PROP]"
    handle = Entrez.esearch(db="nucleotide", term=query, retmax=100)
    record = Entrez.read(handle)
    handle.close()

    id_list = record["IdList"]
    print(f"Found {len(id_list)} sequences. Downloading...")

    # Fetch the sequences
    handle = Entrez.efetch(db="nucleotide", id=id_list, rettype="fasta", retmode="text")
    with open(betaherpes_db, "w") as out_f:
        out_f.write(handle.read())
    handle.close()

# Creates a local BLAST database for Betaherpesvirinae
def create_local_blast_db():
    print("Creating local BLAST database...")
    subprocess.run([
        "makeblastdb",
        "-in", betaherpes_db,
        "-dbtype", "nucl",
        "-out", local_db_name
    ], check=True)

# Run the workflow
download_betaherpes_sequences()
create_local_blast_db()

contigs_file = "spades_output/contigs.fasta"
longest_contig_file = "longest_contig.fasta"

def get_longest_contig():
    """Retrieve the longest contig from the assembly."""
    longest_contig = None
    max_length = 0

    # Find the longest contig
    for record in SeqIO.parse(contigs_file, "fasta"):
        if len(record.seq) > max_length:
            longest_contig = record
            max_length = len(record.seq)

    # Write the longest contig to a file
    if longest_contig:
        SeqIO.write(longest_contig, longest_contig_file, "fasta")
        print(f"Longest contig is {max_length} bp long and saved to {longest_contig_file}.")
    else:
        print("No contigs found!")

# Run the function
get_longest_contig()

from Bio.Blast import NCBIXML

blast_output = "blast_output.xml"
log_file = "mapping_log.txt"

def run_blast():
    """Run BLASTn with the longest contig against the local Betaherpesvirinae DB."""
    print("Running BLASTn...")
    subprocess.run([
        "blastn",
        "-query", longest_contig_file,
        "-db", local_db_name,
        "-out", blast_output,
        "-outfmt", "5",   # XML output for easy parsing
        "-max_target_seqs", "10",
        "-best_hit_overhang", "0.1",
        "-best_hit_score_edge", "0.1"
    ], check=True)

def parse_blast_results():
    """Parse BLAST XML output and log the top 10 hits."""
    with open(blast_output) as result_handle:
        blast_records = NCBIXML.parse(result_handle)

        with open(log_file, "a") as log:
            log.write("\nTop 10 BLAST hits for the longest contig:\n")
            for record in blast_records:
                for alignment in record.alignments[:10]:
                    for hsp in alignment.hsps:
                        log.write(f"Subject accession: {alignment.accession}\n")
                        log.write(f"Percent identity: {hsp.identities / hsp.align_length * 100:.2f}\n")
                        log.write(f"Alignment length: {hsp.align_length}\n")
                        log.write(f"Start of alignment in query: {hsp.query_start}\n")
                        log.write(f"End of alignment in query: {hsp.query_end}\n")
                        log.write(f"Start of alignment in subject: {hsp.sbjct_start}\n")
                        log.write(f"End of alignment in subject: {hsp.sbjct_end}\n")
                        log.write(f"Bit score: {hsp.bits}\n")
                        log.write(f"E-value: {hsp.expect}\n")
                        log.write(f"Subject Title: {alignment.title}\n")
                        log.write("\n")

    print("BLAST analysis complete. Results written to log file.")

# Run BLAST and parse the results
run_blast()
parse_blast_results()


