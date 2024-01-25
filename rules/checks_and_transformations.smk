# INPUT CHECKS AND INITIAL DATA TRANSFORMATIONS

# Check for presence of input assembly
rule input_assembly:
    output:
        assembly = "data/assemblies/" + config["assembly"] + ".fasta",

# Check for presenc of input reads
# rule input_reads:
#     output:
#         zipped_fastq_reads = "data/reads/" + config["reads"] + ".fastq.gz",
#         # zipped_fasta_reads = "data/reads/" + config["reads"] + ".fasta.gz",
#         # unzipped_fastq_reads = "data/reads/" + config["reads"] + ".fastq",
#         # unzipped_fasta_reads = "data/reads/" + config["reads"] + ".fasta"
#     run:
#         shell("""
#         if [ ! -e {output[zipped_fastq_reads]} ]; then
#             echo "Error: None of the expected input read files found. Please ensure your reads are in the data/reads/ directory, and that they are in one of the accepted formats."
#             exit 1  # Exit with a non-zero status to indicate failure
#         fi
#         """)
# && [ ! -e {output[zipped_fasta_reads]} ] && [ ! -e {output[unzipped_fastq_reads]} ] && [ ! -e {output[unzipped_fasta_reads]} ]



#         # Define the Snakemake rule for downloading the NT database
# rule download_nt_db:
#     output:
#         os.path.join(config["ncbi_nt_path"], "nt.115.nin")
#     params:
#         ncbi_path = config["ncbi_nt_path"],
        # log = f"{config['assembly']}/logs/download_nt_db.log",
#     shell:
#         """
#         cd {params[ncbi_path]}
#         update_blastdb.pl --passive --decompress nt
#         cd - 2> {params[log]}
#         """
# blast database runs to nt.115. around 350 Gb of storage required

#  # Download the taxdb archive
# perl update_blastdb.pl taxdb
# # Install it in the BLASTDB directory
# gunzip -cd taxdb.tar.gz | (cd $BLASTDB; tar xvf - )

# generate scaffold list
rule scaffold_list:
    input:
        assembly = f"data/assemblies/{config['assembly']}.fasta",
    output:
        all_scaffs = f"{config['assembly']}/logs/scaffold_list.txt"
    run:
        try:
            from Bio import SeqIO
        except ImportError as e:
            print(f"Error importing Bio.SeqIO: {e}")
            traceback.print_exc()

        def generate_header_list(fasta_file, output_file):
            with open(fasta_file, "r") as f:
                # Use SeqIO.parse to iterate over the records in the FASTA file
                headers = [record.id for record in SeqIO.parse(f, "fasta")]
            # Write the list of headers to the output file
            with open(output_file, "w") as out:
                for header in headers:
                    out.write(header.replace(">","") + "\n")
            print("Scaffold names collected.")
        
        # Input FASTA file path
        fasta_file = str(input['assembly'])
    
        # Output text file path
        output_file = str(output['all_scaffs'])
    
        # Generate the list of headers and write it to the output file
        generate_header_list(fasta_file, output_file)



# Generate individual files for each scaffold
all_scaffs = []

with open(f"{config['assembly']}/logs/scaffold_list.txt", "r") as file:
    # Read each line from the file
    for line in file:
        # Strip any leading/trailing whitespace and append the line to the list
        all_scaffs.append(line.strip())


rule splinter_assembly:
    input:
        assembly = "data/assemblies/" + config["assembly"] + ".fasta",
        scaffolds = f"{config['assembly']}/logs/scaffold_list.txt"
    output:
        expand(config["assembly"] + "/outputs/scaffolds/{all_scaffs}.fasta", all_scaffs = all_scaffs)
    params:
        scaffolds = config["assembly"] + "/outputs/scaffolds",
    shell:
        """
        awk '/^>/ {{ if (seq) {{ print header ORS seq > "{params[scaffolds]}/" filename; close("{params[scaffolds]}/" filename) }} header = $0; filename = substr($0, 2); sub(" .*", "", filename); filename = filename ".fasta"; seq = "" }} !/^>/ {{ seq = seq $0 }} END {{ if (seq) {{ print header ORS seq > "{params[scaffolds]}/" filename; close("{params[scaffolds]}/" filename) }} }}' {input}
        """


# Unzip reads if zipped
rule zip_fastq_to_fasta:
    input:
        reads = "data/reads/" + config["reads"] + ".fastq.gz",
    output:
        reads = "data/reads/" + config["reads"] + ".fasta",
    params:
        log = f"{config['assembly']}/logs/zip_fastq_to_fasta.log",
    shell:
        "zcat -c {input[reads]} | seqkit fq2fa | cat > {output} 2> {params[log]}"


# # Unzip reads if zipped
# rule zip_fasta_to_fasta:
#     input:
#         reads = "data/reads/" + config["reads"] + ".{ext}",
#     output:
#         reads = "data/reads/" + config["reads"] + ".fasta",
#     run:
#         if "fastq.gz" in input[0].path:
#             shell("echo 'Processing as FASTQ' > output.txt")
#             # Add your command for processing FASTQ files here
#         elif "fasta.gz" in input[0].path:
#             shell("echo 'Processing as FASTA' > output.txt")
#             # Add your command for processing FASTA files here
#         else:
#             raise ValueError("Unsupported reads file extension. Please use either fastq.gz or fasta.gz.")
