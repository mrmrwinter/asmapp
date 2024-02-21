# INPUT CHECKS AND INITIAL DATA TRANSFORMATIONS

# Check for presence of input assembly
rule input_assembly:
    output:
        assembly = "data/assemblies/" + config["assembly"] + ".fasta",


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

# Unzip reads if zipped
rule fastq_to_fasta:
    input:
        reads = "data/reads/" + config["reads"] + ".fastq",
    output:
        reads = "data/reads/" + config["reads"] + ".fasta",
    params:
        log = f"{config['assembly']}/logs/zip_fastq_to_fasta.log",
    shell:
        "seqkit fq2fa {input[reads]} > {output} 2> {params[log]}"

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



# # Generate individual files for each scaffold
# all_scaffs = []

# with open(f"{config['assembly']}/logs/scaffold_list.txt", "r") as file:
#     # Read each line from the file
#     for line in file:
#         # Strip any leading/trailing whitespace and append the line to the list
#         all_scaffs.append(line.strip())


# rule splinter_assembly:
#     input:
#         assembly = "data/assemblies/" + config["assembly"] + ".fasta",
#         scaffolds = f"{config['assembly']}/logs/scaffold_list.txt"
#     output:
#         expand(config["assembly"] + "/outputs/scaffolds/{all_scaffs}.fasta", all_scaffs = all_scaffs)
#     params:
#         scaffolds = config["assembly"] + "/outputs/scaffolds",
#     shell:
#         """
#         awk '/^>/ {{ if (seq) {{ print header ORS seq > "{params[scaffolds]}/" filename; close("{params[scaffolds]}/" filename) }} header = $0; filename = substr($0, 2); sub(" .*", "", filename); filename = filename ".fasta"; seq = "" }} !/^>/ {{ seq = seq $0 }} END {{ if (seq) {{ print header ORS seq > "{params[scaffolds]}/" filename; close("{params[scaffolds]}/" filename) }} }}' {input}
#         """




