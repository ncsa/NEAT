import subprocess


# Rearranges the bacterial chromosome by wrapping it around

def wrapper(seq):
    length = len(seq)

    if (length % 2 == 0):
        half_index = length // 2
    else:
        half_index = (length // 2) + 1
    
    first_half = seq[:half_index]
    second_half = seq[half_index:]

    new_seq = second_half + first_half

    return new_seq


# Writes the newly rearranged chromosome's sequence to a new fasta file

def write_fasta_file(new_seq, bacteria_name, fasta_header):
    fasta_file_name = f"wrapped_{bacteria_name}.fna"
    fasta_file = open(fasta_file_name, "w")

    fasta_file.write(fasta_header + "\n" + new_seq)

    fasta_file.close()

    return fasta_file_name


# Writes a yml configuration file for the newly rearranged chromosome's fasta sequence
# This uses the default parameters for NEAT

def write_config_file(ref_config_file, rearranged_seq_file, bacteria_name):
    new_config_file_name = f"new_{bacteria_name}_config_test.yml"

    with open(ref_config_file, 'r') as ref_file, open(new_config_file_name, 'w') as new_file:
        
        for line in ref_file:
            if line.find("reference:") == -1:
                new_file.write(line)
            else:
                new_file.write(f"reference: {rearranged_seq_file}")


    ref_file.close()
    new_file.close()

    return new_config_file_name


# Runs the NEAT read simulator using the given config file

def run_neat(config_file):
    subprocess.run(["neat", "read-simulator", "-c", config_file, "-o", config_file])


# General function for bacterial wrapper that calls all of the functions defined above

def bacterial_wrapper(reference_file, bacteria_name, ref_config_file):
    
    orig_seq = ""

    f = open(reference_file)
    fasta_header = f.readline().strip()

    for line in f:
        if line[0] != ">":
            orig_seq += line
        elif line.find("plasmid") != -1: # exclude plasmids from the sequence to be rearranged
            break
        
    f.close()

    full_orig_seq = orig_seq.replace("\n", "")

    rearranged_seq = wrapper(full_orig_seq)
    rearranged_seq_file = write_fasta_file(rearranged_seq, bacteria_name, fasta_header)
    new_config_file = write_config_file(ref_config_file, rearranged_seq_file, bacteria_name)

    run_neat(ref_config_file)
    run_neat(new_config_file)



# Read in bacterial reference genome files

salmonella_file = "BacterialWrapperScriptReferenceFiles/Salmonella/GCF_000006945.2_ASM694v2_genomic.fna"
pneumonia_file = "BacterialWrapperScriptReferenceFiles/Pneumonia/GCF_001457635.1_NCTC7465_genomic.fna"
tuberculosis_file = "BacterialWrapperScriptReferenceFiles/Tuberculosis/GCF_000195955.2_ASM19595v2_genomic.fna"



# Set config files for each reference file

salmonella_config_file = "Outputs/Salmonella/salmonella_config_test.yml"
pneumonia_config_file = "Outputs/Pneumonia/pneumonia_config_test.yml"
tuberculosis_config_file = "Outputs/Tuberculosis/tuberculosis_config_test.yml"


# Call the bacterial wrapper script on each reference file

bacterial_wrapper(salmonella_file, "Salmonella", salmonella_config_file)
bacterial_wrapper(pneumonia_file, "Pneumonia", pneumonia_config_file)
bacterial_wrapper(tuberculosis_file, "Tuberculosis", tuberculosis_config_file)