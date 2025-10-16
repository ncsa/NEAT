import subprocess
import gzip
import shutil
import yaml
import pysam

from pathlib import Path
from typing import List
from neat.read_simulator.utils import Options, OutputFileWriter
from Bio import bgzf
from Bio.bgzf import BgzfWriter


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

# bacterial_wrapper(salmonella_file, "Salmonella", salmonella_config_file)
# bacterial_wrapper(pneumonia_file, "Pneumonia", pneumonia_config_file)
# bacterial_wrapper(tuberculosis_file, "Tuberculosis", tuberculosis_config_file)



# Stitching outputs back together - based on Keshav's script

# def stitch_outputs(input_files, bacteria_name):
#     dest_file = Path(bacteria_name + "_stitched_reads.fastq")
#     with dest_file.open("wb") as out_f:
#         for input_file in input_files:
#             input_file = Path(input_file)
            
#             if input_file.suffix in {".gz", "bgz"}:
#                 opener = gzip.open
#             else:
#                 opener = open
#             with opener(input_file, "rb") as in_f:
#                 shutil.copyfileobj(in_f, out_f)


# Stitch output fastq files together for each reference file

# orig_salmonella_input = "Outputs/Salmonella/salmonella_config_test.fastq"
# new_salmonella_input = "Outputs/Salmonella/new_Salmonella_config_test.fastq"
# salmonella_inputs = [orig_salmonella_input, new_salmonella_input]
# stitch_outputs(salmonella_inputs, "salmonella")

# orig_pneumonia_input = "Outputs/Pneumonia/pneumonia_config_test.fastq"
# new_pneumonia_input = "Outputs/Pneumonia/new_Pneumonia_config_test.fastq"
# pneumonia_inputs = [orig_pneumonia_input, new_pneumonia_input]
# stitch_outputs(pneumonia_inputs, "pneumonia")

# orig_tuberculosis_input = "Outputs/Tuberculosis/tuberculosis_config_test.fastq"
# new_tuberculosis_input = "Outputs/Tuberculosis/new_Tuberculosis_config_test.fastq"
# tuberculosis_inputs = [orig_tuberculosis_input, new_tuberculosis_input]
# stitch_outputs(tuberculosis_inputs, "tuberculosis")






# Stitching all outputs together - Keshav's script

def concat_fq(input_files: List[Path], dest: BgzfWriter) -> None:
    
    if not input_files:
        return
    
    for input_file in input_files:
        with bgzf.BgzfReader(input_file) as in_f:
            shutil.copyfileobj(in_f, dest)

# def concat_fq(input_files: List[Path], dest: Path) -> None:
    
#     if not input_files:
#         return
    
#     with dest.open("wb") as out_f:
#         for input_file in input_files:
#             input_file = Path(input_file)
            
#             with input_file.open("rb") as in_f:
#                 shutil.copyfileobj(in_f, out_f)

def merge_bam(bams: List[Path], ofw: OutputFileWriter, threads: int) -> None:
    
    if not bams:
        return

    unsorted = ofw.bam.with_suffix(".unsorted.bam")
    pysam.merge("--no-PG", "-@", str(threads), "-f", str(unsorted), *map(str, bams))
    pysam.sort("-@", str(threads), "-o", str(ofw.bam), str(unsorted))
    unsorted.unlink(missing_ok=True)

# def merge_bam(bams: List[Path], dest: Path) -> None:
    
#     if not bams:
#         return

#     unsorted = dest.with_suffix(".unsorted.bam")
#     pysam.merge("--no-PG", "-@", str(unsorted), *map(str, bams))
#     pysam.sort("-@", str(dest), str(unsorted))
#     unsorted.unlink(missing_ok=True)

# def merge_vcf(vcfs: List[Path], dest: Path) -> None:
#     if not vcfs:
#         return
    
#     first, *rest = vcfs
#     shutil.copy(first, dest)

#     with dest.open("ab") as out_f:
#         for vcf in rest:
#             with vcf.open("rb") as fh:
#                 for line in fh:
#                     if not line.startswith(b"#"):
#                         out_f.write(line)

def stitch_all_outputs(ofw: OutputFileWriter, output_files: list[tuple[int, dict[str, Path]]], 
                       threads: int | None = None) -> None:
    
    fq1_list = []
    fq2_list = []
    bam_list = []

    for (thread_idx, file_dict) in output_files:
        if file_dict["fq1"]:
            fq1_list.append(file_dict["fq1"])
        if file_dict["fq2"]:
            fq2_list.append(file_dict["fq2"])
        if file_dict["bam"]:
            bam_list.append(file_dict["bam"])
    
    concat_fq(fq1_list, ofw.files_to_write[ofw.fq1])
    concat_fq(fq2_list, ofw.files_to_write[ofw.fq2])
    merge_bam(bam_list, ofw, threads)

# def stitch_all_outputs(options: Options, thread_options: list[Options]) -> None:
#     fq1_list = []
#     fq2_list = []
#     vcf_list = []
#     bam_list = []

#     for local_ops in thread_options:
#         if local_ops.fq1:
#             fq1_list.append(local_ops.fq1)
#         if local_ops.fq2:
#             fq2_list.append(local_ops.fq2)
#         if local_ops.vcf:
#             vcf_list.append(local_ops.vcf)
#         if local_ops.bam:
#             bam_list.append(local_ops.bam)
    
#     concat_fq(fq1_list, options.fq1)
#     concat_fq(fq2_list, options.fq2)
#     merge_bam(bam_list, options.bam)
#     merge_vcf(vcf_list, options.vcf)


# orig_salmonella_input = Path("Outputs/Salmonella/salmonella_config_test.fastq")
# new_salmonella_input = Path("Outputs/Salmonella/new_Salmonella_config_test.fastq")
# salmonella_inputs = [orig_salmonella_input, new_salmonella_input]
# dest = Path("Outputs/Salmonella/salmonella_stitched.fastq")
# concat_fq(salmonella_inputs, dest)

# orig_salmonella_input = Path("Outputs/Salmonella/salmonella_config_test_golden.bam")
# new_salmonella_input = Path("Outputs/Salmonella/new_Salmonella_config_test_golden.bam")
# salmonella_inputs = [orig_salmonella_input, new_salmonella_input]
# dest = Path("Outputs/Salmonella/salmonella_stitched_golden.bam")
# merge_bam(salmonella_inputs, dest)

# orig_salmonella_input = Path("Outputs/Salmonella/salmonella_config_test_golden.vcf.gz")
# new_salmonella_input = Path("Outputs/Salmonella/new_Salmonella_config_test_golden.vcf.gz")
# salmonella_inputs = [orig_salmonella_input, new_salmonella_input]
# dest = Path("Outputs/Salmonella/salmonella_stitched_golden.vcf.gz")
# merge_vcf(salmonella_inputs, dest)

# output_path = Path("Outputs/Salmonella/all_stitched/final")
# config_file = salmonella_config_file
# reference = salmonella_file
# fq1 = Path("Outputs/Salmonella/salmonella_config_test.fastq")
# fq2 = Path("")
# vcf = Path("Outputs/Salmonella/salmonella_config_test_golden.vcf.gz")
# bam = Path("Outputs/Salmonella/salmonella_config_test_golden.bam")


# options = Options(reference)
# bam_header = {}

# ofw = OutputFileWriter(options, bam_header)
# output_files = [(1, {"fq1": Path("Outputs/Salmonella/salmonella_config_test.fastq"),
#                      "fq2": Path(""),
#                      "bam": Path("Outputs/Salmonella/salmonella_config_test_golden.bam")})]

# stitch_all_outputs(ofw, output_files, 1)

# output_files: list[tuple[int, str, dict[str, Path]]]