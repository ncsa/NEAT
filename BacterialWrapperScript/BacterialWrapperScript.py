import subprocess
import gzip
import shutil
import yaml
import pysam
import unittest
import os

from pathlib import Path
from typing import List
from Bio import bgzf
from Bio.bgzf import BgzfWriter, BgzfReader


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


# Splits the coverage in half for the reference config file
# Writes a yml configuration file for the newly rearranged chromosome's fasta sequence
# This uses the default parameters for NEAT

def write_config_file(ref_config_file, rearranged_seq_file, bacteria_name):
    new_config_file_name = f"test_new_{bacteria_name}_config_test.yml"
    old_config_file_name = f"test_{bacteria_name}_config_test.yml"

    with open(ref_config_file, 'r') as ref_file, open(new_config_file_name, 'w') as new_file, open(old_config_file_name, 'w') as old_file:
        
        for line in ref_file:
            if line.find("reference:") != -1:
                new_file.write(f"reference: {rearranged_seq_file}")
                old_file.write(line)
            elif line.find("coverage:") != -1:
                if line.strip() == "coverage: .":
                    new_coverage = 5.0
                else:
                    new_coverage = float((line.split(" "))[1].strip()) / 2
        
                new_file.write(f"coverage: {new_coverage}")
                old_file.write(f"coverage: {new_coverage}")
            else:
                new_file.write(line)
                old_file.write(line)


    ref_file.close()
    new_file.close()
    old_file.close()

    ref_file = old_file

    return new_config_file_name


# Runs the NEAT read simulator using the given config file

def run_neat(config_file, output_dir, prefix):
    subprocess.run(["neat", "read-simulator", "-c", config_file, "-o", output_dir + "/" + prefix])


# General function for bacterial wrapper that calls all of the functions defined above

def bacterial_wrapper(reference_file, bacteria_name, ref_config_file):
    
    orig_seq = ""

    f = open(reference_file)
    fasta_header = f.readline().strip()

    for line in f:
        if line[0] != ">":
            orig_seq += line.strip()
        elif line.find("plasmid") != -1: # exclude plasmids from the sequence to be rearranged
            break
        
    f.close()

    output_dir = str(Path(ref_config_file).parent.absolute())

    rearranged_seq = wrapper(orig_seq)
    rearranged_seq_file = write_fasta_file(rearranged_seq, bacteria_name, fasta_header)
    new_config_file = write_config_file(ref_config_file, rearranged_seq_file, bacteria_name)

    run_neat(ref_config_file, output_dir, "Regular")
    run_neat(new_config_file, output_dir, "Wrapped")



# Read in bacterial reference genome files

salmonella_file = "BacterialWrapperScriptReferenceFiles/Salmonella/GCF_000006945.2_ASM694v2_genomic.fna"
pneumonia_file = "BacterialWrapperScriptReferenceFiles/Pneumonia/GCF_001457635.1_NCTC7465_genomic.fna"
tuberculosis_file = "BacterialWrapperScriptReferenceFiles/Tuberculosis/GCF_000195955.2_ASM19595v2_genomic.fna"



# Set config files for each reference file

salmonella_config_file = "Outputs/Salmonella/paired-ended/salmonella_config_test.yml"
pneumonia_config_file = "Outputs/Pneumonia/paired-ended/pneumonia_config_test.yml"
tuberculosis_config_file = "Outputs/Tuberculosis/paired-ended/tuberculosis_config_test.yml"

# Call the bacterial wrapper script on each reference file

# bacterial_wrapper(salmonella_file, "Salmonella", salmonella_config_file)
# bacterial_wrapper(pneumonia_file, "Pneumonia", pneumonia_config_file)
# bacterial_wrapper(tuberculosis_file, "Tuberculosis", tuberculosis_config_file)


# Stitching all outputs together - Keshav's script

def concat_fq(input_files: List[Path], dest: BgzfWriter) -> None:
    
    if not input_files:
        return
    
    for input_file in input_files:
        with bgzf.BgzfReader(input_file) as in_f:
            shutil.copyfileobj(in_f, dest)

def merge_bam(bams: List[Path], dest: Path, threads: int) -> None:
    
    if not bams:
        return

    unsorted = dest.with_suffix(".unsorted.bam")
    pysam.merge("--no-PG", "-@", str(threads), "-f", str(unsorted), *map(str, bams))
    pysam.sort("-@", str(threads), "-o", str(dest), str(unsorted))
    unsorted.unlink(missing_ok=True)

def merge_vcf(vcfs: List[Path], dest: Path) -> None:
    if not vcfs:
        return
    
    first, *rest = vcfs
    shutil.copy(first, dest)

    with dest.open("ab") as out_f:
        for vcf in rest:
            with vcf.open("rb") as fh:
                for line in fh:
                    if not line.startswith(b"#"):
                        out_f.write(line)

def stitch_all_outputs(files: List[Path]) -> None:
    fq1_list = []
    fq2_list = []
    vcf_list = []
    bam_list = []

    for file in files:
        file_name = file.stem # use stem to differentiate fq1 and fq2
        suffixes = file.suffixes # use suffixes to catch vcf and bam files
        
        if "r1.fastq" in file_name:
            fq1_list.append(file)
        elif "r2.fastq" in file_name:
            fq2_list.append(file)
        elif ".vcf" in suffixes:
            vcf_list.append(file)
        elif ".bam" in suffixes:
            bam_list.append(file)
    
    dest_fq1 = bgzf.BgzfWriter("stitched_fq1.bgzf")
    dest_fq2 = bgzf.BgzfWriter("stitched_fq2.bgzf")
    dest_bam = Path("stitched.bam")
    dest_vcf = Path("stitched.vcf")
    
    concat_fq(fq1_list, dest_fq1)
    concat_fq(fq2_list, dest_fq2)
    merge_bam(bam_list, dest_bam, 2)
    merge_vcf(vcf_list, dest_vcf)


# Stitching outputs for salmonella, single-ended reads

# orig_salmonella_fq = Path("Outputs/Salmonella/Regular.fastq.gz")
# orig_salmonella_bam = Path("Outputs/Salmonella/Regular_golden.bam")
# orig_salmonella_vcf = Path("Outputs/Salmonella/Regular_golden.vcf")

# new_salmonella_fq = Path("Outputs/Salmonella/Wrapped.fastq.gz")
# new_salmonella_bam = Path("Outputs/Salmonella/Wrapped_golden.bam")
# new_salmonella_vcf = Path("Outputs/Salmonella/Wrapped_golden.vcf")

# salmonella_inputs = [orig_salmonella_fq, orig_salmonella_bam, orig_salmonella_vcf, 
#                      new_salmonella_fq, new_salmonella_bam, new_salmonella_vcf]

# stitch_all_outputs(salmonella_inputs)


# Stitching outputs for salmonella, paired-ended reads

# orig_salmonella_fq1 = Path("Outputs/Salmonella/paired-ended/Regular_r1.fastq.gz")
# orig_salmonella_fq2 = Path("Outputs/Salmonella/paired-ended/Regular_r2.fastq.gz")
# orig_salmonella_bam = Path("Outputs/Salmonella/paired-ended/Regular_golden.bam")
# orig_salmonella_vcf = Path("Outputs/Salmonella/paired-ended/Regular_golden.vcf")

# new_salmonella_fq1 = Path("Outputs/Salmonella/paired-ended/Wrapped_r1.fastq.gz")
# new_salmonella_fq2 = Path("Outputs/Salmonella/paired-ended/Wrapped_r2.fastq.gz")
# new_salmonella_bam = Path("Outputs/Salmonella/paired-ended/Wrapped_golden.bam")
# new_salmonella_vcf = Path("Outputs/Salmonella/paired-ended/Wrapped_golden.vcf")

# salmonella_inputs = [orig_salmonella_fq1, orig_salmonella_fq2, orig_salmonella_bam, orig_salmonella_vcf, 
#                      new_salmonella_fq1, new_salmonella_fq2, new_salmonella_bam, new_salmonella_vcf]

# stitch_all_outputs(salmonella_inputs)


# Stitching outputs for pneumonia, single-ended reads

# orig_pneumonia_fq = Path("Outputs/Pneumonia/Regular.fastq.gz")
# orig_pneumonia_bam = Path("Outputs/Pneumonia/Regular_golden.bam")
# orig_pneumonia_vcf = Path("Outputs/Pneumonia/Regular_golden.vcf")

# new_pneumonia_fq = Path("Outputs/Pneumonia/Wrapped.fastq.gz")
# new_pneumonia_bam = Path("Outputs/Pneumonia/Wrapped_golden.bam")
# new_pneumonia_vcf = Path("Outputs/Pneumonia/Wrapped_golden.vcf")

# pneumonia_inputs = [orig_pneumonia_fq, orig_pneumonia_bam, orig_pneumonia_vcf, 
#                      new_pneumonia_fq, new_pneumonia_bam, new_pneumonia_vcf]

# stitch_all_outputs(pneumonia_inputs)


# Stitching outputs for pneumonia, paired-ended reads

# orig_pneumonia_fq1 = Path("Outputs/Pneumonia/paired-ended/Regular_r1.fastq.gz")
# orig_pneumonia_fq2 = Path("Outputs/Pneumonia/paired-ended/Regular_r2.fastq.gz")
# orig_pneumonia_bam = Path("Outputs/Pneumonia/paired-ended/Regular_golden.bam")
# orig_pneumonia_vcf = Path("Outputs/Pneumonia/paired-ended/Regular_golden.vcf")

# new_pneumonia_fq1 = Path("Outputs/Pneumonia/paired-ended/Wrapped_r1.fastq.gz")
# new_pneumonia_fq2 = Path("Outputs/Pneumonia/paired-ended/Wrapped_r2.fastq.gz")
# new_pneumonia_bam = Path("Outputs/Pneumonia/paired-ended/Wrapped_golden.bam")
# new_pneumonia_vcf = Path("Outputs/Pneumonia/paired-ended/Wrapped_golden.vcf")

# pneumonia_inputs = [orig_pneumonia_fq1, orig_pneumonia_fq2, orig_pneumonia_bam, orig_pneumonia_vcf, 
#                      new_pneumonia_fq1, new_pneumonia_fq2, new_pneumonia_bam, new_pneumonia_vcf]

# stitch_all_outputs(pneumonia_inputs)


# Stitching outputs for tuberculosis, single-ended reads

# orig_tuberculosis_fq = Path("Outputs/Tuberculosis/Regular.fastq.gz")
# orig_tuberculosis_bam = Path("Outputs/Tuberculosis/Regular_golden.bam")
# orig_tuberculosis_vcf = Path("Outputs/Tuberculosis/Regular_golden.vcf")

# new_tuberculosis_fq = Path("Outputs/Tuberculosis/Wrapped.fastq.gz")
# new_tuberculosis_bam = Path("Outputs/Tuberculosis/Wrapped_golden.bam")
# new_tuberculosis_vcf = Path("Outputs/Tuberculosis/Wrapped_golden.vcf")

# tuberculosis_inputs = [orig_tuberculosis_fq, orig_tuberculosis_bam, orig_tuberculosis_vcf, 
#                      new_tuberculosis_fq, new_tuberculosis_bam, new_tuberculosis_vcf]

# stitch_all_outputs(tuberculosis_inputs)


# Stitching outputs for tuberculosis, paired-ended readssource ~/anaconda3/bin/activate

# orig_tuberculosis_fq1 = Path("Outputs/Tuberculosis/paired-ended/Regular_r1.fastq.gz")
# orig_tuberculosis_fq2 = Path("Outputs/Tuberculosis/paired-ended/Regular_r2.fastq.gz")
# orig_tuberculosis_bam = Path("Outputs/Tuberculosis/paired-ended/Regular_golden.bam")
# orig_tuberculosis_vcf = Path("Outputs/Tuberculosis/paired-ended/Regular_golden.vcf")

# new_tuberculosis_fq1 = Path("Outputs/Tuberculosis/paired-ended/Wrapped_r1.fastq.gz")
# new_tuberculosis_fq2 = Path("Outputs/Tuberculosis/paired-ended/Wrapped_r2.fastq.gz")
# new_tuberculosis_bam = Path("Outputs/Tuberculosis/paired-ended/Wrapped_golden.bam")
# new_tuberculosis_vcf = Path("Outputs/Tuberculosis/paired-ended/Wrapped_golden.vcf")

# tuberculosis_inputs = [orig_tuberculosis_fq1, orig_tuberculosis_fq2, orig_tuberculosis_bam, orig_tuberculosis_vcf, 
#                      new_tuberculosis_fq1, new_tuberculosis_fq2, new_tuberculosis_bam, new_tuberculosis_vcf]

# stitch_all_outputs(tuberculosis_inputs)


# Testing functions

class TestWrapper(unittest.TestCase):
    def test_even(self):
        self.assertEqual(wrapper("ABBCBB"), "CBBABB")

    def test_odd(self):
        self.assertEqual(wrapper("ABBCBBC"), "BBCABBC")