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

def write_fasta_file(new_seq, bacteria_name, fasta_header, output_dir_path, type):
    fasta_file_name = f"{type}_{bacteria_name}.fna"
    fasta_file_path = output_dir_path / fasta_file_name
    fasta_file = open(fasta_file_path, "w")

    fasta_file.write(fasta_header + "\n" + new_seq)

    fasta_file.close()

    return fasta_file_path


# Writes a yml configuration file for the newly rearranged chromosome's fasta sequence
# Splits the coverage in half for the reference and new config files
# These use default values for all other parameters for NEAT

def write_config_file(ref_config_file, rearranged_seq_file, orig_seq_file, bacteria_name, output_dir_path):
    new_config_file_name = f"new_{bacteria_name}_config_test.yml"
    old_config_file_name = f"{bacteria_name}_config_test.yml"

    new_config_file_path = output_dir_path / new_config_file_name
    old_config_file_path = output_dir_path / old_config_file_name

    with open(ref_config_file, 'r') as ref_file, open(new_config_file_path, 'w') as new_file, open(old_config_file_path, 'w') as old_file:
        for line in ref_file:
            if line.find("reference:") != -1:
                new_file.write(f"reference: {rearranged_seq_file}\n")
                old_file.write(f"reference: {orig_seq_file}\n")
            elif line.find("coverage:") != -1:
                if line.strip() == "coverage: .":
                    new_coverage = 5.0
                else:
                    new_coverage = float((line.split(" "))[1].strip()) // 2
        
                new_file.write(f"coverage: {new_coverage}\n")
                old_file.write(f"coverage: {new_coverage}\n")
            else:
                new_file.write(line)
                old_file.write(line)


    ref_file.close()
    new_file.close()
    old_file.close()

    return old_config_file_path, new_config_file_path


# Runs the NEAT read simulator using the given config file

def run_neat(config_file, output_dir, prefix):
    subprocess.run(["neat", "read-simulator", "-c", config_file, "-o", output_dir + "/" + prefix])


# General function for bacterial wrapper that calls all of the functions defined above

def bacterial_wrapper(reference_file, bacteria_name, ref_config_file, output_dir):
    
    orig_seq = ""

    f = open(reference_file)
    fasta_header = f.readline().strip()

    plasmids = False

    for line in f:
        if line[0] != ">":
            orig_seq += line.strip()
        elif line.lower().find("plasmid") != -1:  # exclude plasmids from the sequence to be rearranged
            plasmids = True
            break
        
    f.close()

    output_dir_path = Path(output_dir)

    rearranged_seq = wrapper(orig_seq)
    rearranged_seq_file = write_fasta_file(rearranged_seq, bacteria_name, fasta_header, output_dir_path, "wrapped")

    orig_seq_file = reference_file
    if plasmids:
        orig_seq_file = write_fasta_file(orig_seq, bacteria_name, fasta_header, output_dir_path, "orig")
    
    config_files = write_config_file(ref_config_file, rearranged_seq_file, orig_seq_file, bacteria_name, output_dir_path)
    old_config_file = config_files[0]
    new_config_file = config_files[1]

    run_neat(old_config_file, output_dir, "Regular")
    run_neat(new_config_file, output_dir, "Wrapped")


# Stitching all outputs together - Keshav's script

def concat_fq(input_files: List[Path], dest: Path) -> None:
    
    if not input_files:
        # Nothing to do, and no error to throw
        return

    with gzip.open(dest, 'wt') as out_f:
        for input_file in input_files:
            with gzip.open(input_file, 'rt') as in_f:
                shutil.copyfileobj(in_f, out_f)

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

def stitch_all_outputs(files: List[Path], output_dir) -> None:
    fq1_list = []
    fq2_list = []
    vcf_list = []
    bam_list = []

    for file in files:
        file_name = file.stem  # use stem to differentiate fq1 and fq2
        suffixes = file.suffixes  # use suffixes to catch vcf and bam files

        if "r2.fastq" in file_name:
            fq2_list.append(file)
        elif "r1.fastq" in file_name or ".fastq" in suffixes:
            fq1_list.append(file)
        elif ".vcf" in suffixes and ".tbi" not in suffixes:
            vcf_list.append(file)
        elif ".bam" in suffixes and ".bai" not in suffixes:
            bam_list.append(file)

    dest_fq1 = Path(f"{output_dir}/stitched_fq1.gz")
    dest_bam = Path(f"{output_dir}/stitched.bam")
    dest_vcf = Path(f"{output_dir}/stitched.vcf")
    
    concat_fq(fq1_list, dest_fq1)
    
    if (fq2_list):
        dest_fq2 = Path(f"{output_dir}/stitched_fq2.gz")
        concat_fq(fq2_list, dest_fq2)

    merge_bam(bam_list, dest_bam, 2)
    merge_vcf(vcf_list, dest_vcf)


# Testing functions

class TestWrapper(unittest.TestCase):
    def test_even(self):
        self.assertEqual(wrapper("ABBCBB"), "CBBABB")

    def test_odd(self):
        self.assertEqual(wrapper("ABBCBBC"), "BBCABBC")