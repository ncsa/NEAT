"""
Functions and classes for writing out the output.
"""

__all__ = [
    "OutputFileWriter"
]

import gzip
import re
import shutil
from struct import pack
import logging
import pysam
import numpy as np
import pickle

from Bio import bgzf
from Bio import SeqIO
from Bio.Seq import Seq, MutableSeq
from typing import Iterator, TextIO
from pathlib import Path
from numpy.random import Generator

from ...common import validate_output_path, open_output, open_input, ALLOWED_NUCL
from .read import Read
from .options import Options
from .neat_cigar import CigarString

_LOG = logging.getLogger(__name__)


# Some Constants
# TODO make bam compression a configurable option
BAM_COMPRESSION_LEVEL = 6
CIGAR_PACKED = {'M': 0, 'I': 1, 'D': 2, 'N': 3, 'S': 4, 'H': 5, 'P': 6, '=': 7, 'X': 8}
SEQ_PACKED = {'=': 0, 'A': 1, 'C': 2, 'M': 3, 'G': 4, 'R': 5, 'S': 6, 'V': 7,
              'T': 8, 'W': 9, 'Y': 10, 'H': 11, 'K': 12, 'D': 13, 'B': 14, 'N': 15}
# TODO figure out an optimum batch size or get rid of this idea
BUFFER_BATCH_SIZE = 8000  # write out to file after this many reads


def reg2bin(beg: int, end: int):
    """
    Samtools reg2bin function.

    Finds the largest superset bin of region. Numeric values taken from hts-specs
    Note: description of this function taken from source code for bamnostic.bai
        (https://bamnostic.readthedocs.io/en/latest/_modules/bamnostic/bai.html)
    :param beg: inclusive beginning position of region
    :param end: exclusive end position of region
    :return: distinct bin ID or largest superset bin of region
    """
    end -= 1
    if beg >> 14 == end >> 14:
        return ((1 << 15) - 1) // 7 + (beg >> 14)
    if beg >> 17 == end >> 17:
        return ((1 << 12) - 1) // 7 + (beg >> 17)
    if beg >> 20 == end >> 20:
        return ((1 << 9) - 1) // 7 + (beg >> 20)
    if beg >> 23 == end >> 23:
        return ((1 << 6) - 1) // 7 + (beg >> 23)
    if beg >> 26 == end >> 26:
        return ((1 << 3) - 1) // 7 + (beg >> 26)
    return 0


class OutputFileWriter:
    """
    This class sets up the output files and has methods for writing out records
    in the various formats.

    :param options: Options for the current run.
    :param bam_header: A dictionary of lengths of each contig from the reference, keyed by contig id.
    """
    def __init__(self,
                 options: Options,
                 bam_header: dict = None):

        # set the booleans
        self.write_fastq = options.produce_fastq
        self.write_fasta = options.produce_fasta
        self.write_bam = options.produce_bam
        self.write_vcf = options.produce_vcf
        self.paired = options.paired_ended
        self.temporary_dir = options.temp_dir_path

        self.bam_header = bam_header

        # Set the file names
        self.fasta_fns = None

        self.fastq_fns = None
        self.fastq1_fn = None
        self.fastq2_fn = None

        self.bam_fn = None
        self.vcf_fn = None

        # Set up filenames based on booleans
        files_to_write = []
        if self.write_fasta:
            if options.ploidy > 1 and options.fasta_per_ploid:
                self.fasta_fns = [options.output.parent / f'{options.output.stem}_ploid{i+1}.fasta.gz'
                                  for i in range(options.ploidy)]
            else:
                self.fasta_fns = [options.output.parent / f'{options.output.stem}.fasta.gz']
            files_to_write.extend(self.fasta_fns)
        if self.paired and self.write_fastq:
            self.fastq1_fn = options.output.parent / f'{options.output.stem}_r1.fastq.gz'
            self.fastq2_fn = options.output.parent / f'{options.output.stem}_r2.fastq.gz'
            self.fastq_fns = (self.fastq1_fn, self.fastq2_fn)
            files_to_write.extend(self.fastq_fns)
        elif self.write_fastq:
            self.fastq1_fn = options.output.parent / f'{options.output.stem}.fastq.gz'
            self.fastq_fns = self.fastq1_fn
            files_to_write.append(self.fastq_fns)
        if self.write_bam:
            self.bam_fn = options.output.parent / f'{options.output.stem}_golden.bam'
            self.bam_keys = list(bam_header)
            files_to_write.append(self.bam_fn)
        if self.write_vcf:
            self.vcf_fn = options.output.parent / f'{options.output.stem}_golden.vcf.gz'
            files_to_write.append(self.vcf_fn)

        self.files_to_write = files_to_write

        # Create files as applicable
        for file in self.files_to_write:
            validate_output_path(file, True, options.overwrite_output)

        mode = 'xt'
        if options.overwrite_output:
            mode = 'wt'
        # Initialize the vcf and write the header, if applicable
        if self.write_vcf:
            # Writing the vcf header.
            with open_output(self.vcf_fn, mode=mode) as vcf_file:
                vcf_file.write(f'##fileformat=VCFv4.1\n')
                vcf_file.write(f'##reference={Path(options.reference).resolve()}\n')
                vcf_file.write(f'##Generated by NEAT with RNG value: {options.rng_seed}\n')
                vcf_file.write(f'##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">\n')
                vcf_file.write(f'##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n')
                vcf_file.write(f'##INFO=<ID=VMX,Number=1,Type=String,'
                               f'Description="SNP is Missense in these Read Frames">\n')
                vcf_file.write(f'##INFO=<ID=VNX,Number=1,Type=String,'
                               f'Description="SNP is Nonsense in these Read Frames">\n')
                vcf_file.write(f'##INFO=<ID=VFX,Number=1,Type=String,Description="Indel Causes Frameshift">\n')
                vcf_file.write(f'##ALT=<ID=DEL,Description="Deletion">\n')
                vcf_file.write(f'##ALT=<ID=DUP,Description="Duplication">\n')
                vcf_file.write(f'##ALT=<ID=INS,Description="Insertion of novel sequence">\n')
                vcf_file.write(f'##ALT=<ID=INV,Description="Inversion">\n')
                vcf_file.write(f'##ALT=<ID=CNV,Description="Copy number variable region">\n')
                vcf_file.write(f'##ALT=<ID=TRANS,Description="Translocation">\n')
                vcf_file.write(f'##ALT=<ID=INV-TRANS,Description="Inverted translocation">\n')
                vcf_file.write(f'##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
                # Add a neat sample column
                vcf_file.write(f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNEAT_simulated_sample\n')

    def merge_temp_vcfs(self, temporary_files: list):
        """
        This function takes in a list of temporary vcf files and combines them into a final output

        :param temporary_files: The list of temporary files to combine
        """
        with open_output(self.vcf_fn) as vcf_out:
            for temp_file in temporary_files:
                with open_input(temp_file) as infile:
                    vcf_out.write(infile.read())

    def merge_temp_fastas(self, temporary_files: list):
        """
        Takes a list of temporary fasta files and puts them into a final file

        :param temporary_files: A list of temporary fastas to combine
        """
        for file in self.fasta_fns:
            with open_output(file) as vcf_out:
                for temp_file in temporary_files:
                    with open_input(temp_file) as infile:
                        vcf_out.write(infile.read())

    def merge_temp_fastqs(
            self, fastq_files: list, rand_num_gen: Generator
    ):
        """
        Takes a list of fastqs and combines them into a final output. This is the most complicated one, because we need
        to randomize the fastq to make it more realistic. If we are not producing a fastq, we still need to perform
        this step to get the final names of the reads for the sam file. In that case, it should go quickly since there
        will be no IO time.

        :param fastq_files: The temporary fastq files to combine in the final output, or to set the order for the
            final bam
        :param rand_num_gen: the random number generator for the run
        """
        fastq_indexes = []
        read1_keys = []
        read2_keys = []

        # Index the temp fastqs
        for file_pair in fastq_files:
            file1_index = SeqIO.index(str(file_pair[0]), 'fastq')
            file2_index = None
            if file_pair[1]:
                file2_index = SeqIO.index(str(file_pair[1]), 'fastq')
                read2_keys.append(list(file2_index))
            else:
                read2_keys.append(file2_index)
            fastq_indexes.append((file1_index, file2_index))
            read1_keys.append(list(file1_index))

        # Flatten keys lists to shuffle the order
        all_read1_keys = [key for sublist in read1_keys for key in sublist]
        all_read2_keys = [key for sublist in read2_keys for key in sublist]

        # Shuffle the keys
        rand_num_gen.shuffle(all_read1_keys)

        read2_singletons = [key for key in all_read2_keys if key not in all_read1_keys]

        with (
            open_output(self.fastq1_fn) as fq1,
            open_output(self.fastq2_fn) as fq2
        ):
            for i in range(len(all_read1_keys)):

                current_key = all_read1_keys[i]

                # This gives us the first index of the 'row' in the keys table where this read is located,
                #  giving us the list index of the file index where I can find this read
                current_index = [read1_keys.index(x) for x in read1_keys if current_key in x][0]
                # This is the index of the read1 file, at the current key, which should be the read we want
                read1 = fastq_indexes[current_index][0][current_key]
                SeqIO.write(read1, fq1, 'fastq')

                if current_key in fastq_indexes[current_index][1]:
                    read2 = fastq_indexes[current_index][1][current_key]
                    SeqIO.write(read2, fq2, 'fastq')

            for j in range(len(read2_singletons)):

                current_key = read2_singletons[j]

                current_index = [read2_keys.index(x) for x in read2_keys if current_key in x][0]
                read = fastq_indexes[current_index][1][current_key]
                SeqIO.write(read, fq2, 'fastq')

    def output_bam_file(self, reads_files: list, contig_dict: dict):
        """
        This section is for producing a CIGAR string using a temp sam file(sam file with
        original sequence instead of a cigar string)

        :param reads_files: The list of temp sams to combine
        :param contig_dict: A dictionary with the keys as contigs from the reference,
            and the values the index of that contig
        """

        bam_out = bgzf.BgzfWriter(self.bam_fn, 'w', compresslevel=BAM_COMPRESSION_LEVEL)
        bam_out.write("BAM\1")
        header = "@HD\tVN:1.4\tSO:coordinate\n"
        for item in self.bam_header:
            header += f'@SQ\tSN:{item}\tLN:{str(self.bam_header[item])}\n'
        header += "@RG\tID:NEAT\tSM:NEAT\tLB:NEAT\tPL:NEAT\n"
        header_bytes = len(header)
        num_refs = len(self.bam_header)
        bam_out.write(pack('<i', header_bytes))
        bam_out.write(header)
        bam_out.write(pack('<i', num_refs))

        for item in self.bam_header:
            name_length = len(item) + 1
            bam_out.write(pack('<i', name_length))
            bam_out.write(f'{item}\0')
            bam_out.write(pack('<i', self.bam_header[item]))

        for file in reads_files:
            contig_reads_data = pickle.load(gzip.open(file))
            for key in contig_reads_data:
                read1 = contig_reads_data[key][0]
                read2 = contig_reads_data[key][0]
                if read1:
                    self.write_bam_record(read1, contig_dict[read1.reference_id], bam_out)

                if read2:
                    self.write_bam_record(read2, contig_dict[read2.reference_id], bam_out)
        bam_out.close()

    def write_bam_record(self, read: Read, contig_id: int, bam_handle: bgzf.BgzfWriter):
        """
        Takes a read object and writes it out as a bam record

        :param read: a read object containing everything we need to write it out.
        :param contig_id: the index of the reference for this
        :param bam_handle: the handle of the file object to write to.
        """
        read_bin = reg2bin(read.position, read.end_point)

        mate_position = read.get_mpos()
        flag = read.calculate_flags(self.paired)
        template_length = read.get_tlen()
        alt_sequence = read.read_sequence

        cigar = read.make_cigar()

        cig_letters = re.split(r"\d+", cigar)[1:]
        cig_numbers = [int(n) for n in re.findall(r"\d+", cigar)]
        cig_ops = len(cig_letters)

        next_ref_id = contig_id

        if mate_position is None:
            next_pos = 0
            template_length = 0
        else:
            next_pos = mate_position

        encoded_cig = bytearray()
        for i in range(cig_ops):
            encoded_cig.extend(pack('<I', (cig_numbers[i] << 4) + CIGAR_PACKED[cig_letters[i]]))
        encoded_seq = bytearray()
        encoded_len = (len(alt_sequence) + 1) // 2
        seq_len = len(alt_sequence)
        if seq_len & 1:
            alt_sequence += '='
        for i in range(encoded_len):
            # if self.debug:
            #     # Note: trying to remove all this part
            encoded_seq.extend(
                pack('<B',
                     (SEQ_PACKED[alt_sequence[2 * i].capitalize()] << 4) +
                     SEQ_PACKED[alt_sequence[2 * i + 1].capitalize()]))

        # apparently samtools automatically adds 33 to the quality score string...
        encoded_qual = ''.join([chr(ord(n) - 33) for n in read.read_quality_score])

        """
        block_size = 4 +		# refID 		int32
                     4 +		# pos			int32
                     4 +		# bin_mq_nl		uint32
                     4 +		# flag_nc		uint32
                     4 +		# l_seq			int32
                     4 +		# next_ref_id	int32
                     4 +		# next_pos		int32
                     4 +		# tlen			int32
                     len(readName)+1 +
                     len(encoded cigar) +
                     encoded_len +
                     len(seq)
        """

        block_size = 32 + len(read.name) + 1 + len(encoded_cig) + len(encoded_seq) + len(read.read_quality_score)

        bam_handle.write((pack('<i', block_size) +
                          pack('<i', contig_id) +
                          pack('<i', read.position + 1) +
                          pack('<I', (read_bin << 16)
                               + (read.mapping_quality << 8)
                               + len(read.name)
                               + 1) +
                          pack('<I', (flag << 16) + cig_ops) +
                          pack('<i', seq_len) +
                          pack('<i', next_ref_id) +
                          pack('<i', next_pos) +
                          pack('<i', template_length) +
                          read.name.encode('utf-8') + b'\0' +
                          encoded_cig +
                          encoded_seq +
                          encoded_qual.encode('utf-8')))
