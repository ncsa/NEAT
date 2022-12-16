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

from Bio import bgzf
from Bio import SeqIO
from pathlib import Path
from numpy.random import Generator

from ...common import validate_output_path, open_output, open_input
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
BUFFER_BATCH_SIZE = 8000  # write out to file_list after this many reads


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


# takes list of strings, returns numerical flag
def sam_flag(string_list: list) -> int:
    out_val = 0
    strings_to_check = []
    for val in string_list:
        if val not in strings_to_check:
            strings_to_check.append(val)
    for n in string_list:
        if n == 'paired':
            out_val += 1
        elif n == 'proper':
            out_val += 2
        elif n == 'unmapped':
            out_val += 4
        elif n == 'mate_unmapped':
            out_val += 8
        elif n == 'reverse':
            out_val += 16
        elif n == 'mate_reverse':
            out_val += 32
        elif n == 'first':
            out_val += 64
        elif n == 'second':
            out_val += 128
        elif n == 'not_primary':
            out_val += 256
        elif n == 'low_quality':
            out_val += 512
        elif n == 'duplicate':
            out_val += 1024
        elif n == 'supplementary':
            out_val += 2048
    return out_val


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

        # Set the file_list names
        self.fasta_fns = None

        self.fastq_fns = None
        self.fastq1_fn = None
        self.fastq2_fn = None

        self.bam_fn = None
        self.vcf_fn = None
        self.sam_fn = None

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
            self.sam_fn = self.temporary_dir / f'{options.output.stem}_temp.sam'
            self.bam_keys = list(bam_header.keys())
            files_to_write.append(self.sam_fn)
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

        if self.write_bam:
            # Write the temp sam file_list
            with open_output(self.sam_fn, mode=mode) as s:
                s.write('@HD\tVN:1.5\tSO:coordinate\n')
                for key in self.bam_keys:
                    s.write(f'@SQ\tSN:{str(key)}\tLN:{str(self.bam_header[key])}\n')
                s.write('@RG\tID:NEAT\tSM:NEAT\tLB:NEAT\tPL:NEAT\n')

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
        Takes a list of temporary fasta files and puts them into a final file_list

        :param temporary_files: A list of temporary fastas to combine
        """
        for file in self.fasta_fns:
            with open_output(file) as vcf_out:
                for temp_file in temporary_files:
                    with open_input(temp_file) as infile:
                        vcf_out.write(infile.read())

    def merge_temp_fastqs_and_sort_bams(
            self, fastq_files: list, produce_fastq: bool, sam_file_order: list, rand_num_gen: Generator
    ):
        """
        Takes a list of fastqs and combines them into a final output. This is the most complicated one, because we need
        to randomize the fastq to make it more realistic. If we are not producing a fastq, we still need to perform
        this step to get the final names of the reads for the sam file_list. In that case, it should go quickly since there
        will be no IO time.

        :param fastq_files: The temporary fastq files to combine in the the final output, or to set the order for the
            final bam
        :param produce_fastq: Whether to actually write the fastq files
        :param sam_file_order: The order of the sam file_list records
        :param rand_num_gen: the random number generator for the run
        """
        fastq_indexes = []
        read1_keys = []
        read2_keys = []

        # Step 1 randomize the final fastq and rename the records in order added to the final
        for file_pair in fastq_files:
            # TODO: use biopython seqIO for fastqs. Randomize the order based on the number of records we have, then
            #   output them in a shuffled order
            file1_index = SeqIO.index(str(file_pair[0]), 'fastq')
            file2_index = None
            if file_pair[1]:
                file2_index = SeqIO.index(str(file_pair[1]), 'fastq')
                read2_keys.append(list(file2_index.keys()))
            else:
                read2_keys.append(file2_index)
            fastq_indexes.append((file1_index, file2_index))
            read1_keys.append(list(file1_index.keys()))

        all_read1_keys = [key for sublist in read1_keys for key in sublist]
        # Shuffle the keys
        rand_num_gen.shuffle(all_read1_keys)
        all_read2_keys = [key for sublist in read2_keys for key in sublist]
        read2_singletons = [key for key in all_read2_keys if key not in all_read1_keys]

        flattened_sam_file_order = [x for sublist in sam_file_order for x in sublist]
        tsam_remap = {x: "" for x in flattened_sam_file_order}
        base_name = "neat-generated-read:"
        with (
            open_output(self.fastq1_fn) as fq1,
            open_output(self.fastq2_fn) as fq2
        ):
            for i in range(len(all_read1_keys)):
                current_name = f'{base_name} {i}'
                current_key = all_read1_keys[i]
                tsam_remap[current_key] = current_name
                if produce_fastq:
                    # This gives us the first index of the 'row' in the keys table where this read is located,
                    #  giving us the list index of the file_list index where I can find this read
                    current_index = [read1_keys.index(x) for x in read1_keys if current_key in x][0]
                    # This is the index of the read1 file_list, at the current key, which should be the read we want
                    read1 = fastq_indexes[current_index][0][current_key]
                    read1.id = f'{current_name}/1'
                    SeqIO.write(read1, fq1, 'fastq')
                if current_key in fastq_indexes[current_index][1] and produce_fastq:
                    read2 = fastq_indexes[current_index][1][current_key]
                    read2.id = f'{current_name}/2'
                    SeqIO.write(read2, fq2, 'fastq')

            for j in range(len(read2_singletons)):
                # This names them starting at the last name we used
                current_name = f'{base_name} {j + len(all_read1_keys)}'
                current_key = read2_singletons[j]
                tsam_remap[current_key] = current_name
                if produce_fastq:
                    current_index = [read2_keys.index(x) for x in read2_keys if current_key in x][0]
                    read = fastq_indexes[current_index][1][current_key]
                    read.id = f'{current_name}/2'
                    SeqIO.write(read, fq2, 'fastq')

        return tsam_remap

    def combine_tsam_into_bam(self, temp_files: list, renaming_dictionary: dict, sam_order: list):
        _LOG.info("Bam creation coming soon!")
        pass


    def write_sam_record(self, chromosome_index, read_name, pos_0, cigar, seq, qual, output_sam_flag, rnext="=",
                         mate_pos=None, aln_map_quality: int = 70):
        """
        okay, this might be a tricky bit because we have to keep track of when the last record for the chromosome
        has been reached. Once it has, we have to update the RNEXT field to something else.  Have to think on this.

        I'll look more closely at how Zach originally handled this.
        I don't see any real reason this can't be calculated in the moin body of the code then just input into
        this function. Outside will have info about the chromosome and stuff.
        """

        my_t_len = 0
        if mate_pos:
            if mate_pos > pos_0:
                my_t_len = mate_pos - pos_0 + len(seq)
            else:
                my_t_len = mate_pos - pos_0 - len(seq)

        record = f'{read_name}\t{output_sam_flag}\t{chromosome_index}\t{pos_0 + 1}\t' \
                 f'{aln_map_quality}\t{cigar}\t{rnext}\t{mate_pos}\t{my_t_len}\t{seq}\t{qual}\n'

        with open(self.sam_fn, 'a') as f:
            f.write(record)

    def write_bam_record(self, chromosome_index, read_name, pos_0, cigar, seq, qual, output_sam_flag,
                         mate_pos=None, aln_map_quality: int = 70):
        """
        We're going to completely change how this works. Instead of trying to keep track of the cigar during
        the simulation, we'll do a post-hoc effort and determine the bam that way. It will be less perfect,
        but will be much faster, and the ways it will differ from reality will not matter to the alignment.

        So we need to write the bam stuff from the fastq to a temporary file_list. That either needs to be handled here
        or in the main code.
        """
        if self.bam_first_write:
            # First time we open the file_list, let's write it. This is due to a limitation in bgzf where it cannot append
            # items like you can with open.
            self.bam_first_write = False
            self.bam_file = bgzf.BgzfWriter(self.bam_fn, 'w', compresslevel=BAM_COMPRESSION_LEVEL)
            self.bam_file.write("BAM\1")
            self.bam_file.write(pack('<i', self.bamfile_header_bytes))
            self.bam_file.write(self.bamfile_header)
            self.bam_file.write(pack('<i', self.bamfile_num_refs))

            for key in self.bam_keys:
                l_name = len(key) + 1
                self.bam_file.write(pack('<i', l_name))
                self.bam_file.write(key + '\0')
                self.bam_file.write(pack('<i', len(self.bam_header[key])))

        my_bin = reg2bin(pos_0, pos_0 + len(seq))
        # my_bin     = 0	# or just use a dummy value, does this actually matter?

        my_map_quality = aln_map_quality
        cigar_string = CigarString.list_to_string(cigar)
        cig_letters = re.split(r"\d+", cigar_string)[1:]
        cig_numbers = [int(n) for n in re.findall(r"\d+", cigar_string)]
        cig_ops = len(cig_letters)
        next_ref_id = chromosome_index
        if mate_pos is None:
            next_pos = 0
            my_t_len = 0
        else:
            next_pos = mate_pos
            if next_pos > pos_0:
                my_t_len = next_pos - pos_0 + len(seq)
            else:
                my_t_len = next_pos - pos_0 - len(seq)

        encoded_cig = bytearray()
        for i in range(cig_ops):
            encoded_cig.extend(pack('<I', (cig_numbers[i] << 4) + CIGAR_PACKED[cig_letters[i]]))
        encoded_seq = bytearray()
        encoded_len = (len(seq) + 1) // 2
        seq_len = len(seq)
        if seq_len & 1:
            seq += '='
        for i in range(encoded_len):
            if self.debug:
                # Note: trying to remove all this part
                _LOG.debug(f'{seq[2 * i]}, {seq[2 * i + 1]}')
            encoded_seq.extend(
                pack('<B', (SEQ_PACKED[seq[2 * i].capitalize()] << 4) + SEQ_PACKED[seq[2 * i + 1].capitalize()]))

        # apparently samtools automatically adds 33 to the quality score string...
        encoded_qual = ''.join([chr(ord(n) - 33) for n in qual])

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
                     4*cig_ops +
                     encoded_len +
                     len(seq)
        """
        # block_size = 32 + len(readName)+1 + 4*cig_ops + encoded_len + len(seq)
        block_size = 32 + len(read_name) + 1 + len(encoded_cig) + len(encoded_seq) + len(encoded_qual)

        """
        Not sure what the point of the following lines are
        # self.bam_file.write(pack('<i',block_size))
        # self.bam_file.write(pack('<i',refID))
        # self.bam_file.write(pack('<i',pos_0))
        # self.bam_file.write(pack('<I',(my_bin<<16) + (my_map_quality<<8) + len(readName)+1))
        # self.bam_file.write(pack('<I',(samFlag<<16) + cig_ops))
        # self.bam_file.write(pack('<i',seq_len))
        # self.bam_file.write(pack('<i',next_ref_id))
        # self.bam_file.write(pack('<i',next_pos))
        # self.bam_file.write(pack('<i',my_tlen))
        # self.bam_file.write(readName+'\0')
        # self.bam_file.write(encoded_cig)
        # self.bam_file.write(encoded_seq)
        # self.bam_file.write(encoded_qual)
        """

        # a horribly compressed line, I'm sorry.
        # (ref_index, position, data)
        self.bam_file.write((pack('<i', block_size) + pack('<i', chromosome_index) + pack('<i', pos_0) +
                             pack('<I', (my_bin << 16) + (my_map_quality << 8) + len(read_name) + 1) +
                             pack('<I', (output_sam_flag << 16) + cig_ops) + pack('<i', seq_len) +
                             pack('<i', next_ref_id) +
                             pack('<i', next_pos) + pack('<i', my_t_len) + read_name.encode('utf-8') +
                             b'\0' + encoded_cig + encoded_seq + encoded_qual.encode('utf-8')))

    def sort_and_convert_sam_file(self):
        """
        Convert sam to bam
        """
        samout_fn = self.bam_fn.parent / (self.bam_fn.stem + "_fromsam" + self.bam_fn.suffix)
        if self.debug:
            shutil.copyfile(self.sam_fn, '/home/joshfactorial/Documents/neat_outputs/test.sam')
        pysam.sort("-o", str(samout_fn), str(self.sam_fn))
        self.sam_temp_dir.cleanup()

    def sort_bam(self):
        if self.bam_file:
            self.bam_file.close()
        output = self.temporary_dir / 'sorted.bam'
        pysam.sort("-o", str(output), str(self.bam_fn))
        shutil.copyfile(output, self.bam_fn)

    def close_bam_file(self):
        if self.bam_file:
            self.bam_file.close()
