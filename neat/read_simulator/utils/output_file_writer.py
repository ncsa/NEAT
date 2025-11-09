"""
Functions and classes for writing out the output.

TODO This file is in serious need of refactoring.
"""

__all__ = [
    "OutputFileWriter"
]

import os
import re
from struct import pack
import logging
from typing import Any

from Bio import bgzf
from pathlib import Path

#gzip for temp outs, bgzip for final outs
import gzip
from Bio.bgzf import BgzfWriter

from .read import Read
from .options import Options

_LOG = logging.getLogger(__name__)


# Some Constants
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
    :param vcf_format: The format to compress the vcf file, since we need to speed up the intermediate files with gzip
    """
    def __init__(self,
                 options: Options,
                 vcf_format: str = "bgzip",
                 bam_header: dict = None):

        self.paired_ended = options.paired_ended
        self.bam_header = bam_header
        self.vcf_format = vcf_format

        file_handles: dict[Path, Any] = {}

        # Set up filenames based on booleans
        if options.fq1 is not None:
            fq1 = options.fq1
            file_handles[fq1] = gzip.open(fq1, 'wt')
        else:
            fq1 = None
        if options.fq2 is not None:
            fq2 = options.fq2
            file_handles[fq2] = gzip.open(fq2, 'wt')
        else:
            fq2 = None
        if options.vcf is not None:
            vcf = options.vcf
            _open = gzip.open if self.vcf_format == "gzip" else bgzf.BgzfWriter
            file_handles[vcf] = _open(vcf, 'wt')
        else:
            vcf = None
        if options.bam is not None:
            bam = options.bam
            if bam_header:
                file_handles[bam] = bgzf.BgzfWriter(bam, 'w', compresslevel=6)
        else:
            bam = None

        self.fq1 = fq1
        self.fq2 = fq2
        self.vcf = vcf
        self.bam = bam

        if not file_handles:
            _LOG.error("output_file_writer received no files!")
            raise ValueError

        self.files_to_write = file_handles

        # Initialize the vcf and write the header, if applicable
        if options.produce_vcf:
            # Writing the vcf header.
            vcf_header = f'##fileformat=VCFv4.1\n' \
                         f'##reference={Path(options.reference).resolve()}\n' \
                         f'##source=NEAT\n' \
                         f'##RNG_seed={options.rng_seed}\n' \
                         f'##ALT=<ID=DEL,Description="Deletion">\n' \
                         f'##ALT=<ID=INS,Description="Insertion of novel sequence">\n' \
                         f'##ALT=<ID=SNP,Description="Single Nucleotide Polymorphism">\n' \
                         f'##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n' \
                         f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNEAT_simulated_sample\n'
            self.files_to_write[self.vcf].write(vcf_header)

        if options.produce_bam and bam_header:
            # bam header
            bam_handle = self.files_to_write[self.bam]
            bam_handle.write("BAM\1")
            # Without a header, we can't write these as bams.
            bam_header = "@HD\tVN:1.4\tSO:coordinate\n"
            for item in self.bam_header:
                bam_header += f'@SQ\tSN:{item}\tLN:{str(self.bam_header[item])}\n'
            bam_header += "@RG\tID:NEAT\tSM:NEAT\tLB:NEAT\tPL:NEAT\n"
            header_bytes = len(bam_header)
            num_refs = len(self.bam_header)
            bam_handle.write(pack('<i', header_bytes))
            bam_handle.write(bam_header)
            bam_handle.write(pack('<i', num_refs))
            # Contigs and lengths. If we can skip writing this out for intermediate files, great
            for item in self.bam_header:
                name_length = len(item) + 1
                bam_handle.write(pack('<i', name_length))
                bam_handle.write(f'{item}\0')
                bam_handle.write(pack('<i', self.bam_header[item]))

    def write_fastq_record(self, filename: Path, record: str):
        if filename in self.files_to_write:
            self.files_to_write[filename].write(record)
        else:
            _LOG.error(f"Tried to write fastq record to unknown file {filename}")
            raise ValueError

    def flush_and_close_files(self, skip_bam: bool = True):
        """
        flushes and closes all open files. On the thread level, we will want to go ahead and close the bam,
        but on the main level, we close the bam manually, so we don't need to do that here.

        :param skip_bam: Whether to skip trying to close the bam. This prevents a problem where the file was
            being "closed" twice and an extra line ending was getting written. Basically, the block level files
            need to be closed before they can be used by pysam.sort to make sure they are in the proper order.
        """
        for file_name in self.files_to_write:
            if file_name.suffix == "bam" and skip_bam:
                continue
            file_handle = self.files_to_write[file_name]
            try:
                file_handle.flush()
                os.fsync(file_handle.fileno())
                file_handle.close()
            except ValueError:
                _LOG.debug(f"file {file_name} already closed")

    def write_vcf_record(self, line: str):
        """
        This function takes in a list of temporary vcf files and combines them into a final output

        :param line: The format for a vcf line is:

            CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, SAMPLE
        """
        if self.vcf in self.files_to_write:
            self.files_to_write[self.vcf].write(line)
        else:
            _LOG.error(f"Tried to write to unknown vcf file {line}")
            raise ValueError

    def write_bam_record(
            self,
            read: Read,
            contig_id: int,
            bam_handle: BgzfWriter,
            read_length: int
    ):
        """
        Takes a read object and writes it out as a bam record

        :param read: Read object to write to file
        :param contig_id: the index of the reference for this
        :param bam_handle: the handle to write data to
        :param read_length: the length of the read to output
        """
        read_bin = reg2bin(read.position, read.end_point)

        mate_position = read.get_mpos()
        flag = read.calculate_flags(self.paired_ended)
        template_length = read.get_tlen()
        alt_sequence = read.read_sequence

        cigar = read.make_cigar()

        cig_letters = re.split(r"\d+", cigar)[1:]
        cig_numbers = [int(n) for n in re.findall(r"\d+", cigar)]
        cig_ops = len(cig_letters)

        next_ref_id = contig_id

        if not mate_position:
            next_pos = 0
            template_length = 0
        else:
            next_pos = mate_position

        encoded_cig = bytearray()

        for i in range(cig_ops):
            encoded_cig.extend(pack('<I', (cig_numbers[i] << 4) + CIGAR_PACKED[cig_letters[i]]))

        encoded_seq = bytearray()
        encoded_len = (read_length + 1) // 2
        seq_len = read_length
        if seq_len & 1:
            alt_sequence += '='
        for i in range(encoded_len):
            # if self.debug:
            #     # Note: trying to remove all this part
            encoded_seq.extend(
                pack('<B',
                     (SEQ_PACKED[alt_sequence[2 * i].capitalize()] << 4) +
                     SEQ_PACKED[alt_sequence[2 * i + 1].capitalize()]))

        # In NEAT 2.0, this was `encodedQual = ''.join([chr(ord(n)-33) for n in qual])` but this converts the char back into
        # the original quality score, which we saved, so we'll try just using that.
        encoded_qual = "".join([chr(n) for n in read.quality_array])

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

        block_size = 32 + len(read.name) + 1 + len(encoded_cig) + len(encoded_seq) + len(encoded_qual)

        bam_handle.write((
                pack('<i', block_size) +
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
                encoded_qual.encode('utf-8')
        ))