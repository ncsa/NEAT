"""
Functions and classes for writing out the output.

TODO This file is in serious need of refactoring.
"""

__all__ = [
    "OutputFileWriter"
]

import os
import shutil
import time
from struct import pack
import logging
from typing import Any

import numpy as np
from Bio import bgzf
from Bio import SeqIO
from pathlib import Path
from numpy.random import Generator

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

# 256-entry byte→4-bit lookup tables: write_bam_record encodes the read sequence
# 2 bases per byte and packs each CIGAR op into one uint32. A flat numpy array
# indexed by byte value is ~30x faster than a per-character dict lookup in a Python
# loop, which dominated write_bam_record self-time before vectorization.
SEQ_PACKED_LOOKUP = np.zeros(256, dtype=np.uint8)
for _c, _v in SEQ_PACKED.items():
    SEQ_PACKED_LOOKUP[ord(_c)] = _v
CIGAR_PACKED_LOOKUP = np.zeros(256, dtype=np.uint8)
for _c, _v in CIGAR_PACKED.items():
    CIGAR_PACKED_LOOKUP[ord(_c)] = _v

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
    :param vcf_header: optional parameters for vcf header.
    :param bam_header: A dictionary of lengths of each contig from the reference, keyed by contig id.
    :param vcf_format: The format to compress the vcf file, since we need to speed up the intermediate files with gzip
    """
    def __init__(self,
                 options: Options,
                 vcf_header: dict = None,
                 vcf_format: str = "bgzip",
                 bam_header: dict = None):

        self.paired_ended = options.paired_ended
        self.vcf_header = vcf_header
        self.bam_header = bam_header
        self.vcf_format = vcf_format
        self.tmp_dir = options.temp_dir_path

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

        if not file_handles and self.bam is None:
            _LOG.error("output_file_writer received no files!")
            raise ValueError

        self.files_to_write = file_handles

        # Initialize the vcf and write the header, if applicable
        if options.produce_vcf and self.vcf_header:
            # Writing the vcf header.
            ref_line = ""
            for (contig, length) in self.vcf_header.items():
                ref_line += f"##contig=<ID={contig}, length={length}>\n"
            vcf_header = f'##fileformat=VCFv4.1\n' \
                         f'##reference={Path(options.reference).resolve()}\n' \
                         f'##source=NEAT\n' \
                         f'##RNG_seed={options.rng_seed}\n' \
                         f'##ALT=<ID=DEL,Description="Deletion">\n' \
                         f'##ALT=<ID=INS,Description="Insertion of novel sequence">\n' \
                         f'##ALT=<ID=SNP,Description="Single Nucleotide Polymorphism">\n' \
                         f'##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n' \
                         f'{ref_line}' \
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
        if read.is_reverse:
            # having weird issues with the bam
            alt_sequence = read.read_sequence.reverse_complement()
        else:
            alt_sequence = read.read_sequence

        cigar = read.make_cigar()

        # Parse the CIGAR string in one linear pass instead of two regex scans
        # (re.split + re.findall) — equivalent output, no regex overhead.
        cig_letters = []
        cig_numbers = []
        cig_pos = 0
        cig_len = len(cigar)
        while cig_pos < cig_len:
            num_start = cig_pos
            while cig_pos < cig_len and cigar[cig_pos].isdigit():
                cig_pos += 1
            cig_numbers.append(int(cigar[num_start:cig_pos]))
            cig_letters.append(cigar[cig_pos])
            cig_pos += 1
        cig_ops = len(cig_letters)
        next_ref_id = contig_id

        if not mate_position:
            next_pos = 0
            template_length = 0
        else:
            next_pos = mate_position

        # Pack all CIGAR ops in a single struct.pack call. Typical CIGAR has 1-3 ops,
        # so a numpy round-trip would lose to struct.pack on setup overhead.
        encoded_cig = pack(
            f'<{cig_ops}I',
            *((n << 4) | CIGAR_PACKED[l] for n, l in zip(cig_numbers, cig_letters))
        )

        # Vectorized sequence encoding via numpy lookup table. Replaces a per-byte
        # Python loop with 2x dict lookups + struct.pack per iteration, which was the
        # dominant cost in write_bam_record (~28 M Seq.__getitem__ calls on a 185 k-
        # read run). One `str()` flattens the Biopython Seq to bytes once; from there
        # the rest is vectorized.
        seq_len = read_length
        seq_str = str(alt_sequence)
        if seq_len & 1:
            seq_str = seq_str + '='
        alt_bytes = np.frombuffer(seq_str.encode('ascii'), dtype=np.uint8)
        codes = SEQ_PACKED_LOOKUP[alt_bytes]
        encoded_seq = ((codes[0::2] << 4) | codes[1::2]).tobytes()

        # Quality encoding: previously a per-base list comp `chr(n) for n in
        # quality_array` plus .encode('utf-8'). The quality_array is small ints
        # in [0, 255]; converting in one numpy cast → bytes call is much faster.
        # np.asarray coerces list inputs (used in tests) to ndarray without copy
        # when already an ndarray of the right dtype.
        encoded_qual = np.asarray(read.quality_array, dtype=np.uint8).tobytes()

        name_bytes = read.name.encode('utf-8') + b'\0'

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

        block_size = 32 + len(name_bytes) + len(encoded_cig) + len(encoded_seq) + len(encoded_qual)

        # One struct.pack for the 9 fixed-size header fields instead of nine separate
        # pack() calls plus bytes concatenation.
        header = pack(
            '<iiiIIiiii',
            block_size,
            contig_id,
            read.position,
            (read_bin << 16) + (read.mapping_quality << 8) + len(name_bytes),
            (flag << 16) + cig_ops,
            seq_len,
            next_ref_id,
            next_pos,
            template_length,
        )

        bam_handle.write(header + name_bytes + encoded_cig + encoded_seq + encoded_qual)