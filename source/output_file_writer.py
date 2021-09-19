from struct import pack
import Bio.bgzf as bgzf
import pathlib
import re
import sys
import logging

from source.neat_cigar import CigarString

# TODO make this a configurable option
BAM_COMPRESSION_LEVEL = 6


def reverse_complement(dna_string) -> str:
    """
    Return the reverse complement of a string from a DNA strand. Found this method that is slightly faster than
    biopython. Thanks to this stack exchange post:
    https://bioinformatics.stackexchange.com/questions/3583/what-is-the-fastest-way-to-get-the-reverse-complement-of-a-dna-sequence-in-pytho
    :param dna_string: string of DNA, either in string or Seq format
    :return: the reverse complement of the above string in either string or MutableSeq format
    """
    if type(dna_string) != str:
        dna_string.reverse_complement()
        return dna_string
    else:
        tab = str.maketrans("ACTGN", "TGACN")

        return dna_string.translate(tab)[::-1]


# SAMtools reg2bin function
def reg2bin(beg: int, end: int):
    """
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
    string_list = list(set(string_list))
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


CIGAR_PACKED = {'M': 0, 'I': 1, 'D': 2, 'N': 3, 'S': 4, 'H': 5, 'P': 6, '=': 7, 'X': 8}
SEQ_PACKED = {'=': 0, 'A': 1, 'C': 2, 'M': 3, 'G': 4, 'R': 5, 'S': 6, 'V': 7,
              'T': 8, 'W': 9, 'Y': 10, 'H': 11, 'K': 12, 'D': 13, 'B': 14, 'N': 15}

# TODO figure out an optimum batch size
BUFFER_BATCH_SIZE = 8000  # write out to file after this many reads


# TODO find a better way to write output files
class OutputFileWriter:
    def __init__(self, out_prefix, paired=False, bam_header=None, vcf_header=None,
                 write_fastq=True, write_fasta=False, write_bam=False, write_vcf=False):

        # TODO Eliminate paired end as an option for fastas. Plan is to create a write fasta method.
        self.write_fastq = write_fastq
        self.write_fasta = write_fasta
        self.write_bam = write_bam
        self.write_vcf = write_vcf
        output1 = None
        output2 = None
        if self.write_fasta:
            output1 = out_prefix.parent / (out_prefix.name + '.fasta.gz')
        elif paired and self.write_fastq:
            output1 = out_prefix.parent / (out_prefix.name + '_read1.fq.gz')
            output2 = out_prefix.parent / (out_prefix.name + '_read2.fq.gz')
        elif self.write_fastq:
            output1 = out_prefix.parent / (out_prefix.name + '.fq.gz')
        if self.write_bam:
            bam = out_prefix.parent / (out_prefix.name + '_golden.bam')
        if self.write_vcf:
            vcf = out_prefix.parent / (out_prefix.name + '_golden.vcf.gz')

        # TODO Make a fasta-specific method
        self.output1_file = None
        self.output2_file = None
        if output1:
            self.output1_file = bgzf.open(output1, 'w')
        if output2:
            self.output2_file = bgzf.open(output2, 'w')

        # VCF OUTPUT
        self.vcf_file = None
        # Assuming we wanted the vcf and there's nothing wrong with the header, then we will proceed with creating the
        # vcf file. If the header is empty and we wanted it, that is a bug we need to catch right here.
        if self.write_vcf and not vcf_header:
            print("ERROR: Something wrong with VCF header.")
            logging.error("Something wrong with VCF header.")
            sys.exit(1)

        if self.write_vcf:
            self.vcf_file = bgzf.open(vcf, 'wb')

            # WRITE VCF HEADER
            self.vcf_file.write('##fileformat=VCFv4.1\n'.encode('utf-8'))
            reference = '##reference=' + vcf_header[0] + '\n'
            self.vcf_file.write(reference.encode('utf-8'))
            self.vcf_file.write('##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">\n'.encode('utf-8'))
            self.vcf_file.write(
                '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n'.encode('utf-8'))
            self.vcf_file.write(
                '##INFO=<ID=VMX,Number=1,Type=String,Description="SNP is Missense in these Read Frames">\n'.encode(
                    'utf-8'))
            self.vcf_file.write(
                '##INFO=<ID=VNX,Number=1,Type=String,Description="SNP is Nonsense in these Read Frames">\n'.encode(
                    'utf-8'))
            self.vcf_file.write(
                '##INFO=<ID=VFX,Number=1,Type=String,Description="Indel Causes Frameshift">\n'.encode('utf-8'))
            self.vcf_file.write(
                '##INFO=<ID=WP,Number=A,Type=Integer,Description="NEAT-GenReads ploidy indicator">\n'.encode(
                    'utf-8'))
            self.vcf_file.write('##ALT=<ID=DEL,Description="Deletion">\n'.encode('utf-8'))
            self.vcf_file.write('##ALT=<ID=DUP,Description="Duplication">\n'.encode('utf-8'))
            self.vcf_file.write('##ALT=<ID=INS,Description="Insertion of novel sequence">\n'.encode('utf-8'))
            self.vcf_file.write('##ALT=<ID=INV,Description="Inversion">\n'.encode('utf-8'))
            self.vcf_file.write('##ALT=<ID=CNV,Description="Copy number variable region">\n'.encode('utf-8'))
            self.vcf_file.write('##ALT=<ID=TRANS,Description="Translocation">\n'.encode('utf-8'))
            self.vcf_file.write('##ALT=<ID=INV-TRANS,Description="Inverted translocation">\n'.encode('utf-8'))
            # TODO add sample to vcf output
            self.vcf_file.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'.encode('utf-8'))

        # BAM OUTPUT
        self.bam_file = None
        # Assuming we wanted the bam and there's nothing wrong with the header, then we will proceed with creating the
        # bam file. If the header is empty and we wanted it, that is a bug we need to catch right here.
        if self.write_bam and not bam_header:
            print("Something wrong with VCF header.")
            logging.error("")
            sys.exit(1)

        if self.write_bam:
            self.bam_file = bgzf.BgzfWriter(bam, 'w', compresslevel=BAM_COMPRESSION_LEVEL)

            # WRITE BAM HEADER
            self.bam_file.write("BAM\1")
            header = '@HD\tVN:1.5\tSO:coordinate\n'
            for key in bam_header.keys():
                header += '@SQ\tSN:' + str(key) + '\tLN:' + str(len(bam_header[key])) + '\n'
            header += '@RG\tID:NEAT\tSM:NEAT\tLB:NEAT\tPL:NEAT\n'
            header_bytes = len(header)
            num_refs = len(bam_header)
            self.bam_file.write(pack('<i', header_bytes))
            self.bam_file.write(header)
            self.bam_file.write(pack('<i', num_refs))

            for key in bam_header.keys():
                l_name = len(key) + 1
                self.bam_file.write(pack('<i', l_name))
                self.bam_file.write(key + '\0')
                self.bam_file.write(pack('<i', len(bam_header[key])))

        # buffers for more efficient writing
        self.output1_buffer = []
        self.output2_buffer = []
        self.bam_buffer = []

    def write_fasta_record(self, read_name, read):
        """
        Needs to take the input given and make a fasta record. We need to figure out how to get it into standard fasta
        format for this to work
        """
        # TODO This will make a technically correct but not standard fasta file.
        #  We either need to recombine this in another step, or think of a better way to do this
        self.output1_buffer.append('>' + read_name + '/1\n' + str(read) + '\n')

    def write_fastq_record(self, read_name, read1, quality1, read2=None, quality2=None, orientation=None):
        # Since read1 and read2 are Seq objects from Biopython, they have reverse_complement methods built-in
        (read1, quality1) = (read1, quality1)
        if read2 is not None and orientation is True:
            (read2, quality2) = (read2.reverse_complement(), quality2[::-1])
        elif read2 is not None and orientation is False:
            read2_tmp = read2
            qual2_tmp = quality2
            (read2, quality2) = (read1, quality1)
            (read1, quality1) = (read2_tmp.reverse_complement(), qual2_tmp[::-1])

        if self.write_fasta:
            self.output1_buffer.append('>' + read_name + '/1\n' + str(read1) + '\n')
            if read2 is not None:
                self.output2_buffer.append('>' + read_name + '/2\n' + str(read2) + '\n')
        else:
            self.output1_buffer.append('@' + read_name + '/1\n' + str(read1) + '\n+\n' + quality1 + '\n')
            if read2 is not None:
                self.output2_buffer.append('@' + read_name + '/2\n' + str(read2) + '\n+\n' + quality2 + '\n')

    def write_vcf_record(self, chrom, pos, id_str, ref, alt, qual, filt, info):
        self.vcf_file.write(
            str(chrom) + '\t' + str(pos) + '\t' + str(id_str) + '\t' + str(ref) + '\t' + str(alt) + '\t' + str(
                qual) + '\t' + str(filt) + '\t' + str(info) + '\n')

    def write_bam_record(self, chromosome_index, read_name, pos_0, cigar, seq, qual, output_sam_flag,
                         mate_pos=None, aln_map_quality=70):

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
            # print(seq[2*i], seq[2*i+1])
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
        self.bam_buffer.append(
            (chromosome_index, pos_0, pack('<i', block_size) + pack('<i', chromosome_index) + pack('<i', pos_0) +
             pack('<I', (my_bin << 16) + (my_map_quality << 8) + len(read_name) + 1) +
             pack('<I', (output_sam_flag << 16) + cig_ops) + pack('<i', seq_len) +
             pack('<i', next_ref_id) +
             pack('<i', next_pos) + pack('<i', my_t_len) + read_name.encode('utf-8') +
             b'\0' + encoded_cig + encoded_seq + encoded_qual.encode('utf-8')))

    def flush_buffers(self, bam_max=None, last_time=False):
        if (len(self.output1_buffer) >= BUFFER_BATCH_SIZE or len(self.bam_buffer) >= BUFFER_BATCH_SIZE) or (
                len(self.output1_buffer) and last_time) or (len(self.bam_buffer) and last_time):
            # fasta
            # TODO this is a potential place to reorg the fasta file before it's written
            if self.write_fasta:
                self.output1_file.write(''.join.output1_buffer)

            # fastq
            elif self.write_fastq:
                self.output1_file.write(''.join(self.output1_buffer))
                if len(self.output2_buffer):
                    self.output2_file.write(''.join(self.output2_buffer))

            # bam
            if len(self.bam_buffer):
                bam_data = sorted(self.bam_buffer)
                if last_time:
                    self.bam_file.write(b''.join([n[2] for n in bam_data]))
                    self.bam_buffer = []
                else:
                    ind_to_stop_at = 0
                    for i in range(0, len(bam_data)):
                        # if we are from previous reference, or have coordinates lower
                        # than next window position, it's safe to write out to file
                        if bam_data[i][0] != bam_data[-1][0] or bam_data[i][1] < bam_max:
                            ind_to_stop_at = i + 1
                        else:
                            break
                    self.bam_file.write(b''.join([n[2] for n in bam_data[:ind_to_stop_at]]))
                    # Debug statement
                    # print(f'BAM WRITING: {ind_to_stop_at}/{len(bam_data)}')
                    if ind_to_stop_at >= len(bam_data):
                        self.bam_buffer = []
                    else:
                        self.bam_buffer = bam_data[ind_to_stop_at:]
            self.output1_buffer = []
            self.output2_buffer = []

    def close_files(self):
        self.flush_buffers(last_time=True)
        if self.output1_file:
            self.output1_file.close()
        if self.output2_file:
            self.output2_file.close()
        if self.vcf_file:
            self.vcf_file.close()
        if self.bam_file:
            self.bam_file.close()
