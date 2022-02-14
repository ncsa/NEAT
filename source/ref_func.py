import sys
import time
import os
import gzip
import pathlib
import random
import logging
from Bio.Seq import Seq
from Bio.Seq import MutableSeq

from source.error_handling import premature_exit
from source.constants_and_models import ALLOWED_NUCL, OK_CHR_ORD


def index_ref(reference_path: str) -> list:
    """
    Index reference fasta
    :param reference_path: string path to the reference
    :return: reference index in list from
    """
    tt = time.time()

    absolute_reference_location = pathlib.Path(reference_path).resolve()

    # sanity check
    if not absolute_reference_location.is_file():
        print("\nProblem reading the reference fasta file.\n")
        logging.error("Problem reading the reference fasta file.")
        premature_exit(1)

    index_filename = None

    # check if the reference file already exists
    if absolute_reference_location.with_suffix('.fai').is_file():
        print('found index ' + str(absolute_reference_location.with_suffix('.fai')))
        index_filename = absolute_reference_location.with_suffix('.fai')
    elif absolute_reference_location.with_suffix(absolute_reference_location.suffix + '.fai').is_file():
        print('found index ' +
              str(absolute_reference_location.with_suffix(absolute_reference_location.suffix + '.fai')))
        index_filename = absolute_reference_location.with_suffix(absolute_reference_location.suffix + '.fai')
    else:
        pass

    ref_indices = []
    if index_filename is not None:
        fai = open(index_filename, 'r')
        for line in fai:
            splt = line.strip().split('\t')
            # Defined as the number of bases in the contig
            seq_len = int(splt[1])
            # Defined as the byte index where the contig sequence begins
            offset = int(splt[2])
            # Defined as bases per line in the Fasta file
            line_ln = int(splt[3])
            n_lines = seq_len // line_ln
            if seq_len % line_ln != 0:
                n_lines += 1
            # Item 3 in this gives you the byte position of the next contig, I believe
            ref_indices.append((splt[0], offset, offset + seq_len + n_lines, seq_len))
        fai.close()
        return ref_indices

    print('Index not found, creating one... ')
    if absolute_reference_location.suffix == ".gz":
        ref_file = gzip.open(absolute_reference_location, 'rt')
    else:
        ref_file = open(absolute_reference_location, 'r')
    prev_r = None
    prev_p = None
    seq_len = 0

    while True:
        data = ref_file.readline()
        if not data:
            ref_indices.append((prev_r, prev_p, ref_file.tell() - len(data), seq_len))
            break
        elif data[0] == '>':
            if prev_p is not None:
                ref_indices.append((prev_r, prev_p, ref_file.tell() - len(data), seq_len))
            seq_len = 0
            prev_p = ref_file.tell()
            prev_r = data[1:].strip().split()

        else:
            seq_len += len(data) - 1
    ref_file.close()

    print('{0:.3f} (sec)'.format(time.time() - tt))
    return ref_indices


def read_ref(ref_path, ref_inds_i, n_handling, quiet=False):
    tt = time.time()
    if not quiet:
        print('reading ' + ref_inds_i[0] + '... ')

    absolute_reference_path = pathlib.Path(ref_path)
    if absolute_reference_path.suffix == '.gz':
        ref_file = gzip.open(absolute_reference_path, 'rt')
    else:
        ref_file = open(absolute_reference_path, 'r')

    # TODO convert to SeqIO containers
    # for seq_record in SeqIO.parse(ref_file, "fasta"):
    #     pass

    ref_file.seek(ref_inds_i[1])
    my_dat = ''.join(ref_file.read(ref_inds_i[2] - ref_inds_i[1]).split('\n'))
    my_dat = Seq(my_dat.upper())
    # Mutable seqs have a number of disadvantages. I'm going to try making them immutable and see if that helps
    # my_dat = MutableSeq(my_dat)

    # find N regions
    # data explanation: my_dat[n_atlas[0][0]:n_atlas[0][1]] = solid block of Ns
    prev_ni = 0
    n_count = 0
    n_atlas = []
    for i in range(len(my_dat)):
        if my_dat[i] == 'N' or my_dat[i] not in ALLOWED_NUCL:
            if n_count == 0:
                prev_ni = i
            n_count += 1
            if i == len(my_dat) - 1:
                n_atlas.append((prev_ni, prev_ni + n_count))
        else:
            if n_count > 0:
                n_atlas.append((prev_ni, prev_ni + n_count))
            n_count = 0

    # handle N base-calls as desired
    # TODO this seems to randomly replace an N with a base. Is this necessary? How to do this in an immutable seq?
    n_info = {'all': [], 'big': [], 'non_N': []}
    if n_handling[0] == 'random':
        for region in n_atlas:
            n_info['all'].extend(region)
            if region[1] - region[0] <= n_handling[1]:
                for i in range(region[0], region[1]):
                    temp = MutableSeq(my_dat)
                    temp[i] = random.choice(ALLOWED_NUCL)
                    my_dat = Seq(temp)
            else:
                n_info['big'].extend(region)
    elif n_handling[0] == 'allChr' and n_handling[2] in OK_CHR_ORD:
        for region in n_atlas:
            n_info['all'].extend(region)
            if region[1] - region[0] <= n_handling[1]:
                for i in range(region[0], region[1]):
                    temp = MutableSeq(my_dat)
                    temp[i] = n_handling[2]
                    my_dat = Seq(temp)
            else:
                n_info['big'].extend(region)
    elif n_handling[0] == 'ignore':
        for region in n_atlas:
            n_info['all'].extend(region)
            n_info['big'].extend(region)
    else:
        print('\nERROR: UNKNOWN N_HANDLING MODE\n')
        logging.error("UNKNOWN N_HANDLING MODE")
        premature_exit(1)

    habitable_regions = []
    if not n_info['big']:
        n_info['non_N'] = [(0, len(my_dat))]
    else:
        for i in range(0, len(n_info['big']), 2):
            if i == 0:
                habitable_regions.append((0, n_info['big'][0]))
            else:
                habitable_regions.append((n_info['big'][i - 1], n_info['big'][i]))
        habitable_regions.append((n_info['big'][-1], len(my_dat)))
    for n in habitable_regions:
        if n[0] != n[1]:
            n_info['non_N'].append(n)

    ref_file.close()

    if not quiet:
        print('{0:.3f} (sec)'.format(time.time() - tt))

    return my_dat, n_info


def find_n_regions(input_sequence: Seq, n_handling: tuple, n_length_allowed: int = 25) -> dict:
    """
    Finds N regions in the sequence
    :param input_sequence: Biopython Seq object containing the sequence to scan.
    :param n_handling: tuple describing the n handling parameters
    :return:

    >>> myseq = Seq("NNNNNAAACCCTTTNAAAACCCCNNNNN")
    >>> n_handling = ('')
    """
    # data explanation: my_dat[n_atlas[0][0]:n_atlas[0][1]] = solid block of Ns
    prev_ni = 0
    n_count = 0
    n_atlas = []
    my_dat = input_sequence.replace(x, "N")
    for i in range(len(my_dat)):
        if my_dat[i] == 'N' or my_dat[i] not in ALLOWED_NUCL:
            if n_count == 0:
                prev_ni = i
            n_count += 1
            if i == len(my_dat) - 1:
                n_atlas.append((prev_ni, prev_ni + n_count))
        else:
            if n_count > 0:
                n_atlas.append((prev_ni, prev_ni + n_count))
            n_count = 0

    # handle N base-calls as desired
    # random is for paired end, ignore for single end allChr is never used.
    n_info = {'all': [], 'big': [], 'non_N': []}
    if n_handling[0] == 'random':
        for region in n_atlas:
            n_info['all'].extend(region)
            if region[1] - region[0] <= n_handling[1]:
                for i in range(region[0], region[1]):
                    my_dat[i] = random.choice(ALLOWED_NUCL)
            else:
                n_info['big'].extend(region)

    # this block of code is never used because 'allChr' is not used.
    # The current code options are n_handling = ('random', fragment_size) and
    # n_handling = ('ignore', read_len). Looking at this block it would replace
    # all N's with whatever base is in the third position of n_handling, so it  would
    # have to look something like n_handling = ('allChr', fragment_size, 'A').
    # Bottom line, this part may be deletable
    elif n_handling[0] == 'allChr' and n_handling[2] in ALLOWED_NUCL:
        for region in n_atlas:
            n_info['all'].extend(region)
            if region[1] - region[0] <= n_handling[1]:
                for i in range(region[0], region[1]):
                    my_dat[i] = n_handling[2]
            else:
                n_info['big'].extend(region)

    elif n_handling[0] == 'ignore':
        for region in n_atlas:
            n_info['all'].extend(region)
            n_info['big'].extend(region)
    else:
        print('\nERROR: UNKNOWN N_HANDLING MODE')
        logging.error("UNKNOWN N_HANDLING MODE")
        premature_exit(1)

    habitable_regions = []
    if not n_info['big']:
        n_info['non_N'] = [(0, len(my_dat))]
    else:
        for i in range(0, len(n_info['big']), 2):
            if i == 0:
                habitable_regions.append((0, n_info['big'][0]))
            else:
                habitable_regions.append((n_info['big'][i - 1], n_info['big'][i]))
        habitable_regions.append((n_info['big'][-1], len(my_dat)))
    for n in habitable_regions:
        if n[0] != n[1]:
            n_info['non_N'].append(n)

    return n_info
