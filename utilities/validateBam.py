#!/usr/bin/env python
#
#
#   validateBam.py
#   Checks the BAM file for valid alignment data
#
#   Takes an input BAM file   
#
#   Usage: python validateBam.py -i /path/to/BAM_file
#
#
# Python 3 ready

import sys
import os
import _io
import gzip
import argparse
from struct import unpack
from contextlib import redirect_stdout


def get_bytes(fmt: str, amt: int, f: _io.BufferedReader):
    """
    Reads the associated number of bytes for different elements
    """
    if fmt == '<i' or fmt == '<I':
        my_size = 4
    elif fmt == '<c' or fmt == '<b' or fmt == '<B':
        my_size = 1
    else:
        print('\nError, unknown format:', fmt, '\n')
        exit(1)
    if amt == 1:
        f_read = f.read(my_size)
        if not f_read:
            return None
        return unpack(fmt, f_read)[0]

    else:
        f_read = f.read(my_size * amt)
        if not f_read:
            return None
        return unpack(fmt, f_read)


def examine_alignemnt(f: _io.BufferedReader):
    """
    Examines each block from the input file for various elements of the sequences 

    :param f: BufferedReader reading the input file  

    :return: None    
    """
    print('\nEXAMINING ALIGNMENT DATA...\n')
    aln_N = 0
    while True:
        aln_N += 1
        block_size = get_bytes('<i', 1, f)
        if block_size == None:
            break
        with open('valBAM_out.txt', 'a') as op:
            with redirect_stdout(op):
                print('[' + str(aln_N) + ']:', 'block_size:', block_size)
                print('-- refID:', get_bytes('<i', 1, f))
                print('-- pos:  ', get_bytes('<i', 1, f))
                bmqnl = get_bytes('<I', 1, f)
                binv = (bmqnl >> 16) & 65535
                mapq = (bmqnl >> 8) & 255
                lrn = bmqnl & 255
                print('-- bmqnl:', bmqnl, '(bin=' + str(binv) + ', mapq=' + str(mapq) + ', l_readname+1=' + str(lrn) + ')')
                flgnc = get_bytes('<I', 1, f)
                flag = (flgnc >> 16) & 65535
                ncig = flgnc & 65535
                print('-- flgnc:', flgnc, '(flag=' + str(flag) + ', ncig=' + str(ncig) + ')')
                print('-- l_seq:', get_bytes('<i', 1, f))
                print('-- nxtID:', get_bytes('<i', 1, f))
                print('-- nxtPo:', get_bytes('<i', 1, f))
                print('-- tlen: ', get_bytes('<i', 1, f))
                print('-- rname:', str([f.read(lrn)])[1:-1])

        f.read(block_size - 32 - lrn)


def func_parser() -> argparse.Namespace:
    """
    Defines what arguments the program requires, and argparse will figure out how to parse those out of sys.argv

    :return: an instance of the argparse class that can be used to access command line arguments
    """

    parser = argparse.ArgumentParser(description='validateBam.py')
    parser.add_argument('-i', type=str, required=True, metavar='<str>', help="input_file.bam")
    args = parser.parse_args()

    return args


def main():
    """
    Validates a given BAM file for the alignment

    :return: Generates an output file - valBAM_out.txt containing references read from BAM file
    """

    args = func_parser()

    BAM_EOF = ['1f', '8b', '08', '04', '00', '00', '00', '00', '00', 'ff', '06', '00', '42', '43', '02', '00', '1b', '00',
        '03', '00', '00', '00', '00', '00', '00', '00', '00', '00']

    # check eof
    IN_BAM = args.i
    f = _io.open(IN_BAM, 'rb')
    

    f.seek(os.path.getsize(IN_BAM) - 28)
    EOF = [format(n, '02x') for n in f.read()]
    print('EOF_MARKER:  ', ' '.join(EOF))
    if EOF != BAM_EOF:
        print('\nWARNING: BAM EOF DOES NOT MATCH EXPECTED STRING.\n')
    f.close()

    # check other stuff
    f = gzip.open(IN_BAM, 'rb')


    with open('valBAM_out.txt', 'w') as opf:
        with redirect_stdout(opf):
            print('MAGIC STRING:', f.read(4))
            l_text = get_bytes('<i', 1, f)
            print('l_text:      ', l_text)
            print('text:      \n', f.read(l_text))
            n_ref = get_bytes('<i', 1, f)
            print('n_ref:       ', n_ref)

    for i in range(n_ref):
        l_name = get_bytes('<i', 1, f)
        with open('valBAM_out.txt', 'a') as op:
            with redirect_stdout(op):
                print('ref' + str(i) + ' - l_name:', l_name)
                print('ref' + str(i) + ' - name:  ', f.read(l_name))
                print('ref' + str(i) + ' - l_ref: ', get_bytes('<i', 1, f))

    examine_alignemnt(f)
    print('Verified!')

    f.close()


if __name__ == '__main__':
    main()
