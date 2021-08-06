#!/usr/bin/env python
#
#
#   compute_fraglen.py   
#   Compute Fragment Length Model for gen_reads.py
#                  
#   Takes as input path to a BAM or SAM file and outputs a pickle file with fragment length statistics for input
#   into gen_reads.py
#
#   Usage: python compute_fraglen.py -i /path/to/file.bam(.sam) -o output_prefix
#
#
# Python 3 ready

import pickle
import argparse
import platform
import sys
import os.path

ops = platform.system()
if ops != 'Windows':
    import pysam
else:
    print("Pysam is a required module for this utility and currently does not have a Windows version. Please re-run"
          "with a Linux-based OS (e.g., MacOS or Ubuntu).")


def median(datalist: list) -> float:
    """
    Finds the median of a list of data. For this function, the data are expected to be a list of
    numbers, either float or int.
    :param datalist: the list of data to find the median of. This should be a set of numbers.
    :return: The median of the set
    >>> median([2])
    2
    >>> median([2183, 2292, 4064, 4795, 7471, 12766, 14603, 15182, 16803, 18704, 21504, 21677, 23347, 23586, 24612, 24878, 25310, 25993, 26448, 28018, 28352, 28373, 28786, 30037, 31659, 31786, 33487, 33531, 34442, 39138, 39718, 39815, 41518, 41934, 43301])
    25993
    >>> median([1,2,4,6,8,12,14,15,17,21])
    10.0
    """
    # using integer division here gives the index of the midpoint, due to zero-based indexing.
    midpoint = len(datalist)//2

    # Once we've found the midpoint, we calculate the median, which is just the middle value if there are an
    # odd number of values, or the average of the two middle values if there are an even number
    return_median = datalist[midpoint]
    if len(datalist) % 2 == 0:
        return_median = (datalist[midpoint] + datalist[midpoint-1])/2

    return return_median


def median_absolute_deviation(datalist: list) -> list:
    """
    Calculates the absolute value of the median deviation from the median for each element of of a datalist. 
    Then returns the median of these values.
    :param datalist: A list of data to find the MAD of
    :return: index of median of the deviations
    >>> median_absolute_deviation([2183, 2292, 4064, 4795, 7471, 12766, 14603, 15182, 16803, 18704, 21504, 21677, 23347, 23586, 24612, 24878, 25310, 25993, 26448, 28018, 28352, 28373, 28786, 30037, 31659, 31786, 33487, 33531, 34442, 39138, 39718, 39815, 41518, 41934, 43301])
    7494
    >>> median_absolute_deviation([1,2,4,6,8,12,14,15,17,21])
    5.5
    >>> median_absolute_deviation([0,2])
    1.0
    """
    my_median = median(datalist)
    deviations = []
    for item in datalist:
        # We take the absolute difference between the value and the median
        X_value = abs(item - my_median)
        # This creates a dataset that is the absolute deviations about the median
        deviations.append(X_value)
    # The median of the absolute deviations is the median absolute deviation
    return median(sorted(deviations))


def count_frags(file: str) -> list:
    """
    Takes a sam file input and creates a list of the number of reads that are paired,
    first in the pair, confidently mapped and whose pair is mapped to the same reference
    :param file: A sam input file
    :return: A list of the tlens from the bam/sam file
    """

    filter_mapqual = 10  # only consider reads that are mapped with at least this mapping quality
    count_list = []

    try:
        file_to_parse = pysam.AlignmentFile(file, 'rb')
    except ValueError:
        print("File type not recognized by pysam.")
        sys.exit(1)

    for item in file_to_parse.fetch():

        my_template_length = abs(item.template_length)

        # if read is paired, and is first in pair, and is confidently mapped...
        if item.flag & 1 and item.flag & 64 and item.mapping_quality > filter_mapqual:
            # and mate is mapped to same reference
            if item.next_reference_id == '=' or item.next_reference_id == item.reference_id:
                count_list.append(my_template_length)
    count_list = sorted(count_list)
    file_to_parse.close()
    return count_list


def compute_probs(datalist: list) -> (list, list):
    """
    Computes the probabilities for fragments with at least 100 pairs supporting it and that are at least 10 median
    deviations from the median.
    :param datalist: A list of fragments with counts
    :return: A list of values that meet the criteria and a list of their associated probabilities
    """
    FILTER_MINREADS = 100  # only consider fragment lengths that have at least this many read pairs supporting it
    FILTER_MEDDEV_M = 10  # only consider fragment lengths this many median deviations above the median
    values = []
    probabilities = []
    med = median(datalist)
    mad = median_absolute_deviation(datalist)

    for item in list(set(datalist)):
        if 0 < item <= med + FILTER_MEDDEV_M * mad:
            data_count = datalist.count(item)
            if data_count >= FILTER_MINREADS:
                values.append(item)
                probabilities.append(data_count)
    count_sum = float(sum(probabilities))
    probabilities = [n / count_sum for n in probabilities]
    return values, probabilities


def func_parser() -> argparse.Namespace:
    """
    Defines what arguments the program requires, and argparse will figure out how to parse those out of sys.argv

    :return: an instance of the argparse class that can be used to access command line arguments
    """

    parser = argparse.ArgumentParser(description="compute_fraglen.py",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter,)
    parser.add_argument('-i', type=str, metavar="input", required=True, default=None,
                        help="Input file name - file.bam(.sam)")
    parser.add_argument('-o', type=str, metavar="output", required=True, default=None, help="Prefix for output")

    args = parser.parse_args()

    return args


def main():
    """
    Main function takes 2 arguments:
        input - a path to a sam or bam file input. Note that sam files can be formed by applying samtools to a bam file
        in the following way: samtools view nameof.bam > nameof.sam
        
        output - the string prefix of the output. The actual output will be the prefix plus ".p" at the end
        for pickle file. The list of values and list of probabilities are dumped as a list of lists
        into a pickle file on completion of the analysis

    :return: None
    """

    args = func_parser() 

    checkinp = args.i

     
    if checkinp[-4:] == '.bam':
        if not os.path.isfile(checkinp+'.bai'):
            command = 'samtools index ' + checkinp
            os.system(command)
    
    input_file = checkinp

    output_prefix = args.o
    output = output_prefix + '.p'

    all_tlens = count_frags(input_file)
    print('\nSaving model...')
    out_vals, out_probs = compute_probs(all_tlens)
    # Print statements for debugging:
    # print(f'Out vals: {out_vals}')
    # print(f'out probs: {out_probs}')
    pickle.dump([out_vals, out_probs], open(output, 'wb'))
    print('\nModel successfully saved.')


if __name__ == "__main__":
    main()
