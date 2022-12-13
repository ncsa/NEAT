import logging

from Bio.Seq import Seq

from ...models import MutationModel
from ...common import ALLOWED_NUCL, TRI_IND

_LOG = logging.getLogger(__name__)

__all__ = [
    "find_file_breaks",
    "map_chromosome"
]


def find_file_breaks(threads: int, mode: str, reference_index: dict) -> dict:
    """
    Returns a dictionary with the chromosomes as keys
    For the chrom method, the value for each key will just  be "all"
    whereas for subdivison, the value for each key should be a list of indices that partition the
    sequence into roughly equal sizes.
    :param threads: number of threads for this run
    :param mode: partition mode for this run (chrom or subdivision)
    :param reference_index: a dictionary with chromosome keys and sequence values
    :return: a dictionary containing the chromosomes as keys and either "all" for values, or a list of indices
    """
    partitions = {}
    if mode.lower() == "chrom" or threads == 1:
        for contig in reference_index.keys():
            partitions[contig] = [(0, len(reference_index[contig]))]
    elif mode.lower() == "subdivision":
        # Add items one at a time to partition list until the total length is greater than delta.
        for contig in reference_index:
            contig_length = len(reference_index[contig])
            delta = contig_length // threads

            if contig not in partitions:
                partitions[contig] = []

            breakpoints = list(range(0, contig_length, delta))

            # Since we know the first index will be zero, we can skip the first item in the breakpoints list
            # And since we want the last partition to grab the rest, we'll stop short and add it manually
            for index in breakpoints:
                if index + delta <= contig_length:
                    partitions[contig].append((index, index + delta))
                else:
                    # Have to extend the last one so we don't get a tiny read we can't process
                    partitions[contig][-1] = (partitions[contig][-1][0], contig_length)

        _LOG.debug(f'breaks = {partitions}')
    else:
        raise ValueError("Invalid partition mode. Must be either chrom or subdivision.")

    return partitions


def map_chromosome(sequence: Seq, mut_model: MutationModel):
    """
    This will create a trinucleotide probability map for the chromosome of interest. Creating this on the fly
    was too slow, so we'll try frontloading the work

    :param sequence: the sequence to model
    :param mut_model: the models for this simulation
    :return: A map of the trinuc probs
    """
    # Set up the model dictionary
    trinuc_models = [0.0] * len(sequence)
    # Start at +1 so we get a trinucleotide to start, and end one shy for the same reason
    for i in range(1, len(sequence) - 1):
        trinuc = sequence[i-1: i + 2].seq
        # Let's double check to make sure we didn't pick up a stray N
        if any([j for j in trinuc if j not in ALLOWED_NUCL]):
            # leave the probability at 0
            continue
        trinuc_models[i] = mut_model.trinuc_trans_bias[TRI_IND[trinuc]]

    return trinuc_models

