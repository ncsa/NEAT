import bisect
from source.error_handling import log_mssg, premature_exit
from source.Options import Options
from source.Models import Models
from Bio.Seq import Seq
from Bio import SeqRecord


def generate_coverage_model(coverage_dat: list, sequence: Seq = None, options: dict = None,
                            models: Models = None, degenerate: bool = False):

    # Incorporate force coverage HERE
    if degenerate:
        degenerate = True
        coverage_avg = coverage_dat[0]
    else:
        window_size, gc_scalars, target_cov_scalars = coverage_dat
        gc_cov_vals = []
        target_cov_vals = []
        avg_out = []
        coverage_distribution = []

        if len(sequence) > options.read_len:
            max_coord = len(sequence) - options.read_len
        else:
            max_coord = len(sequence)

        gc_cov_vals = compute_gc_bias(window_size, sequence, gc_scalars)

        target_cov_vals.append(target_cov_scalars[0])




def compute_gc_bias(input_sequence: Seq, size_of_window: int, gc_weight_list: list):
    """
    Takes a Biopython Sequence, a window size and a list of gc weights for the window
    and calculates the GC-bias for each base in the window.

    :param input_sequence: A biopython Seq object, which contains read data.
    :param size_of_window: An int representing the size of the window to count
    :param gc_weight_list: A list of weights to calculate any GC bias for each base
    :return:
    """
    j = 0
    return_list = []
    while j + size_of_window < len(input_sequence):
        gc_count = input_sequence[j:j + size_of_window].count('G') + \
                   input_sequence[j:j + size_of_window].count('C')
        return_list.extend([gc_weight_list[gc_count]] * size_of_window)
        j += size_of_window

    # count the remainder
    gc_count = input_sequence[-size_of_window:].count('G') + input_sequence[-size_of_window:].count('C')
    return_list.extend([gc_weight_list[gc_count]] * (len(input_sequence) - len(return_list)))

    return return_list


def compute_coverage(window: (int, int), reference: SeqRecord, chromosome: str, targeted_regions: list,
                     overlap: float, min_window_size: int,
                     models: Models, options: Options):
    """
    Compute the coverage for a sequence in a given reference, divided by window sizes.

    :param window: The window to analyze, a tuple of start and end coordinates. Chromosome is assumed at this point.
    :param reference: A SeqRecord object containing the read data for a single chromosome.
    :param chromosome: Which chromosome we are analyzing
    :param targeted_regions: A dictionary of regions in the chromosome contained in the bed file. Unless otherwise
           specified in parameters, NEAT will split coverage 98/2 for on-target/off-target regions
    :param overlap: This is the overlap size for coverage
    :param min_window_size: The minimum size we need a sequence to be to try to analyze it
    :param models: Models representing GC Bias and other coverage parameters
    :param options: NEAT gen_reads options object.
    :return:
    """
    gc_scale_value = models.gc_model
    window_size = len(gc_scale_value) - 1
    start = window[0]
    end = window[1]
    # We'll be looking for at least 50% coverage within a window
    window_average_coverage = 0.5

    coverage_data = [window_size, gc_scale_value, []]
    target_hits = 0
    # TODO work in force_coverage here
    if not options.target_bed:
        coverage_data[2] = [1.0] * (end - start)
    else:
        if not targeted_regions:
            coverage_data[2] = [options.off_target_scalar] * (end - start)
        else:
            targeted_regions.insert(0, -1)
            for j in range(start, end):
                if not (bisect.bisect(targeted_regions, j) % 2):
                    coverage_data[2].append(options.on_target_scalar)
                    target_hits += 1
                else:
                    coverage_data[2].append(options.off_target_scalar)

    """
    First check:
    Discard regions that are too small. Overlap is either the read length or, if paired ended, the fragment mean
    Added 10 to pad this out a little. Might find a parameter such as standard deviation to make it more meaningful
    
    Second Check:
    If the average coverage scalar is less than our window threshold, we'll skip this one
    This prevents areas wih for example 1 read in a window with 100 coverage and the rest 
    with 0 or something equally bad
    
    Third Check:
    If it is possible to select a fragment size that is bigger than the window, then skip the window to avoid
    erroring out.
    
    If any of these are true, return a degenerate coverage_model with avg 0.0
    """
    if (options.off_target_discard and target_hits <= overlap + 10) or \
            ((sum(coverage_data[2])/len(coverage_data[2])) < window_average_coverage) or \
            ((end - start) < min_window_size):
        log_mssg(f"Skipping segment: {chromosome}: {window}, as the region is not large enough", 'debug')
        return generate_coverage_model([0.0], degenerate=True), False

    return generate_coverage_model(coverage_data, reference.seq, options, models), True

