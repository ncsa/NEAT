import sys

from time import sleep
from tqdm import trange

from source.error_handling import print_and_log


class Tracker:
    """
    Container to keep track of the reads.
    """

    def __init__(self, reference_idx, current_progress):
        """
        This will keep track of how many reads we've sampled and overall progress
        """
        self.reference_idx = reference_idx
        self.current_progress = current_progress
        # Quickly total up the lengths of the contigs in the reference
        self.total_bp_spanned = sum([len(value) for key, value in self.reference_idx.items()])
        print(f'total bp spanned: {self.total_bp_spanned}')

        # self.finalize_progress_info()

    def update_current_progress(self):
        self.current_progress += 1

    @staticmethod
    def finalize_progress_info():
        print_and_log('Simulating reads ... FINISHED!', 'info')
        print_and_log('Simulation complete, writing files.', 'info')
