#!/usr/bin/env source

import gen_reads_runner
import argparse
import pathlib
import gzip
import shutil
import sys
import copy
from time import time


class Bacterium:
    def __init__(self, reference: str, name: str, chrom_names: list):
        """
        Class Bacterium for keeping track of all the elements of a bacterium for the purposes of this simulation
        :param reference: The str representing the location of the reference file
        :param name: The name of this particular bacterium.
        :param chrom_names: The list of chromosome names from the progenitor bacterium
        """
        self.reference = pathlib.Path(reference)
        self.name = name
        self.chroms = chrom_names
        # Temporarily set the reference as the bacterium's file, until it is analyzed
        self.file = pathlib.Path(reference)
        self.analyze()

    def __repr__(self):
        return str(self.name)

    def __str__(self):
        return str(self.name)

    def get_file(self):
        return self.file

    def get_chroms(self):
        return self.chroms

    def analyze(self):
        """
        This function is supposed to just run genreads for the bacterium, but doing so requires some file
        manipulation to unzip the file and fix genreads horribly formatted fasta file.
        :return: None
        """
        args = ['-r', str(self.reference), '-R', '101', '-o', self.name, '--fa', '-c', '1']
        gen_reads_runner.main(args)
        self.file = pathlib.Path().absolute() / (self.name + ".fasta.gz")

        # The following workaround is due to the fact that genReads cannot handle gzipped
        # fasta files, so we have to unzip it for it to actually work.
        unzipped_path = pathlib.Path().absolute() / (self.name + ".fasta")
        unzip_file(self.file, unzipped_path)
        pathlib.Path.unlink(pathlib.Path().absolute() / (self.name + ".fasta.gz"))  # deletes unused zip file
        self.file = unzipped_path
        # end workaround

        # Now we further have to fix the fasta file, which outputs in a form that doesn't make much sense,
        # so that it can be properly analyzed in the next generation by genreads.
        temp_name_list = copy.deepcopy(self.chroms)
        temp_file = self.file.parents[0] / 'neat_temporary_fasta_file.fa'
        temp_file.touch()
        chromosome_name = ""
        sequence = ""
        with self.file.open() as f:
            for line in f:
                if line.startswith(">"):
                    for name in temp_name_list:
                        if name in line:
                            if chromosome_name != ">" + name + "\n":
                                if sequence:
                                    temp_file.open('a').write(chromosome_name + sequence)
                                    sequence = ""
                                chromosome_name = ">" + name + "\n"
                                temp_name_list.remove(name)
                            else:
                                continue
                    if not chromosome_name:
                        print("Something went wrong with the generated fasta file.\n")
                        sys.exit(1)
                else:
                    sequence = sequence + line
        temp_file.open('a').write(chromosome_name + sequence)
        shutil.copy(temp_file, self.file)
        pathlib.Path.unlink(temp_file)

    def sample(self, coverage_value: int, fragment_size: int, fragment_std: int):
        """
        This function simple runs genreads on the file associated with this bacterium
        :param coverage_value: What depth of coverage to sample the reads at.
        :param fragment_size: The mean insert size for the resultant fastqs.
        :param fragment_std: The standard deviation of the insert size
        :return: None
        """
        args = ['-r', str(self.file), '-M', '0', '-R', '101', '-o', self.name,
                '-c', str(coverage_value), '--pe', str(fragment_size), str(fragment_std), '--vcf', '--bam']

        gen_reads_runner.main(args)

    def remove(self):
        """
        This function simple deletes the file associated with this bacterium, or raises an error if there is a problem
        :return: None
        """
        try:
            pathlib.Path.unlink(self.file)
        except FileExistsError:
            print('\nThere was a problem deleting a file\n')
            raise FileExistsError()


def unzip_file(zipped_file: pathlib, unzipped_file: pathlib):
    """
    This unzips a gzipped file, then saves the unzipped file as a new file.
    :param zipped_file: pathlib object that points to the zipped file
    :param unzipped_file: pathlib object that points to the unzipped file
    :return: None
    """
    with gzip.open(zipped_file, 'rb') as f_in:
        with open(unzipped_file, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)


def cull(population: list, percentage: float = 0.5) -> list:
    """
    The purpose of this function will be to cull the bacteria created in the model
    :param percentage: percentage of the population to eliminate
    :param population: the list of members to cull
    :return: The list of remaining members
    """
    cull_amount = round(len(population) * percentage)
    print("Culling {} members from population".format(cull_amount))
    for _ in range(cull_amount):
        selection = rng.choice(population)
        population.remove(selection)
        selection.remove()
    return population


def initialize_population(reference: str, pop_size: int, chrom_names: list) -> list:
    """
    The purpose of this function is to evolve the initial population of bacteria. All bacteria are stored as
    Bacterium objects.
    :param chrom_names: A list of contigs from the original fasta
    :param reference: string path to the reference fasta file
    :param pop_size: size of the population to initialize.
    :return population: returns a list of bacterium objects.
    """
    names = []
    for j in range(pop_size):
        names.append("bacterium_0_{}".format(j+1))
    population = []
    for i in range(pop_size):
        new_member = Bacterium(reference, names[i], chrom_names)
        population.append(new_member)
    return population


def evolve_population(population: list, generation: int) -> list:
    """
    This evolves an existing population by doubling them (binary fission), then introducing random mutation to
    each member of the population.
    :param generation: Helps determine the starting point of the numbering system so the bacteria have unique names
    :param population: A list of fasta files representing the bacteria.
    :return: None
    """
    children_population = population + population
    names = []
    new_population = []
    for j in range(len(children_population)):
        names.append("bacterium_{}_{}".format(generation, j+1))
    for i in range(len(children_population)):
        child = Bacterium(children_population[i].get_file(), names[i], children_population[i].get_chroms())
        new_population.append(child)
    return new_population


def sample_population(population: list, target_coverage: int, fragment_size: int, fragment_std: int):
    """
    This will create a fastq based on each member of the population.
    :param target_coverage: The target coverage value for the sample.
    :param population: a list of bacteria
    :return: None
    """
    for bacterium in population:
        bacterium.sample(target_coverage, fragment_size, fragment_std)


def extract_names(reference: str) -> list:
    """
    This function attempts to extract the chromosome names from a fasta file
    :param reference: The fasta file to analyze
    :return: A list of chromosome names
    """
    ref_names = []
    absolute_reference_path = pathlib.Path(reference)
    if absolute_reference_path.suffix == '.gz':
        with gzip.open(absolute_reference_path, 'rt') as ref:
            for line in ref:
                if line.startswith(">"):
                    ref_names.append(line[1:].rstrip())
    else:
        with open(absolute_reference_path, 'r') as ref:
            for line in ref:
                if line.startswith(">"):
                    ref_names.append(line[1:].rstrip())
    if not ref_names:
        print("Malformed fasta file. Missing properly formatted chromosome names.\n")
        sys.exit(1)

    return ref_names


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="Wrapper for gen_reads_runner.py that simulates multiple generations"
                                                 "of bacteria.")
    parser.add_argument('-r', type=str, required=True, metavar='reference.fasta',
                        help="Reference file for organism in fasta format")
    parser.add_argument('-g', type=int, required=True, metavar='generations', help="Number of generations to run")
    parser.add_argument('-i', type=int, required=True, metavar='initial pop', help="Initial population size")
    parser.add_argument('-k', type=float, required=False, metavar='cull pct',
                        help="Percentage of population to cull each cycle "
                             "(The default of 0.5 will keep population relatively stable)",
                        default=0.5)
    parser.add_argument('-c', type=int, required=False, default=10, metavar='coverage value',
                        help='Target coverage value for final set of sampled fastqs')
    parser.add_argument('--pe', nargs=2, type=int, required=False, metavar=('<int>', '<int>'), default=(500, 50),
                        help='Paired-end fragment length mean and std.')
    args = parser.parse_args()

    (ref_fasta, init_population_size, generations) = (args.r, args.i, args.g)
    cull_percentage = args.k
    # Need to add a numpy rng generator and redo how this is run.
    coverage = args.c
    (fragment_size, fragment_std) = args.pe

    chrom_names = extract_names(ref_fasta)

    population = initialize_population(ref_fasta, init_population_size, chrom_names)

    for i in range(generations):
        new_population = evolve_population(population, i+1)

        new_population = cull(new_population, cull_percentage)

        # If all elements get culled, then break the loop
        if not new_population:
            break

        population = new_population

    sample_population(population, coverage, fragment_size, fragment_std)


if __name__ == '__main__':
    start_time = time()
    main()
    print(f'Finished bacterial wrapper in {time() - start_time} seconds.')