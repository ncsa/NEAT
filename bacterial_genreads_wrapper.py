#!/usr/bin/env python

import genReads
import argparse
from Bio import SeqIO


def cull(population: list) -> list:
    """
    The purpose of this function will be to cull the bacteria created in the model
    :param population:
    :return:
    """
    return population


def crossover(population: list) -> list:
    """
    This function will take a list of individuals and do crossover type mixing
    :param population:
    :return new_population:
    """
    new_population = population
    return new_population


def initialize_population(reference: str, pop_size):
    """
    The purpose of this function is to evolve the bacteria
    :param reference:
    :param pop_size:
    :return population:
    """
    names = []
    for j in range(pop_size):
        names.append("name" + str(j))
    population = []
    for i in range(pop_size):
        args = ['-r', reference, '-R', '151', '-o', names[i], '--fa']
        new_member = genReads.main(args)
        population.append(new_member)
    return population


def evolve(new_population, param):
    pass


def main():
    parser = argparse.ArgumentParser(description='bacterial_genreads_wrapper.py')
    parser.add_argument('-r', type=str, required=True, metavar='/path/to/reference.fasta',
                        help="Reference file for organism in fasta format")
    parser.add_argument('-C', type=int, required=True, metavar='<str>', help="Number of cycles to run")
    parser.add_argument('-i', type=int, required=True, metavar='<str>', help="Initial population size")
    parser.add_argument('-c', type=float, required=False, metavar='<float>',
                        help="Percentage of population to cull each cycle (0.5 will keep population relatively stable)",
                        default=0.5)
    args = parser.parse_args()

    (ref_fasta, init_population_size, cycles) = (args.r, args.i, args.C)
    cull_percentage = args.c

    population = initialize_population(ref_fasta, init_population_size)
    for i in range(cycles):
        new_population = cull(population)
        # If all elements get culled, then break and quit
        if not new_population:
            break
        new_population = crossover(new_population)
        population = evolve(new_population, 2)
