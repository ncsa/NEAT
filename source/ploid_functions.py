from numpy.random import Generator


def pick_ploids(rand_num_gen: Generator, ploidy, homozygous_frequency, number_alts=1) -> list:
    """
    Applies alts to random ploids. Picks at least one, maybe more.
    :param ploidy: how many copies of each chromosome this organism has
    :param homozygous_frequency:
    :param number_alts: If there is more than one alt, this will assign a random alt to the ploids
    :return: a list of strings representing the genotype of each ploid.
    """
    # number of ploids to make this mutation on (always at least 1)
    if rand_num_gen.random() < homozygous_frequency:
        if ploidy <= 2:
            # If it's homozygous. it's on both ploids
            how_many = ploidy
        else:
            # if it's polyploid, we'll just pick some.
            # TODO may need to improve the modeling for polyploid
            how_many = 1
            for i in range(ploidy):
                # Not totally sure how to model this, so I'm counting each
                # ploid as a separate homozygous event. That doesn't exactly make
                # sense though so we'll improve this later.
                if rand_num_gen.random() < homozygous_frequency:
                    how_many += 1
                else:
                    break
    else:
        how_many = 1

    # wp is just the temporary genotype list
    wp = [0] * ploidy
    while how_many > 0:
        x = rand_num_gen.choice(range(ploidy))
        # pick a random alt. in VCF terminology, 0 = REF, 1 = ALT1, 2 = ALT2, etc
        wp[x] = rand_num_gen.choice(range(1, number_alts + 1))
        how_many -= 1

    return wp


def which_ploid(genotype_list):
    """
    Finds the ploid given an input genotype in list form. For example a ploidy of 4 might look like [0, 1, 0, 0]
    :param genotype_list: a list of the genotype for each ploid in the dataset (0 indicates ref, 1 first alt,
                          2 second alt, etc.). Elements should be ints
    """
    return [genotype_list.index(x) for x in genotype_list if x == 0]
