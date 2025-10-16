import numpy as np
from numpy.random import default_rng

from neat.common.ploid_functions import pick_ploids, get_genotype_string


def test_pick_ploids_homozygous_all_copies():
    rng = default_rng(42)
    ploidy = 4
    geno = pick_ploids(ploidy=ploidy, homozygous_frequency=1.0, number_alts=1, rng=rng)
    # All copies should be non-zero (ALT=1)
    assert geno.shape == (ploidy,)
    assert np.all(geno > 0), f"Expected all >0, got {geno}"


def test_pick_ploids_heterozygous_polyploid_half_nonzero():
    rng = default_rng(123)
    ploidy = 4
    geno = pick_ploids(ploidy=ploidy, homozygous_frequency=0.0, number_alts=2, rng=rng)
    # Expect exactly ploidy // 2 non-zero entries according to implementation
    assert int(np.count_nonzero(geno)) == ploidy // 2
    # And all non-zeros should be in 1..number_alts
    assert set(np.unique(geno[geno > 0])).issubset({1.0, 2.0})


def test_pick_ploids_haploid_only_one_alt():
    rng = default_rng(7)
    ploidy = 1
    geno = pick_ploids(ploidy=ploidy, homozygous_frequency=0.0, number_alts=3, rng=rng)
    # Haploid always exactly one allele (ALT 1..3)
    assert geno.shape == (1,)
    assert geno[0] in {1.0, 2.0, 3.0}


def test_get_genotype_string():
    arr = np.array([0, 1, 2])
    s = get_genotype_string(arr)
    assert s == "0/1/2"