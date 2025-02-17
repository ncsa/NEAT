---
title: 'Title'

tags:
  - Python
  - genomics
  - DNA sequencing
  - simulation

authors:
  - name: X Y. Z
    orcid: 0000-0000-0000-0000
    equal-contrib: true
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
  - name: Author Without ORCID
    equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
    affiliation: 2
  - name: Author with no affiliation
    corresponding: true # (This is how to denote the corresponding author)
    affiliation: 3
  - given-names: Ludwig
    dropping-particle: van
    surname: Beethoven
    affiliation: 3

affiliations:
 - name: Lyman Spitzer, Jr. Fellow, Princeton University, United States
   index: 1
   ror: 
 - name: Institution Name, Country
   index: 2
 - name: Independent Researcher, Country
   index: 3

date: XYZ February 2025

bibliography: paper.bib

# Summary

While the field of genomics has advanced significantly with the advent of high-throughput sequencing technologies, challenges related to the proprietary access, availability, complexity, and variability of this data can introduce difficulty in the development and validation of computational tools. Simulated datasets provide ground-truth estimates and reproducible scalability that are important in testing algorithms and benchmarking software. Ideally, these datasets mimic the properties of real sequencing datasets—from introducing specific patterns of sequencing errors to localized regions of mutations that often characterize cancer genomes.

# Statement of need

The NExt-generation sequencing Analysis Toolkit (NEAT) is an open-source Python package that creates simulated next-generation sequencing datasets. NEAT’s simulations account for a wide range of sequencing parameters (e.g., DNA read fragment length, sequencing error rates, mutation frequencies, etc.) and allow users to customize their sequencing data.<sup>1</sup> Since the original release of NEAT in 2016, most scripts have been greatly modified, and NEAT is currently on version 4.2. The code has undergone significant ongoing changes since 2020. Upgrading to Python 3 has enabled NEAT to achieve a flexible and intuitive user interface with minimal dependencies. The toolkit is optimized for both speed and accuracy, and new features have been implemented, such as improved ploidy simulation, mutation modeling, and the ability to model mutational profiles directly from data. A summary of algorithmic changes is provided in Table 1.

NEAT can integrate seamlessly with existing bioinformatics workflows, providing outputs in several common file formats. The toolkit’s ability to simulate gold-standard synthetic datasets with ground truth annotations is useful for testing bioinformatics pipelines. Uses of NEAT continue to be prominently featured—from scientists who have comprehensively sequenced the human Y chromosome<sup>2</sup> to researchers who use NEAT to evaluate and validate the performance of other high-profile bioinformatics tools.<sup>3, 4</sup> Earlier versions of NEAT have also demonstrated utility when benchmarked in comparison to similar tools.<sup>5</sup> The source code for both original and updated versions of NEAT is freely available on GitHub.<sup>1</sup>

# Tables

| Brief description of the problem, issue, and/or method | Describe the old function/algorithm as of 2.0? | Describe the new function/algorithm as of 4.X? | Why were changes made? | What evidence is there of improvement? |
|--------------------------------------------------------|------------------------------------------------|------------------------------------------------|-------------------------|-----------------------------------------|
| modeling quality score along read                     | Originally was also Markov-based, but it moved back to using a binning technique around 3.0. | The Markov process is new. | There were issues with runtime, not achieving the tapering effect on the read edges, and FastQC reports that were not ideal. | There were still long runtimes, but we need to make comparisons. We achieved the tapering effect, and the FastQC reports look good. |
| computing GC bias                                      | A script to calculate GC bias was written.     | N/A                                            | It caused user errors, and advances in sequencing technology make GC bias calculation less important. | There are no bugs related to this issue, and runtime decreased. |
| ploidy                                                | The script only really handled ploidy = 1 or 2. | We now handle potentially unbounded ploidy levels. | Organisms like plants are often polyploid, allowing for different kinds of genomes to be simulated. | Now, ploidy > 2 is properly handled. |
| variant insertion                                      | Had trouble with inserted variants, especially losing key information, such as genotype. Also, NEAT had a strict blacklist that prevented some otherwise valid variants from being inserted (in a way that felt arbitrary). | Now, NEAT preserves data from the insert variants VCF and puts that data into the final VCF. | Attempting to allow a greater number of variants and give users more control over the data that makes it into the final VCF. | Most variants are now inserted properly (long variants still an issue). |
| read generation                                        | Reads generated by NEAT showed a gap type structure, where stretches of the chromosome would have reads, then there would be a ~50bp gap, then a new group of reads would start. This was a consequence of how NEAT generated reads using a sliding window. | NEAT reads each chromosome and selects sets of coordinates for reads, then during read generation, it uses those coordinates to generate the final reads that go into the FASTQ. This eliminates the gap structure at the cost of some speed. | In an effort to make the final dataset look more realistic. | The gaps that once presented themselves in the reads are no longer present. |
| variant types                                          | A strict structure to the code and the variant types did not allow for new types of variants to be introduced. | Improved the code structure for variants, introduced a generic variant type that can read in variants from a VCF and replicate their genotype and some other information. Separated insertions and deletions, so they could be handled separately, allowing for more complex insertions and deletions. | Paving the way to add structural variants and copy number variants. | NEAT can now handle generic variants from insertions better and the code changes are in place to add more variant types. |
| custom rng                                            | NEAT used a complex RNG algorithm.             | Currently working on this!                      | NEAT does not need this complexity, and this change could reduce runtime and promote scientific reproducibility. | Currently working on this! |
| custom alignment                                      | Keep track of all changes in process.          | Isolated BAM file creation.                    | Every step of NEAT was tied to BAM creation, and there was no way to turn this setting “off,” which could increase runtimes. | Isolated BAM creation to prevent refactoring of the entire codebase. |


| Brief description of change | What was the old code’s function as of 2.0? | What was the new code’s function as of 4.X? | Why were changes made? | What evidence is there of improvement? |
|-----------------------------|---------------------------------------------|---------------------------------------------|-------------------------|-----------------------------------------|
| Testing framework           | There was no testing framework             | There are several automated tests in GitHub when pushing or pulling into specific branches | Convenience of development and testing. Covers most flags needed to rigorously test generate_reads. | Helped us spot random bugs and fix user issues (e.g., script that uses BED files). |
| Memory profiling            | There was no memory profiling              | We add a script for memory profiling in R that provides visualizations | We would like to profile the runtime and memory usage. | Around the same amount of memory taken for NEAT 2.0 vs. 4.X. |
| Rewritten codebase to be modular and pip installable using the poetry package | Could not be installable via package         | Can now be easily installed                    | User friendliness, convenience of development and testing, portability (less dependencies than before with higher value for production and usage). | Convenience. |
| Unit tests and code refactoring | Was not originally present               | Refactored with testable, discrete functions | Better to maintain the integrity of the code. | The code itself is testable and multiple people can collaborate on development. Readability. |
| User experience             | Command line interface used to use all flags | Refactored to use config files with a simplified structure | Config files are used as a design principle. | Readability, convenience, accuracy, debugging user tickets more easily, scientific reproducibility. |
| Parallelization             | Was not originally present                 | In progress!                                  | Speed up runtimes now that we have refactored code. Need to “split up” BAM file reading, etc., in a way that can be parallelized. | N/A |

# Acknowledgements

# References
1.	Stephens ZD, Hudson ME, Mainzer LS, Taschuk M, Weber MR, Iyer RK. Simulating Next-Generation Sequencing Datasets from Empirical Mutation and Sequencing Models. PLOS ONE. 2016;11(11). doi:10.1371/journal.pone.0167047
2.	Rhie A, Nurk S, Cechova M, et al. The complete sequence of a human Y chromosome. Nature. 2023;621(7978):344-354. doi:10.1038/s41586-023-06457-y
3.	Lefouili M, Nam K. The evaluation of Bcftools mpileup and GATK HaplotypeCaller for variant calling in non-human species. Sci Rep. 2022;12(1). doi:10.1038/s41598-022-15563-2
4.	Zhao S, Agafonov O, Azab A, Stokowy T, Hovig E. Accuracy and efficiency of germline variant calling pipelines for human genome data. Scientific Reports. 2020;10(1). doi:10.1038/s41598-020-77218-4
5.	Alosaimi S, Bandiang A, van Biljon N, et al. A broad survey of DNA sequence data simulation tools. Brief Funct Genomics. 2020;19(1):49-59. doi:10.1093/bfgp/elz033
