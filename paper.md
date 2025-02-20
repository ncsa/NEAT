---
title: 'Enhancing Next-Generation Sequencing Simulation: Updates to NEAT'

tags:
  - Python
  - genomics
  - DNA sequencing
  - simulation

authors:
  - name: Joshua M. Allen
    equal-contrib: true
    affiliation: 1
  - name: Keshav R. Gandhi
    orcid: 0009-0000-1718-1862
    equal-contrib: true
    corresponding: true
    affiliation: "1, 2"
  - name: Christina E. Fliege
    orcid: 0000-0001-8085-779X
    affiliation: 1

affiliations:
 - name: National Center for Supercomputing Applications, Genomics Group, 
   index: 1
 - name: University of Illinois at Chicago, 
   index: 2

date: XYZ February 2025

bibliography: paper.bib

# Summary

While the field of genomics has advanced significantly with the advent of high-throughput sequencing technologies, challenges related to the proprietary access, availability, complexity, and variability of this data can introduce difficulty in the development and validation of computational tools. Simulated datasets provide ground-truth estimates and reproducible scalability that are important in testing algorithms and benchmarking software. Ideally, these datasets mimic the properties of real sequencing datasets—from introducing specific patterns of sequencing errors to localized regions of mutations that often characterize cancer genomes.

# Statement of need

The NExt-generation sequencing Analysis Toolkit (NEAT) is an open-source Python package that creates simulated next-generation sequencing datasets. NEAT’s simulations account for a wide range of sequencing parameters (e.g., DNA read fragment length, sequencing error rates, mutation frequencies, etc.) and allow users to customize their sequencing data.<sup>1</sup> Since the original release of NEAT in 2016, most scripts have been greatly modified, and NEAT is currently on version 4.2. The code has undergone significant ongoing changes since 2020. Upgrading to Python 3 has enabled NEAT to achieve a flexible and intuitive user interface with minimal dependencies. The toolkit is optimized for both speed and accuracy, and new features have been implemented, such as improved ploidy simulation, mutation modeling, and the ability to model mutational profiles directly from data. A summary of algorithmic changes is provided in Table 1.

NEAT can integrate seamlessly with existing bioinformatics workflows, providing outputs in several common file formats. The toolkit’s ability to simulate gold-standard synthetic datasets with ground truth annotations is useful for testing bioinformatics pipelines. Uses of NEAT continue to be prominently featured—from scientists who have comprehensively sequenced the human Y chromosome<sup>2</sup> to researchers who use NEAT to evaluate and validate the performance of other high-profile bioinformatics tools.<sup>3, 4</sup> Earlier versions of NEAT have also demonstrated utility when benchmarked in comparison to similar tools.<sup>5</sup> The source code for both original and updated versions of NEAT is freely available on GitHub.<sup>1</sup>

# Tables

## Algorithmic Improvements and Methodological Changes

## Table 1. Enhancements in Algorithmic Performance

| Feature | Prior Implementation (v2.0) | Updated Implementation (v4.X) | Rationale for Change | Demonstrated Improvement |
|---------|------------------------------|--------------------------------|----------------------|-------------------------|
| **Read Quality Modeling** | Markov-based model | Binning method with an option to also implement a revised Markov-based model | Did not achieve a tapering effect on a simulated read's edges | The tapering effect was achieved with the revised Markov model |
| **Guanine-Cytosine (GC) Bias Computation** | Used a custom script for GC bias calculation | Feature deprecated | The script has decreased relevance due to advances in sequencing technology | Reduced runtime and absence of reported bugs |
| **Ploidy Simulation** | Limited to diploid organisms in practice | Supports unbounded ploidy levels | This is essential for the accurate modeling of tumors or polyploid organisms, such as plants | Inputs of ploidy greater than two and fractional ploidies will correctly simulate reads |
| **Variant Insertion** | Issues with inserted variants (loss of genotype data, prevented certain valid variants from insertion) | Preserves genotype data in the final variant call format (VCF) file | Allows greater flexibility and user control over variant inclusion | Improved accuracy of inserted variants (long variant support still in progress) |
| **Read Generation** | The sliding-window approach to generate reads resulted in artificial gaps in sequencing reads (~50 base pairs) | A new form of coordinate-based read selection eliminates these gaps | Aimed to produce datasets more representative of real sequencing patterns | Elimination of artificial gaps |
| **Variant Type Handling** | The code structure limited the introduction of new variant types | A modular design supports generic variant handling and the separation of insertions and deletions | Paves the way for structural and copy number variant support | More flexible insertion handling and future extensibility |
| **Binary Alignment Map (BAM) File Generation** | File generation was tightly integrated with all NEAT processes | BAM creation was isolated from core functions | Improves runtime and modularity | BAM generation can now be toggled independently |

## Table 2. Performance Enhancements and User-Centric Modifications

| Feature | Prior Implementation (v2.0) | Updated Implementation (v4.X) | Rationale for Change | Demonstrated Improvement |
|---------|------------------------------|--------------------------------|----------------------|-------------------------|
| **Modular Codebase & Installation** | Not installable as a package | Fully modular and pip-installable via Poetry | Facilitates ease of development, portability, and deployment | Reduced dependencies, improved maintainability |
| **Code Refactoring & Unit Testing** | Monolithic, unstructured codebase | Rewritten with testable, discrete functions | Enhances maintainability and collaborative development | Improved code readability and integrity |
| **User Experience: Configuration Management** | Required explicit command-line flags | Introduced structured configuration files | Improves usability, debugging, and reproducibility | Simplified interface, increased accessibility |
| **Automated Testing Framework** | No formal testing framework | Implemented continuous integration with GitHub-based automated tests | Improves development efficiency and debugging capabilities | Enhanced detection of random bugs and user issues (e.g., file handling) |

Parallelization and memory profiling tools will be updated shortly.

# Acknowledgements

We thank the original creators of NEAT. We also thank Raghid, Yash, Varenya, and Karen.

# References
1.	Stephens ZD, Hudson ME, Mainzer LS, Taschuk M, Weber MR, Iyer RK. Simulating Next-Generation Sequencing Datasets from Empirical Mutation and Sequencing Models. PLOS ONE. 2016;11(11). doi:10.1371/journal.pone.0167047
2.	Rhie A, Nurk S, Cechova M, et al. The complete sequence of a human Y chromosome. Nature. 2023;621(7978):344-354. doi:10.1038/s41586-023-06457-y
3.	Lefouili M, Nam K. The evaluation of Bcftools mpileup and GATK HaplotypeCaller for variant calling in non-human species. Sci Rep. 2022;12(1). doi:10.1038/s41598-022-15563-2
4.	Zhao S, Agafonov O, Azab A, Stokowy T, Hovig E. Accuracy and efficiency of germline variant calling pipelines for human genome data. Scientific Reports. 2020;10(1). doi:10.1038/s41598-020-77218-4
5.	Alosaimi S, Bandiang A, van Biljon N, et al. A broad survey of DNA sequence data simulation tools. Brief Funct Genomics. 2020;19(1):49-59. doi:10.1093/bfgp/elz033
