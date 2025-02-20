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
 - name: National Center for Supercomputing Applications, Genomics Group, Urbana, IL, USA, 61801, 
   index: 1
 - name: University of Illinois at Chicago, Chicago, IL, USA, 60607
   index: 2

date: 20 February 2025

bibliography: paper.bib
---

# Summary

While the field of genomics has advanced significantly with the advent of high-throughput sequencing technologies, challenges related to the availability, complexity, and variability of this data can introduce difficulty in the development and validation of computational tools. Simulated deoxyribonucleic acid (DNA) sequencing datasets provide ground-truth estimates and reproducible scalability that are important in testing algorithms and benchmarking software. Ideally, these datasets mimic the properties of real sequencing datasets—from introducing specific patterns of sequencing errors to modeling localized regions of mutations.

# Statement of need

The NExt-generation sequencing Analysis Toolkit (NEAT) is an open-source Python package that creates simulated next-generation sequencing datasets. NEAT’s simulations account for a wide range of sequencing parameters (e.g., DNA read fragment length, sequencing error rates, mutation frequencies, etc.) and allow users to customize their sequencing data.[@Stephens:2016] Since the original release of NEAT in 2016, most scripts have been greatly modified, and NEAT is currently on version 4.2. The code has undergone significant ongoing changes since 2020. Upgrading to Python 3 has enabled NEAT to achieve a flexible and intuitive user interface with minimal dependencies. The toolkit is optimized for both speed and accuracy, and new features have been implemented, such as improved ploidy simulation, mutation modeling, and the ability to model mutational profiles directly from data. A summary of algorithmic changes is provided in Table 1.

NEAT can integrate seamlessly with existing bioinformatics workflows, providing outputs in several common file formats. The toolkit’s ability to simulate gold-standard synthetic datasets with ground truth annotations is useful for testing bioinformatics pipelines. Uses of NEAT continue to be prominently featured—from scientists who have comprehensively sequenced the human Y chromosome<sup>2</sup> to researchers who use NEAT to evaluate and validate the performance of other high-profile bioinformatics tools.<sup>3, 4</sup> Earlier versions of NEAT have also demonstrated utility when benchmarked in comparison to similar tools.<sup>5</sup> The source code for both original and updated versions of NEAT is freely available on GitHub.<sup>1</sup>

\newpage

# Tables

### Table 1. Algorithmic Improvements and Methodological Changes

| # | Feature Name | Prior Implementation (v2.0) | Updated Implementation (v4.X) |
|---|-------------|------------------------------|--------------------------------|
| 1 | **BAM File Generation** | File generation was tightly integrated with all NEAT processes | BAM creation was isolated from core functions |
|---|-------------|------------------------------|--------------------------------|
| 2 | **GC Bias Computation** | Used a custom script for GC bias calculation | Feature deprecated |
|---|-------------|------------------------------|--------------------------------|
| 3 | **Ploidy Simulation** | Limited to diploid organisms in practice | Supports unbounded ploidy levels |
|---|-------------|------------------------------|--------------------------------|
| 4 | **Read Generation** | The sliding-window approach to generate reads resulted in artificial gaps in sequencing reads (~50 base pairs) | A new form of coordinate-based read selection eliminates these gaps |
|---|-------------|------------------------------|--------------------------------|
| 5 | **Read Quality Modeling** | Markov-based model | Binning method with an option to also implement a revised Markov-based model |
|---|-------------|------------------------------|--------------------------------|
| 6 | **Variant Insertion** | Issues with inserted variants (loss of genotype data, prevented certain valid variants from insertion) | Preserves genotype data in the final variant call format (VCF) file |
|---|-------------|------------------------------|--------------------------------|
| 7 | **Variant Handling** | The code structure limited the introduction of new variant types | A modular design supports generic variant handling and the separation of insertions and deletions |

The creation of simulated Binary Alignment Map (BAM) files (**1**) in NEAT 2.0 was tightly integrated with all NEAT functions. The new update isolates BAM creation, improving runtime and modularity. Guanine-cytosine (GC) bias computation (**2**) was removed due to redundancy, and its removal reduced runtime. Ploidy simulation (**3**) has been extended to improve accurate simulation of tumor genomes and polyploid organisms (e.g., plants), and ploidy inputs greater than two and fractional ploidies are now handled. Previously, NEAT 2.0's read generation (**4**) algorithm introduced read gaps (~50 base pairs) due to its sliding-window approach. The updated coordinate-based selection eliminates these gaps. Modeling of sequencing quality scores for each nucleotide base (**5**) was updated by incorporating a revised Markov model alongside a binning method. We accurately account for a tapering effect that reduces sequencing quality scores along a simulated sequence's edges. Variant insertion (**6**) was updated to preserve genotype data in the final simulated Variant Call Format (VCF) file, improving accuracy and giving users greater control over the insertion of variants. Finally, variant handling (**7**) has been modularized to support structural and copy number variants, increasing flexibility and ensuring future extensibility for handling more complex variants.

\newpage

### Table 2. Improvements in User Experience

| # | Feature Name | Prior Implementation (v2.0) | Updated Implementation (v4.X) |
|---|-------------|------------------------------|--------------------------------|
| 1 | **Automated Testing** | No formal testing framework | Implemented continuous integration with GitHub-based automated tests |
|---|-------------|------------------------------|--------------------------------|
| 2 | **Configuration Files** | Required explicit command-line flags | Introduced structured configuration files |
|---|-------------|------------------------------|--------------------------------|
| 3 | **Friendly Installation** | Not installable as a package | Fully modular and pip-installable via Poetry |
|---|-------------|------------------------------|--------------------------------|
| 4 | **Refactored Unit Testing** | Not originally present | Rewritten with testable, discrete functions |

Our new continuous integration pipeline (**1**) detects bugs early, streamlining development and enhancing error detection (e.g., handling of multiple genomic file formats as inputs and outputs). Configuration files in NEAT v4.X (**2**) and package installation (**3**) facilitate user friendliness and portability. NEAT v4.X features testable, discrete functions (**4**) that allows users to debug more easily. Parallelization of NEAT v4.X is in progress.

\newpage

# Acknowledgements

We thank the original creators of NEAT: Zachary D. Stephens, Matthew E. Hudson, Liudmila S. Mainzer, Morgan Taschuk, Matthew R. Weber, and Ravishankar K. Iyer. 

We also thank Raghid Alhamzy, Yash Wasnik, Varenya Jain, and Karen H. Xiong.

# References
1.	Stephens ZD, Hudson ME, Mainzer LS, Taschuk M, Weber MR, Iyer RK. Simulating Next-Generation Sequencing Datasets from Empirical Mutation and Sequencing Models. PLOS ONE. 2016;11(11). doi:10.1371/journal.pone.0167047
2.	Rhie A, Nurk S, Cechova M, et al. The complete sequence of a human Y chromosome. Nature. 2023;621(7978):344-354. doi:10.1038/s41586-023-06457-y
3.	Lefouili M, Nam K. The evaluation of Bcftools mpileup and GATK HaplotypeCaller for variant calling in non-human species. Sci Rep. 2022;12(1). doi:10.1038/s41598-022-15563-2
4.	Zhao S, Agafonov O, Azab A, Stokowy T, Hovig E. Accuracy and efficiency of germline variant calling pipelines for human genome data. Scientific Reports. 2020;10(1). doi:10.1038/s41598-020-77218-4
5.	Alosaimi S, Bandiang A, van Biljon N, et al. A broad survey of DNA sequence data simulation tools. Brief Funct Genomics. 2020;19(1):49-59. doi:10.1093/bfgp/elz033
