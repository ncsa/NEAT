---
title: 'Enhancing short-read sequencing simulation: Updates to NEAT'

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
    email: krg3@uic.edu
    equal-contrib: true
    corresponding: true
    affiliation: "1, 2"
  - name: Raghid Alhamzy
    affiliation: 1
  - name: Yash Wasnik
    affiliation: 3    
  - name: Christina E. Fliege
    affiliation: 1
    email: cfliege2@illinois.edu
    corresponding: true

affiliations:
 - name: National Center for Supercomputing Applications, Genomics Group, Urbana, IL, USA, 61801, 
   index: 1
   
 - name: University of Illinois at Chicago, Chicago, IL, USA, 60607,
   index: 2

 - name: Blue Health Intelligence, Chicago, United States, 60611
   index: 3
   
date: 25 March 2026

bibliography: paper.bib
---

# Summary

While the field of genomics has advanced significantly with the advent of high-throughput sequencing technologies, challenges related to the availability, complexity, and variability of this data can introduce difficulty to the development and validation of computational tools. Simulated short-read sequencing datasets provide researchers a way to obtain reproducible, verified data to test algorithms and benchmark software. Simulations also avoid the limitations of working with real data, including the cost of genomic sequencing, time to process sequencing data, accessibility of data, and protection of privacy for human datasets. Ideally, simulated datasets mimic the properties of real sequencing datasets. From introducing specific patterns of sequencing errors to modeling localized regions of mutations, realistic datasets are necessary to evaluate the accuracy and robustness of downstream alignment, variant-calling, and analysis pipelines. Sequencing parameters such as guanine-cytosine (GC) content, repeat structure, and genome complexity vary across species, and it is important for simulated reads to reflect these species-specific characteristics. The NExt-generation sequencing Analysis Toolkit (NEAT) is an open-source Python package that creates such simulated sequencing datasets and was originally released as NEAT 2.0 in 2016. Here, we describe ongoing updates released collectively as NEAT 4.3, including new features and increased speed, accuracy, and usability.

# Statement of need

Developing and validating methods for read alignment, variant calling, structural variant detection, and other analyses depends on access to genomic data with known ground truth. As a result, many research groups rely on simulated reads, and it can be useful to vary the parameters of the sequencing process itself in order to demonstrate how analysis pipelines may vary in performance across different coverages, read lengths, ploidies, and more. NEAT addresses this need as an open-source Python package that can be easily integrated within existing bioinformatics workflows—its simulations account for a wide range of sequencing parameters (e.g., coverage, fragment length, sequencing error models, mutational frequencies, ploidy, etc.) and allow users to customize their sequencing data [@Stephens:2016].

NEAT is designed to simulate short reads. NEAT is adaptable to different sequencing platforms with user-specific, custom error models and the capability to handle single-base substitutions, insertions, deletions, and larger structural variants. Unlike simulators that rely solely on fixed error profiles, NEAT can learn empirical mutation and sequencing models from real datasets and use these models to generate realistic sequencing data, providing outputs in several common bioinformatics file formats, including FASTQ, Binary Alignment/Map (BAM), and Variant Call Format (VCF) files. NEAT can generate BAM and VCF outputs, facilitating benchmarking of full bioinformatics pipelines.

# Software design

After the original release of NEAT 2.0, most scripts have undergone significant changes. Upgrading to Python 3 brought NEAT 4.3 up to modern coding standards and allowed it to use standard Python libraries in order to streamline the code and improve its maintainability. The toolkit is optimized for both speed and accuracy, and new features have been implemented. A summary of algorithmic changes is provided in **Table 1**, and a summary of updates to NEAT’s software robustness and user experience is provided in **Table 2**.

### Table 1. Algorithmic Improvements and Methodological Changes

+---+-----------------------+----------------------------------------------------------------+------------------------------------------------------+
| # | Feature Name          | Prior Implementation (2.0)                                     | Updated Implementation (4.3)                         |
+===+=======================+================================================================+======================================================+
| 1 | BAM File Generation   | File generation was tightly integrated with all NEAT processes | BAM creation isolated from core functions            |
+---+-----------------------+----------------------------------------------------------------+------------------------------------------------------+
|   |                       |                                                                |                                                      |
+---+-----------------------+----------------------------------------------------------------+------------------------------------------------------+
| 2 | GC Bias Computation   | Used a custom script for GC bias calculation                   | Removed during refactoring                           |
+---+-----------------------+----------------------------------------------------------------+------------------------------------------------------+
|   |                       |                                                                |                                                      |
+---+-----------------------+----------------------------------------------------------------+------------------------------------------------------+
| 3 | Ploidy Simulation     | Limited to diploid organisms in practice                       | Supports polyploidy                                  |
+---+-----------------------+----------------------------------------------------------------+------------------------------------------------------+
|   |                       |                                                                |                                                      |
+---+-----------------------+----------------------------------------------------------------+------------------------------------------------------+
| 4 | Read Generation       | Sliding-window approach to generate reads                      | A new form of coordinate-based read selection        |
+---+-----------------------+----------------------------------------------------------------+------------------------------------------------------+
|   |                       |                                                                |                                                      |
+---+-----------------------+----------------------------------------------------------------+------------------------------------------------------+
| 5 | Read Quality Modeling | Markov-based model                                             | Revised binning method (optional Markov-based model) |
+---+-----------------------+----------------------------------------------------------------+------------------------------------------------------+
|   |                       |                                                                |                                                      |
+---+-----------------------+----------------------------------------------------------------+------------------------------------------------------+
| 6 | Variant Handling      | Limited introduction of new variant types                      | A modular design with generic variant handling       |
+---+-----------------------+----------------------------------------------------------------+------------------------------------------------------+
|   |                       |                                                                |                                                      |
+---+-----------------------+----------------------------------------------------------------+------------------------------------------------------+
| 7 | Variant Insertion     | Issues with inserted variants (loss of genotype data)          | Preserves genotype data in the final VCF file        |
+---+-----------------------+----------------------------------------------------------------+------------------------------------------------------+

Below, we summarize methodological updates (**Table 1**)  present in NEAT 4.3:

- In NEAT, simulated reads can be accompanied by ground-truth alignment and variant outputs to support benchmarking pipelines. BAM file generation (**1**) produces aligned reads in a standard format for mapping and downstream processing. In NEAT 2.0, BAM creation was tightly integrated with all functions, which increased runtimes for users who only simulated FASTQ files. This update isolates BAM creation from the core simulation, improving modularity when alignment is not needed.


- GC bias (**2**) is the tendency of sequencing data to have uneven coverage as a function of a local genomic region’s GC content, often underrepresenting regions of particularly high or low GC levels [@Benjamini:2012; @Ross:2013]. NEAT 4.3 removes NEAT 2.0’s module that models GC-dependent coverage bias. Evidence suggests that much of GC bias arises from library preparation and amplification steps rather than sequencing itself [@Benjamini:2012]. Additionally, recent benchmarks across seven sequencing platforms reported relatively uniform coverage across the moderate GC range (20% to 60%), which comprises about 95% of the human genome [@Kim:2021].


- Ploidy simulation (**3**) mediates how many chromosomal copies are present in simulated sequencing data. It has been extended to improve the simulation of tumor genomes and polyploid organisms, such as plants, as ploidy inputs greater than two are now supported.


- Read generation (**4**) determines how read start positions are selected across the genome to simulate realistic coverage. Previously, NEAT 2.0's read generation algorithm introduced gaps (~50 base pairs) due to its sliding-window approach. The updated coordinate-based selection eliminates these gaps.


- Read quality modeling (**5**) assigns quality scores to each base, impacting downstream error profiles and variant calling. A revised binning method was implemented (alongside an optional Markov model that uses per-position transition matrices) to account for a slight tapering effect that reduces sequencing quality scores at a simulated sequence's edges.


- Variant handling (**6**) defines how different mutation classes are represented. This logic has been modularized to separate insertions from deletions and facilitate the future handling of structural and copy-number variants.


- Variant insertion (**7**) handles how user-provided variants are incorporated into simulated reads and reported in the final VCF output. NEAT allows users to input their own VCF file of variants to include in the simulated reads. In the final simulated VCF, NEAT 4.3 preserves per-sample genotypes and allele-frequency information.


### Table 2. Improvements to Software Robustness and User Experience

+---+-----------------------+----------------------------------------------------------------+------------------------------------------------------+
| # | Feature Name          | Prior Implementation (2.0)                                     | Updated Implementation (4.3)                         |
+===+=======================+================================================================+======================================================+
| 1 | Automated Testing     | No formal testing framework                                    | Implemented continuous integration with GitHub-based |
|   |                       |                                                                | automated tests                                      |
+---+-----------------------+----------------------------------------------------------------+------------------------------------------------------+
|   |                       |                                                                |                                                      |
+---+-----------------------+----------------------------------------------------------------+------------------------------------------------------+
| 2 | Configuration Files   | Required explicit command-line flags                           | Introduced structured configuration files            |
+---+-----------------------+----------------------------------------------------------------+------------------------------------------------------+
|   |                       |                                                                |                                                      |
+---+-----------------------+----------------------------------------------------------------+------------------------------------------------------+
| 3 | Detailed Logging      | Minimal error logging                                          | Extensive logs to recreate runs and provide          |
|   |                       |                                                                | descriptions of errors                               |
+---+-----------------------+----------------------------------------------------------------+------------------------------------------------------+
|   |                       |                                                                |                                                      |
+---+-----------------------+----------------------------------------------------------------+------------------------------------------------------+
| 4 | Friendly Installation | Not installable as a package                                   | Easy installation via Bioconda or Poetry             |
+---+-----------------------+----------------------------------------------------------------+------------------------------------------------------+
|   |                       |                                                                |                                                      |
+---+-----------------------+----------------------------------------------------------------+------------------------------------------------------+
| 5 | Parallelization       | Single-threaded only                                           | Multi-threaded runs are now possible                 |
+---+-----------------------+----------------------------------------------------------------+------------------------------------------------------+
|   |                       |                                                                |                                                      |
+---+-----------------------+----------------------------------------------------------------+------------------------------------------------------+
| 6 | Refactored Unit       | Not originally present                                         | Rewritten with testable, discrete functions          |
|   | Testing               |                                                                |                                                      |
+---+-----------------------+----------------------------------------------------------------+------------------------------------------------------+

**Table 2** outlines improvements to the software robustness and user experience. Our new continuous integration pipeline (**1**) detects bugs early, streamlining development and enhancing error detection (e.g., handling of multiple genomic file formats as inputs and outputs). Configuration files in NEAT 4.3 (**2**), detailed error logging (**3**), and package installation (**4**) facilitate user friendliness and portability. NEAT 4.3's read simulator is also parallelized (**5**), promoting faster runtimes and ease of use. NEAT 4.3 features testable, discrete functions (**6**) that allow users to debug more easily. Documentation of these improvements, alongside benchmarks of parallelization-related performance speedups, is maintained as part of the project’s online repository and updated regularly. Quality-of-life development continues into the present.

# State of the field

A large number of bulk sequencing read simulators exist. Early comparative reviews cataloged over twenty simulators of genomic data and highlighted trade-offs between realism, speed, and flexibility across sequencing platforms and use cases such as whole genomes, exomes, and metagenomes [@Escalona:2016; @Zhao:2017]. More recent benchmarks of read simulators have focused on systematic performance comparisons and whether empirical error-profile learning improves realism [@Alosaimi:2020; @Milhaven:2023; @Schmeing:2021]. Among short-read simulators benchmarked in these studies, commonly used options in addition to NEAT include ART, CuReSim, DWGSIM, GemSIM, InSilicoSeq, Mason, pIRS, ReSeq, SInC, and wgsim [@Huang:2012; @Caboche:2014; @Homer:2010; @McElroy:2012; @Gourle:2019; @Holtgrewe:2010; @Hu:2012; @Schmeing:2021; @Pattnaik:2014; @Li:2011]. However, only a subset of existing read simulators directly produce short-read sequencing data with explicit ground truth suitable for benchmarking alignment and variant-calling pipelines [@Alosaimi:2020]. Milhaven and Pfeifer (2023) provide a detailed feature matrix comparing six short-read simulators (ART, DWGSIM, InSilicoSeq, Mason, NEAT, and wgsim), noting NEAT for its ability to generate such benchmarking-ready ground-truth data while also documenting each tool’s supported platforms, error and quality-score modeling, mutation types, and available outputs.

Additionally, in their benchmark of twenty DNA read simulators that use reference genomes and produce FASTQ files, Alosaimi et al. (2020) found that NEAT has demonstrated utility in comparison to similar tools; although NEAT’s runtimes were the second-longest, it achieved the second-highest mapping sensitivity and precision on a human chromosome 22 test set [@Alosaimi:2020]. Among these tools, NEAT is noteworthy for its ability to combine mutation and sequencing models in a single framework and accept user-specified mutation models [@Stephens:2016]. More recently, all six read simulators benchmarked by Milhaven and Pfeifer (2023) were found to have strengths that can be selected based on a user’s desired set of features, noting NEAT 3.0 in particular for its realism in sequencing—but also its slow simulation runtimes [@Milhaven:2023]. NEAT 4.3 addresses weaknesses mentioned in these reviews by implementing algorithmic updates and a parallelized design that enables faster, multi-threaded simulation runs relative to NEAT 3.0.

The rise of single-cell DNA sequencing has also impacted read simulators. Single-cell sequencing may introduce additional sources of variability (e.g., amplification bias, allelic dropout, and coverage heterogeneity) that most bulk simulators do not model explicitly. These tools may wrap existing bulk simulators to generate correlated bulk and single-cell genomes. Tools such as SCSIM simulate correlated bulk and single-cell samples on top of underlying variant processes and could complement NEAT’s bulk-focused design [@Giguere:2020]. It is also possible to integrate NEAT 4.3 into workflows that produce single-cell sequencing datasets.

# Research impact statement

Within this landscape, NEAT is notable for integrating ground-truth BAM and VCF outputs suitable for systematic benchmarking pipelines and combining user-trainable mutation and sequencing models, supporting primarily short-read error profiles via custom models [@Alosaimi:2020; @Schmeing:2021; @Stephens:2016; @Milhaven:2023]. Uses of NEAT continue to be featured—from assisting benchmarks in the sequencing of the human Y chromosome [@Rhie:2023] to evaluating the performance of other bioinformatics tools [@Lefouili:2022; @Zhao:2020]. NEAT has supported benchmarking and method development across a range of bioinformatics tasks, including optimization of high-throughput variant-calling workflows [@Ahmed:2019; @Kendig:2019], feasibility studies of exome sequencing [@RuizSchultz:2021], pan-genome mapping [@Jandrasits:2019], Bayesian approaches for resolving ambiguously mapped reads [@Shah:2021], and ultra-sensitive multi-sample variant calling [@Delhomme:2020].

The updates described here (NEAT 4.3) focus on improving NEAT’s capabilities so that it remains a practical tool for future work in the field alongside other simulators. The source code for both original and updated versions of NEAT is freely available on GitHub [@Stephens:2016].

# Acknowledgements

We thank the original creators of NEAT: Zachary D. Stephens, Matthew E. Hudson, Liudmila S. Mainzer, Morgan Taschuk, Matthew R. Weber, and Ravishankar K. Iyer. 

We also thank Varenya Jain, Meredith Pudlewski, and Karen H. Xiong for their work on updating NEAT as well as all of the NEAT users who have contributed to the codebase.

Finally, we thank the National Center for Supercomputing Applications’ Students Pushing Innovation (SPIN) program for funding a portion of this research. This research was supported in part by the Illinois Computes project through the University of Illinois Urbana-Champaign.

# References
