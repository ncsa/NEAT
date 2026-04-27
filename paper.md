---
title: 'Enhancing short-read sequencing simulation: Updates to NEAT'

tags:
  - Python
  - genomics
  - DNA sequencing
  - simulation

authors:
  - name: Joshua M. Allen
    orcid: 0009-0008-0002-5239
    equal-contrib: true
    affiliation: 1
  - name: Keshav R. Gandhi
    orcid: 0009-0000-1718-1862
    email: krg3@uic.edu
    equal-contrib: true
    corresponding: true
    affiliation: "1, 2"
  - name: Raghid Alhamzy
    orcid: 0009-0005-6614-7050
    affiliation: 1
  - name: Yash Wasnik
    orcid: 0009-0006-0108-7445
    affiliation: 3    
  - name: Christina E. Fliege
    orcid: 0000-0001-8085-779X
    affiliation: 1
    email: cfliege2@illinois.edu
    corresponding: true

affiliations:
 - name: National Center for Supercomputing Applications, Genomics Group, Urbana, IL, United States, 61801, 
   index: 1
   
 - name: University of Illinois at Chicago, Chicago, IL, United States, 60607,
   index: 2

 - name: Blue Health Intelligence, Chicago, IL, United States, 60611
   index: 3
   
date: 15 April 2026

bibliography: paper.bib
---

# Summary

While the field of genomics has advanced significantly with high-throughput sequencing technologies, challenges related to the availability, complexity, and variability of data introduce difficulty when developing and validating computational tools. Simulated short-read sequencing datasets provide researchers with reproducible, verified data to test algorithms and benchmark software. Simulations also avoid the limitations of working with real data, such as the cost of genomic sequencing, processing time, accessibility, and protection of privacy for human datasets. Ideally, simulated datasets mimic the properties of real sequencing data, which is necessary to evaluate the accuracy and robustness of downstream alignment, variant-calling, and analysis pipelines. Sequencing parameters such as ploidy, repeat structure, and genome complexity vary across species, and it is important for simulated reads to reflect these species-specific characteristics. The NExt-generation sequencing Analysis Toolkit (NEAT) is an open-source Python package that creates simulated sequencing datasets, originally released in 2016. The current release, v4.4, introduces increased speed, accuracy, and usability.

# Statement of need

Developing and validating methods for read alignment, variant calling, and other analyses requires genomic data with known ground truth. Varying the parameters of the sequencing process allows us to test analysis pipelines across different coverages, read lengths, ploidies, and more. NEAT addresses this need as an open-source Python package that allows users to customize their sequencing data [@Stephens:2016].

NEAT is designed to simulate short reads from sequencing platforms with machine-specific custom error models and can account for single-base substitutions, insertions, deletions, and larger structural variants. Unlike simulators that rely on fixed statistical profiles, NEAT learns empirical mutation and sequencing statistics from real data. NEAT models realistic, benchmarking-ready sequencing data, providing outputs in common bioinformatics file formats, including FASTQ, binary alignment map (BAM), and variant call format (VCF) files.

# Software design

NEAT has changed significantly since the release of version 2.0 in 2016. The codebase was updated to Python 3 to take advantage of the latest libraries and then was optimized for speed and accuracy alongside several new features. A summary of changes to NEAT is provided in **Table 1**.

### Table 1. Methodological, software robustness, and user experience updates in NEAT 4.4

+--------+-----------------------+---------------------------------------------+------------------------------------------------------+
| #      | Feature Name          | Prior Implementation (2.0)                  | Updated Implementation (4.4)                         |
+========+=======================+=============================================+======================================================+
| **1**  | GC Bias Computation   | Used custom model of GC bias                | Removed, pending further investigation               |
+--------+-----------------------+---------------------------------------------+------------------------------------------------------+
|        |                       |                                             |                                                      |
+--------+-----------------------+---------------------------------------------+------------------------------------------------------+
| **2**  | Read Generation       | Sliding-window approach                     | Coordinate-based read selection                      |
+--------+-----------------------+---------------------------------------------+------------------------------------------------------+
|        |                       |                                             |                                                      |
+--------+-----------------------+---------------------------------------------+------------------------------------------------------+
| **3**  | Variant Input         | Partial data loss from input variants       | Preserves all input data                             |
+--------+-----------------------+---------------------------------------------+------------------------------------------------------+
|        |                       |                                             |                                                      |
+--------+-----------------------+---------------------------------------------+------------------------------------------------------+
| **4**  | Variant Modeling      | Two variant types supported                 | Framework to expand variant types                    |
+--------+-----------------------+---------------------------------------------+------------------------------------------------------+
|        |                       |                                             |                                                      |
+--------+-----------------------+---------------------------------------------+------------------------------------------------------+
| **5**  | Automated Testing     | No formal testing framework                 | Unit tests and automated continuous integration      |
+--------+-----------------------+---------------------------------------------+------------------------------------------------------+
|        |                       |                                             |                                                      |
+--------+-----------------------+---------------------------------------------+------------------------------------------------------+
| **6**  | Configuration Files   | Command line interface only                 | Structured configuration files for reproducibility   |
+--------+-----------------------+---------------------------------------------+------------------------------------------------------+
|        |                       |                                             |                                                      |
+--------+-----------------------+---------------------------------------------+------------------------------------------------------+
| **7**  | Detailed Logging      | Minimal error logging                       | Extensive logs to recreate runs and describe errors  |
+--------+-----------------------+---------------------------------------------+------------------------------------------------------+
|        |                       |                                             |                                                      |
+--------+-----------------------+---------------------------------------------+------------------------------------------------------+
| **8**  | Friendly Installation | Clone repository and install                | pip installable                                      |
+--------+-----------------------+---------------------------------------------+------------------------------------------------------+
|        |                       |                                             |                                                      |
+--------+-----------------------+---------------------------------------------+------------------------------------------------------+
| **9**  | Parallelization       | Single-threaded                             | Multi-threaded                                       |
+--------+-----------------------+---------------------------------------------+------------------------------------------------------+
|        |                       |                                             |                                                      |
+--------+-----------------------+---------------------------------------------+------------------------------------------------------+
| **10** | Refactored Unit       | None                                        | Added for all major functions                        |
|        | Testing               |                                             |                                                      |
+--------+-----------------------+---------------------------------------------+------------------------------------------------------+

Below, we summarize methodological changes (**Table 1**) present in NEAT 4.4:


- [**1**] Guanine-cytosine (GC) bias refers to sequencing machine bias in GC-rich or GC-poor regions, causing uneven coverage results in real data [@Benjamini:2012; @Ross:2013]. Evidence suggests GC bias may arise from library preparation and amplification rather than from sequencing [@Benjamini:2012]. This is an area of ongoing research and development for NEAT.


- [**2**] NEAT 4.4's read generation algorithm eliminates artificial gaps in the output and now facilitates parallelization.


- [**3**] Variants can be incorporated into NEAT-generated reads via user-inputted VCF files. NEAT 2.0 reads only minimal data of the variant, but this has been expanded in NEAT 4.4.


- [**4**] NEAT 4.4 uses the same variants as NEAT 2.0, but the code has been updated to allow for additional variant types in future releases.


Documentation of software robustness, user experience, parallelization performance benchmarks, and more future improvements [**5**–**10**] is available on the project's online repository.

# State of the field

Even as long-read sequencing advances, short-read, bulk sequencing remains prominent due to its comparatively low cost and high throughput. Simulating short-read datasets can replicate sequencing pipelines used in a wide variety of research settings. Investigations of read simulators have analyzed use cases across whole genomes, exomes, and metagenomes [@Escalona:2016; @Zhao:2017] and whether empirical error-profile learning improves realism [@Alosaimi:2020; @Milhaven:2023; @Schmeing:2021]. While many short-read simulators have appeared in these studies (ART, CuReSim, DWGSIM, GemSIM, InSilicoSeq, Mason, NEAT, pIRS, ReSeq, SInC, and wgsim [@Huang:2012; @Caboche:2014; @Homer:2010; @McElroy:2012; @Gourle:2019; @Holtgrewe:2010; @Hu:2012; @Schmeing:2021; @Pattnaik:2014; @Li:2011]), only a subset were found to produce sequencing data with explicit ground truth suitable for benchmarking [@Alosaimi:2020], including NEAT [@Milhaven:2023].

Additionally, in their benchmark of twenty DNA read simulators that use reference genomes and produce FASTQ files, Alosaimi et al. (2020) found that NEAT's suite of features compares favorably to other tools. Although NEAT achieved the second-highest mapping sensitivity and precision on a human chromosome 22 test set, NEAT's runtimes were the second-longest [@Alosaimi:2020]. Milhaven and Pfeifer also noted NEAT's realism in sequencing—but also its slow simulation runtimes [@Milhaven:2023]. NEAT is noteworthy for its ability to combine mutation and sequencing models in a single framework and accept user-specified mutation models [@Stephens:2016]. The latest version of NEAT maintains these strengths and addresses some of these weaknesses with multi-threading and algorithmic updates, as outlined above.

# Research impact statement

NEAT is notable for producing ground-truth BAM and VCF outputs suitable for systematic benchmarking pipelines using custom mutation and sequencing models [@Alosaimi:2020; @Schmeing:2021; @Stephens:2016; @Milhaven:2023]. Researchers have been using NEAT since its initial release—from assisting benchmarks in the sequencing of the human Y chromosome [@Rhie:2023] to evaluating other bioinformatics tools [@Lefouili:2022; @Zhao:2020]. NEAT supports method development across a range of bioinformatics tasks, including optimization of high-throughput variant-calling workflows [@Ahmed:2019; @Kendig:2019], feasibility studies of exome sequencing [@RuizSchultz:2021], pan-genome mapping [@Jandrasits:2019], Bayesian approaches for resolving ambiguously mapped reads [@Shah:2021], and ultra-sensitive multi-sample variant calling [@Delhomme:2020].

The updates described here (NEAT 4.4) highlight NEAT's practicality for future work in the field. The source code for NEAT is freely available on GitHub [@Stephens:2016].

# Acknowledgements

We thank the original creators of NEAT: Zachary D. Stephens, Matthew E. Hudson, Liudmila S. Mainzer, Morgan Taschuk, Matthew R. Weber, and Ravishankar K. Iyer.

We also thank Varenya Jain, Meredith Pudlewski, Karen H. Xiong, and other contributors for their work on updating NEAT.

Portions of this project were funded by the National Center for Supercomputing Applications' Students Pushing Innovation (SPIN) program and the Illinois Computes project through the University of Illinois Urbana-Champaign.

# References
