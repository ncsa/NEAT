# Summary

While the field of genomics has advanced significantly with the advent of high-throughput sequencing technologies, challenges related to the proprietary access, availability, complexity, and variability of this data can introduce difficulty in the development and validation of computational tools. Simulated datasets provide ground-truth estimates and reproducible scalability that are important in testing algorithms and benchmarking software. Ideally, these datasets mimic the properties of real sequencing datasets—from introducing specific patterns of sequencing errors to localized regions of mutations that often characterize cancer genomes.

# Statement of need

The NExt-generation sequencing Analysis Toolkit (NEAT) is an open-source Python package that creates simulated next-generation sequencing datasets. NEAT’s simulations account for a wide range of sequencing parameters (e.g., DNA read fragment length, sequencing error rates, mutation frequencies, etc.) and allow users to customize their sequencing data.<sup>1</sup> Since the original release of NEAT in 2016, most scripts have been greatly modified, and NEAT is currently on version 4.2. The code has undergone significant ongoing changes since 2020. Upgrading to Python 3 has enabled NEAT to achieve a flexible and intuitive user interface with minimal dependencies. The toolkit is optimized for both speed and accuracy, and new features have been implemented, such as improved ploidy simulation, mutation modeling, and the ability to model mutational profiles directly from data. A summary of algorithmic changes is provided in Table 1.

NEAT can integrate seamlessly with existing bioinformatics workflows, providing outputs in several common file formats. The toolkit’s ability to simulate gold-standard synthetic datasets with ground truth annotations is useful for testing bioinformatics pipelines. Uses of NEAT continue to be prominently featured—from scientists who have comprehensively sequenced the human Y chromosome<sup>2</sup> to researchers who use NEAT to evaluate and validate the performance of other high-profile bioinformatics tools.<sup>3, 4</sup> Earlier versions of NEAT have also demonstrated utility when benchmarked in comparison to similar tools.<sup>5</sup> The source code for both original and updated versions of NEAT is freely available on GitHub.<sup>1</sup>

# Tables

# Acknowledgements

# References
1.	Stephens ZD, Hudson ME, Mainzer LS, Taschuk M, Weber MR, Iyer RK. Simulating Next-Generation Sequencing Datasets from Empirical Mutation and Sequencing Models. PLOS ONE. 2016;11(11). doi:10.1371/journal.pone.0167047
2.	Rhie A, Nurk S, Cechova M, et al. The complete sequence of a human Y chromosome. Nature. 2023;621(7978):344-354. doi:10.1038/s41586-023-06457-y
3.	Lefouili M, Nam K. The evaluation of Bcftools mpileup and GATK HaplotypeCaller for variant calling in non-human species. Sci Rep. 2022;12(1). doi:10.1038/s41598-022-15563-2
4.	Zhao S, Agafonov O, Azab A, Stokowy T, Hovig E. Accuracy and efficiency of germline variant calling pipelines for human genome data. Scientific Reports. 2020;10(1). doi:10.1038/s41598-020-77218-4
5.	Alosaimi S, Bandiang A, van Biljon N, et al. A broad survey of DNA sequence data simulation tools. Brief Funct Genomics. 2020;19(1):49-59. doi:10.1093/bfgp/elz033
