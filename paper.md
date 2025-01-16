# Summary

While the field of genomics has advanced with the advent of high-throughput sequencing, the proprietary access, availability, complexity, and variability of this data all pose significant challenges to the development and validation of computational tools in genomics. Simulated datasets provide the ground-truth estimates and reproducible scalability that are vital to test algorithms and benchmark software. Ideally, these datasets could be created to mimic the properties of real sequencing datasets—from widespread sequencing errors to local mutational contexts that often characterize diseases like cancer.

# Statement of need

NEAT (NExt-generation sequencing Analysis Toolkit) is an open-source Python package that creates simulated next-generation sequencing datasets. NEAT’s simulations account for a wide range of sequencing parameters (e.g., read length, error rates, mutation frequencies, etc.) and allow users to customize their sequencing data.<sup>1</sup> Since the original release of NEAT in 2016, most scripts have been greatly modified, and NEAT is currently on version 4.2. Upgrading to Python 3 has enabled NEAT to achieve a flexible and intuitive user interface with minimal dependencies. The toolkit is optimized for both speed and accuracy, and new features have been implemented, such as ploidy simulation, mutation modeling, and the ability to model mutational profiles directly from data.

NEAT can integrate seamlessly with existing bioinformatics workflows, providing outputs in several common file formats. The toolkit’s ability to simulate gold-standard synthetic datasets with ground truth annotations is useful for testing bioinformatics pipelines. Uses of NEAT continue to be prominently featured—from scientists who have comprehensively sequenced the human Y chromosome<sup>2</sup> to researchers who use NEAT to evaluate and validate the performance of other high-profile bioinformatics tools.<sup>3</sup> Earlier versions of NEAT have also demonstrated utility when benchmarked in comparison to similar tools.<sup>4</sup> The source code for both original and updated versions of NEAT is freely available on GitHub.<sup>1</sup>

References:
1.	Stephens ZD, Hudson ME, Mainzer LS, Taschuk M, Weber MR, Iyer RK. Simulating Next-Generation Sequencing Datasets from Empirical Mutation and Sequencing Models. PLOS ONE. 2016;11(11). doi:10.1371/journal.pone.0167047
2.	Rhie A, Nurk S, Cechova M, et al. The complete sequence of a human Y chromosome. Nature. 2023;621(7978):344-354. doi:10.1038/s41586-023-06457-y
3.	Lefouili M, Nam K. The evaluation of Bcftools mpileup and GATK HaplotypeCaller for variant calling in non-human species. Sci Rep. 2022;12(1). doi:10.1038/s41598-022-15563-2
4.	Alosaimi S, Bandiang A, van Biljon N, et al. A broad survey of DNA sequence data simulation tools. Brief Funct Genomics. 2020;19(1):49-59. doi:10.1093/bfgp/elz033
