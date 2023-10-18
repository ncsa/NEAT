import logging

_LOG = logging.getLogger(__name__)


def generate_bam(reference_chrom, options, temporary_vcf, temporary_dir, chrom):
    # Step 4 (optional): Create a golden bam with the reads aligned to the original reference
    print(reference_chrom, options, temporary_vcf, temporary_dir, chrom)

