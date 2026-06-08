"""
Command line interface for NEAT's compare-vcfs subcommand (issue #297).
"""
import argparse

from ...compare_vcfs import compare_vcfs_runner
from .base import BaseCommand


class Command(BaseCommand):
    """
    Compare a downstream variant-caller VCF against a NEAT-simulated truth VCF
    and attribute false negatives to the simulator's configuration.
    """
    name = "compare-vcfs"
    description = (
        "Compare a NEAT-simulated truth VCF (golden) against a downstream variant "
        "caller VCF (called), with NEAT-aware false-negative attribution."
    )

    def add_arguments(self, parser: argparse.ArgumentParser):
        parser.add_argument(
            "golden_vcf",
            type=str, metavar="golden.vcf",
            help="NEAT-simulated truth VCF (typically <prefix>_golden.vcf.gz)."
        )
        parser.add_argument(
            "called_vcf",
            type=str, metavar="called.vcf",
            help="Downstream variant caller's VCF produced from the simulated reads."
        )
        parser.add_argument(
            "--neat-run-dir",
            dest="neat_run_dir",
            type=str, required=True, metavar="DIR",
            help="Directory containing the NEAT simulator output, including simulation_summary.json."
        )
        parser.add_argument(
            "--output-dir",
            dest="output_dir",
            type=str, required=True, metavar="DIR",
            help="Where to write comparison_summary.{json,txt} and FN_with_reasons.vcf. Created if absent."
        )
        parser.add_argument(
            "--reference",
            type=str, default=None, metavar="ref.fa",
            help="Optional reference FASTA forwarded to hap.py."
        )
        parser.add_argument(
            "--target-bed",
            dest="target_bed",
            type=str, default=None, metavar="target.bed",
            help="Optional BED restricting comparison to these regions."
        )
        parser.add_argument(
            "--happy-bin",
            dest="happy_bin",
            type=str, default=None, metavar="PATH",
            help="Explicit path to the hap.py binary. Defaults to looking on $PATH."
        )
        parser.add_argument(
            "--plot",
            dest="plot",
            action="store_true",
            help="Also write fn_attribution.png — a bar chart of FN reason counts."
        )
        parser.add_argument(
            "--chrom-aliases",
            dest="chrom_aliases",
            type=str, default=None, metavar="TSV",
            help="Two-column TSV mapping BED chrom names to reference-canonical names, "
                 "applied to mutation_bed and target_bed at load time. Used when the BED "
                 "uses '1' but the reference uses 'chr1', or similar prefix/mt mismatches."
        )

    def execute(self, arguments: argparse.Namespace):
        compare_vcfs_runner(
            golden_vcf=arguments.golden_vcf,
            called_vcf=arguments.called_vcf,
            neat_run_dir=arguments.neat_run_dir,
            output_dir=arguments.output_dir,
            reference=arguments.reference,
            target_bed=arguments.target_bed,
            happy_bin=arguments.happy_bin,
            plot=arguments.plot,
            chrom_aliases=arguments.chrom_aliases,
        )
