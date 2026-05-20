import argparse

from ...model_quality_score import model_qual_score_runner
from ...quality_score_modeling.presets import QUALITY_PRESETS
from .base import BaseCommand
from .options import output_group


class Command(BaseCommand):
    """
    Generate a quality score model (traditional or Markov) from FASTQ.
    """

    name = "model-qual-score"
    description = "Generate quality score model from FASTQ (optional Markov chain)."

    def add_arguments(self, parser: argparse.ArgumentParser):

        parser.add_argument(
            "-i",
            dest="input_files",
            metavar="FILE",
            nargs="+",
            required=True,
            help="Input FASTQ file(s) (gzipped or plain).",
        )

        parser.add_argument(
            "-q",
            dest="quality_offset",
            type=int,
            default=33,
            help="Quality score offset [33].",
        )

        parser.add_argument(
            "-Q",
            dest="quality_scores",
            type=int,
            nargs="+",
            default=[42],
            help="Maximum quality score (single int) or explicit list of bin values "
                 "(e.g. -Q 2 12 23 37). A list enables Markov binning: all observed "
                 "scores are down-binned to the nearest value and simulation output is "
                 "constrained to those values. Overridden by --quality-preset. [42]",
        )

        parser.add_argument(
            "--quality-preset",
            dest="quality_preset",
            choices=list(QUALITY_PRESETS),
            default=None,
            metavar="PRESET",
            help="Named bin preset for common Illumina instruments. Implies --markov. "
                 f"Choices: {', '.join(QUALITY_PRESETS)}. Overrides -Q when set.",
        )

        parser.add_argument(
            "-m",
            dest="max_num",
            type=int,
            default=-1,
            help="Max number of reads to process [-1 = all].",
        )

        parser.add_argument(
            "--markov",
            dest="use_markov",
            action="store_true",
            default=False,
            help="Use Markov quality model instead of the traditional model.",
        )

        parser.add_argument(
            "--overwrite",
            dest="overwrite",
            action="store_true",
            default=False,
            help="Overwrite existing output file if present.",
        )

        output_group.add_to_parser(parser)

    def execute(self, arguments: argparse.Namespace):

        if arguments.quality_preset:
            qual_scores: int | list[int] = QUALITY_PRESETS[arguments.quality_preset]
            use_markov = True
        elif len(arguments.quality_scores) == 1:
            qual_scores = arguments.quality_scores[0]
            use_markov = arguments.use_markov
        else:
            qual_scores = arguments.quality_scores
            use_markov = arguments.use_markov

        model_qual_score_runner(
            files=arguments.input_files,
            offset=arguments.quality_offset,
            qual_scores=qual_scores,
            max_reads=arguments.max_num,
            overwrite=arguments.overwrite,
            output_dir=arguments.output_dir,
            output_prefix=arguments.prefix,
            use_markov=use_markov,
        )
