import argparse

from ...model_quality_score import model_qual_score_runner
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
            help="Max quality or explicit list of quality scores [42].",
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

        if len(arguments.quality_scores) == 1:
            qual_scores: int | list[int] = arguments.quality_scores[0]

        else:
            qual_scores = arguments.quality_scores

        model_qual_score_runner(
            files=arguments.input_files,
            offset=arguments.quality_offset,
            qual_scores=qual_scores,
            max_reads=arguments.max_num,
            overwrite=arguments.overwrite,
            output_dir=arguments.output_dir,
            output_prefix=arguments.prefix,
            use_markov=arguments.use_markov,
        )
