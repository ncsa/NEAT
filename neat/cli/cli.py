"""Implements command line interface used by the package."""

__all__ = ['Cli', 'main', 'run']

import argparse
import importlib
import logging
import traceback
import time
import pkgutil
import sys
import os

from typing import Final
from datetime import datetime

from ..common import setup_logging
from .commands import BaseCommand

log = logging.getLogger("neat")

COMMANDS_MODULE_PATH: Final = importlib.import_module("neat.cli.commands").__path__


class Cli:
    """
    NEAT command line interface.
    """

    parser: argparse.ArgumentParser
    """
    Main command line parser.
    """

    subparsers: argparse._SubParsersAction
    """
    Object storing subcommand parsers
    """

    def __init__(self):
        self.parser = argparse.ArgumentParser(
            prog="neat", description="Run NEAT components"
        )
        self.parser.add_argument(
            "--no-log",
            default=False,
            action='store_true',
            help="Set to turn off log file creation."
        )
        self.parser.add_argument(
            "--log-dir",
            type=str,
            default=os.getcwd(),
            help="directory where to write the log file (default is working directory)"
        )
        self.parser.add_argument(
            "--log-name",
            type=str,
            default=f"{time.time()}_NEAT.log",
            help="Name of the log file to produce, if producing."
        )
        self.parser.add_argument(
            "--log-level",
            choices=["DEBUG", "INFO", "WARN", "WARNING", "ERROR"],
            default="INFO",
            help="Severity level of the log messages to display"
        )
        self.parser.add_argument(
            "--log-detail",
            choices=["LOW", "MEDIUM", "HIGH"],
            default="MEDIUM",
            help="Level of detail to include in the log message"
        )
        self.parser.add_argument(
            "--silent-mode",
            default=False,
            action="store_true",
            help="If entered, this will suppress messages to stdout"
        )

        self.subparsers = self.parser.add_subparsers()

        # Discover and register existing commands.
        for _, name, _ in pkgutil.iter_modules(COMMANDS_MODULE_PATH):
            module = importlib.import_module(f"neat.cli.commands.{name}")
            try:
                cls = module.Command
            except AttributeError:
                continue
            self.register_command(cls, cls.name or name)

    def register_command(self, command: BaseCommand, name: str | None):
        """
        Register a subcommand

        :param command: The command class to register
        :param name: The name of the subcommand. If not given, `command.name` will be used
        """
        assert self.subparsers
        command.register_to(self.subparsers, name)


def main(parser: argparse.ArgumentParser, arguments: list[str]) -> int:
    """
    The command line entry point.

    :param parser: The argument parser
    :param arguments: The list of command line arguments
    :return: Exit code, 0 if command executed successfully, a positive integer otherwise
    """
    try:
        args = parser.parse_args(arguments)
    except SystemExit:
        return 2

    setup_logging(
        omit_log=args.no_log,
        severity=args.log_level,
        verbosity=args.log_detail,
        directory=args.log_dir,
        filename=args.log_name,
        silent_mode=args.silent_mode
    )

    try:
        cmd, name = args.cmd_handler, args.cmd_name
    except AttributeError:
        parser.print_help()
        return 1
    else:
        start = time.time()
        try:
            cmd(args)
        except Exception as exc:
            log.exception(f"{name} failed, see the traceback below")
            err_msg = f"ERROR: {name} failed, showing the last error"
            print(err_msg)
            traceback.print_exception(exc, chain=False)
            return 1
        else:
            end = time.time()
            log.info(
                f"command finished successfully; execution took {(end - start)/60:.2f} m"
            )
        return 0


def run():
    """
    Console script entry point
    """

    cli = Cli()
    rc = main(cli.parser, sys.argv[1:])
    sys.exit(rc)
