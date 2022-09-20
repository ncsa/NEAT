"""Module with definitions of classes used by CLI."""

__all__ = ["BaseCommand", "Group", "Option"]

import abc
import argparse
from typing import Any


class BaseOption(abc.ABC):
    """
    A common interface for reusable options and groups of options.
    """

    @abc.abstractmethod
    def __init__(self):
        """
        Initialize the object.
        """

    @abc.abstractmethod
    def add_to_parser(self, parser: argparse.ArgumentParser):
        """
        Add an option to a parser

        :param parser: argparse.ArgumentParser object to add argument to
        """


class Option(BaseOption):
    """
    A reusable option

    :param args: Either a name or a list of options strings, e.g., 'foo' or '-f, --foo'
    :param kwargs: Any combination of named arugments accepted by argparse.add_argument().
    """

    def __init__(self, *args: str, **kwargs: Any):
        self.args: tuple[str, ...] = args
        self.kwargs: dict[str, Any] = kwargs

    def add_to_parser(self, parser: argparse.ArgumentParser):
        """
        Add the option to a parser

        :param parser: ArgumentParser object to add argument to
        """
        parser.add_argument(*self.args, **self.kwargs)

    def add_to_group(self, group: argparse._ArgumentGroup):
        """
        Add the opiton to an argument group.

        :param group: The group to add the option to
        """
        group.add_argument(*self.args, **self.kwargs)


class Group(BaseOption):
    """
    A reusable argument group.

    :param name: Name of the group.
    :param description: Description of the group.
    :param is_mutually_exclusive: If set, the arguments in the group will be mutually exclusive, i.e.,
        only one fo the arguments in the group can be present on the command
        line. Defaults to False.
    :param required: If set, at least one of the mutually exclusive arguments is required.
        Defaults to False.
    """
    def __init__(
        self,
        name: str | None = None,
        description: str | None = None,
        is_mutually_exclusive: bool = False,
        required: bool = False,
    ):
        self.name = name
        self.description = description
        self.is_mutually_exclusive = is_mutually_exclusive
        self.required = required
        self.options: list[Option] = []

    def add_argument(self, *args: Any, **kwargs: Any):
        """
        Add an argument to the group

        :param args: Either a name or a list of option strings
        :param kwargs: Any named argument accepted by argparse.add_argument()
        """
        if args and isinstance(args[0], Option):
            self.options.append(args[0])
        else:
            self.options.append(Option(*args, **kwargs))

    def add_to_parser(self, parser: argparse.ArgumentParser) -> None:
        """
        Add the group to an argument parser.

        :param parser: Argument parser to add the group to.
        """
        group: argparse._ArgumentGroup
        if self.is_mutually_exclusive:
            group = parser.add_mutually_exclusive_group(required=self.required)
        else:
            group = parser.add_argument_group(
                title=self.name, description=self.description
            )
        for option in self.options:
            option.add_to_group(group)


class BaseCommand(abc.ABC):
    """
    A CLI subcommand

    All subcommands should inherit from this base class

    :param parser: Main argument parser to which the subcommand options and arguments needs to be added.
    """
    name: str | None = None
    """
    Name of the subcommand.
    """

    description: str | None = None
    """
    The subcommand's help string. If not given, __doc__ will be used.
    """

    options: list[Option] | None = None
    """
    A list of predefined options.
    """
    def __init__(self, parser: argparse.ArgumentParser):
        if self.options:
            for opt in self.options:
                opt.add_to_parser(parser)
        self.add_arguments(parser)

    @abc.abstractmethod
    def add_arguments(self, parser: argparse.ArgumentParser):
        """
        Add arguments to the subcommand's argument parser.

        :param parser: The parser to add arguments to
        """

    @abc.abstractmethod
    def execute(self, arguments: argparse.Namespace):
        """
        Execute the command.

        :param arguments : The namespace with arguments and their values.
        """

    @classmethod
    def register_to(cls, subparsers: argparse._SubParsersAction, name: str | None):
        """
        Register subcommand's parser with the main command argument parser.

        Adds command specific subparser to the list of subparsers known by
        the main command argument parser.

        :param subparsers: argparse object representing subparsers.
        :param name: Name of the subcommand. If not given, defaults to None which means
            that the class attribute 'name' will be used instead.
        """
        cmd_name = name or cls.name
        help_text = cls.description or cls.__doc__
        parser = subparsers.add_parser(cmd_name, description=help_text, help=help_text)
        command = cls(parser)
        parser.set_defaults(cmd_handler=command.execute)
        parser.set_defaults(cmd_name=cmd_name)
