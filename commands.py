from collections import OrderedDict
import argparse
from argparse import ArgumentParser
from gettext import gettext as _

COMMANDS = OrderedDict()

argparse.PARSER =  'A...?'

class ArgumentParserPlus(ArgumentParser):
    """A set of hacks and backports to make ArgParse subparsers:
        - optional,
        - possible to group in more than one group
    See: http://bugs.python.org/issue9253 among others.
    """

    OPTIONAL_SUBPARSER = 'A...?'

    def add_subparsers(self, **kwargs):
        """This function code is a modification of python argparse code,
        python licence appiles."""
        if self._subparsers is None:
            self._subparsers = []

        # add the parser class to the arguments if it's not present
        kwargs.setdefault('parser_class', type(self))

        if 'title' in kwargs or 'description' in kwargs:
            title = _(kwargs.pop('title', 'subcommands'))
            description = _(kwargs.pop('description', None))
            subparsers_group = self.add_argument_group(title, description)
        else:
            subparsers_group = self._positionals

        self._subparsers.append(subparsers_group)

        # prog defaults to the usage message of this parser, skipping
        # optional arguments and with no "usage:" prefix
        if kwargs.get('prog') is None:
            formatter = self._get_formatter()
            positionals = self._get_positional_actions()
            groups = self._mutually_exclusive_groups
            formatter.add_usage(self.usage, positionals, groups, '')
            kwargs['prog'] = formatter.format_help().strip()

        # create the parsers action and add it to the positionals list
        parsers_class = self._pop_action_class(kwargs, 'parsers')
        action = parsers_class(option_strings=[], **kwargs)
        self._subparsers[-1]._add_action(action)

        action.nargs = self.OPTIONAL_SUBPARSER

        # return the created parsers action
        return action

    def _get_nargs_pattern(self, action):
        nargs = action.nargs
        if nargs and nargs == self.OPTIONAL_SUBPARSER:
            return '(-*A?[-AO]*)?'
        return super(ArgumentParserPlus, self)._get_nargs_pattern(action)


    def _get_values(self, action, arg_strings):
        nargs = action.nargs
        if nargs and nargs == self.OPTIONAL_SUBPARSER:
            if (not arg_strings and action.nargs == self.OPTIONAL_SUBPARSER):
                return [action._name_parser_map.keys()[0]]
        return super(ArgumentParserPlus, self)._get_values(action, arg_strings)


class SubparserGroupsRegistry(type):
    registry = []

    def __new__(cls, name, bases, attrs):
        new_class = super(SubparserGroupsRegistry, cls).__new__(cls, name, bases, attrs)
        new_class.registry = OrderedDict()
        SubparserGroupsRegistry.registry.append(new_class)
        return new_class


class Subparser(object):

    __metaclass__ = SubparserGroupsRegistry

    name = 'generic subcommands'

    def __init__(self, *args, **kwargs):

        self.name = args[0]

        self.args = args
        self.kwargs = kwargs
        self.commands = {}

        self.registry[self.name] = self

    def command(self, *args, **kwargs):
        if args in self.commands:
            raise ValueError('%s is already defined!' % args)

        def decorator(handler=None):
            self.commands[args] = (kwargs, handler)

            return handler

        return decorator

    def add_command(self, *args, **kwargs):
        if args in self.commands:
            raise ValueError('%s is already defined!' % args)
        self.commands[args] = (kwargs, None)


class SourceSubparser(Subparser):
    name = 'source'


class AnalysisSubparser(Subparser):
    name = 'analysis'



def append_subparsers(parser):
    """Add subparsers to given parser"""

    for subparser_group in SubparserGroupsRegistry.registry:

        name = subparser_group.name
        subparsers = subparser_group.registry

        if subparsers:
            sps = parser.add_subparsers(
                title=name,
                description=name + ' arguments: ' + ', '.join(subparsers.keys())
            )

            for args, subparser in subparsers.items():
                sp = sps.add_parser(*subparser.args, **subparser.kwargs)
                append_commands(sp, subparser.commands)


def command(*args, **kwargs):
    """Decorator registering commands to be invoked from argparse"""

    if args in COMMANDS:
        raise ValueError('%s is already defined!' % args)

    def decorator(handler=None):
        COMMANDS[args] = (kwargs, handler)

        return handler

    return decorator


def execute_commands(args):
    """Execute commands registered with decorator 'command',
    using provided 'args' arguments - result of argparse parsing.
    """
    for cmd_args, command_tuple in COMMANDS.items():
        _, handler = command_tuple
        arg = getattr(args, cmd_args[0].lstrip('-'))
        if handler:
            handler(arg, args)


def append_commands(parser, commands=COMMANDS):
    """Add arguments to given parser"""

    for args, command_tuple in commands.items():
        kwargs, func = command_tuple
        parser.add_argument(*args, **kwargs)
