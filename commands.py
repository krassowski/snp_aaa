from collections import OrderedDict
import argparse
from argparse import ArgumentParser
from gettext import gettext as _

COMMANDS = OrderedDict()

argparse.PARSER = 'A...?'


class ArgumentParserPlus(ArgumentParser):
    """A set of hacks and backports to make ArgParse subparsers:
        - optional,
        - possible to group in more than one group
    See: http://bugs.python.org/issue9253 among others.
    """

    OPTIONAL_SUBPARSER = 'A...?'

    def __init__(self, *args, **kwargs):
        super(ArgumentParserPlus, self).__init__(*args, **kwargs)
        self.queued_subparser_value = None

    def add_subparsers(self, **kwargs):
        """This function code is a modification of python argparse code,
        python licence applies."""
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
        if nargs == self.OPTIONAL_SUBPARSER:
            return '(-*A?[-AO]*)?'
        return super(ArgumentParserPlus, self)._get_nargs_pattern(action)

    def _get_values(self, action, arg_strings):
        nargs = action.nargs

        if nargs == self.OPTIONAL_SUBPARSER:

            if not arg_strings and self.queued_subparser_value:
                arg_strings = self.queued_subparser_value
                self.queued_subparser_value = None

            if arg_strings:

                if action.option_strings:
                    value = action.const
                else:
                    value = action.default

                if not value:
                    this_group_subparsers = action._name_parser_map.keys()

                    value = [self._get_value(action, v) for v in arg_strings]

                    # the parser choice does not match this subparsers group:
                    if value[0] not in this_group_subparsers:
                        # save the argument strings until the matching subparsers
                        # group will be called. As we consumed the arg_strings
                        # the _get_values calls for next subparsers groups will have
                        # empty arg_strings -> we can then substitute our queued value
                        # and again check if we the values match to give subparser
                        self.queued_subparser_value = arg_strings

                        # mock selection of a valid subparser, but without arguments
                        return [action._name_parser_map.keys()[0]]

                    return value

            else:
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
        self.commands = OrderedDict()

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

    def add_to_argparse(self, subparser_group):
        """subparser_group is of a "parser" type"""
        _subparser = subparser_group.add_parser(*self.args, **self.kwargs)
        append_commands(_subparser, self.commands)
        return _subparser

    def execute_commands(self, args):
        return execute_commands(args, self.commands)


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
            _subparser_group = parser.add_subparsers(
                title=name,
                description=name + ' arguments: ' + ', '.join(subparsers.keys())
            )

            for subparser in subparsers.values():
                subparser.add_to_argparse(_subparser_group)


def execute_subparser_commands(args):
    """Execute commands registered with decorator 'command',
    using provided 'args' arguments - result of argparse parsing.
    """
    for subparser_group in SubparserGroupsRegistry.registry:

        subparsers = subparser_group.registry

        for subparser in subparsers.values():
            subparser.execute_commands(args)


def command(*args, **kwargs):
    """Decorator registering commands to be invoked from argparse"""

    if args in COMMANDS:
        raise ValueError('%s is already defined!' % args)

    def decorator(handler=None):
        COMMANDS[args] = (kwargs, handler)

        return handler

    return decorator


def execute_commands(args, commands=COMMANDS):
    """Execute commands registered with decorator 'command',
    using provided 'args' arguments - result of argparse parsing.
    """
    for cmd_args, command_tuple in commands.items():
        _, handler = command_tuple
        try:
            arg = getattr(args, cmd_args[0].lstrip('-'))
            if handler:
                handler(arg, args)
        except AttributeError:
            # silently skip commands from parsers which are absent
            pass


def append_commands(parser, commands=COMMANDS):
    """Add arguments to given parser"""

    for args, command_tuple in commands.items():
        kwargs, func = command_tuple
        parser.add_argument(*args, **kwargs)
