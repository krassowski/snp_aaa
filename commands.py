from collections import OrderedDict

COMMANDS = OrderedDict()


class SubparserGroupsRegistry(type):
    registry = set()

    def __new__(cls, name, bases, attrs):
        new_class = super().__new__(cls, name, bases, attrs)
        new_class.registry = OrderedDict()
        cls.registry.add(new_class)
        cls.registry -= set(bases)
        return new_class


class Subparser(object, metaclass=SubparserGroupsRegistry):

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
    name = 'source-options'


class AnalysisSubparser(Subparser):
    name = 'analysis-options'


def append_subparsers(parser):
    """Add subparsers to given parser"""

    subcommands = parser.add_subparsers(help='sub-commands')

    for group in SubparserGroupsRegistry.registry:
        group.parser = subcommands.add_parser(help=group.name, name=group.name)


    for group in SubparserGroupsRegistry.registry:

        name = group.name
        subparsers = group.registry

        if subparsers:
            _subparser_group = group.parser.add_subparsers(
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
