COMMANDS = {}


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


def append_commands(parser):
    """Add arguments to given parser"""

    for args, command_tuple in COMMANDS.items():
        kwargs, func = command_tuple
        parser.add_argument(*args, **kwargs)
