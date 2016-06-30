from output_formatter import OutputFormatter

o = OutputFormatter()

def inspect(obj):
    """
    Just for debugging and exploration
    """
    o.unmute()
    o.print(obj)
    o.indent()
    o.print(type(obj))
    for k, v in obj.__dict__.items():
        o.print(k, v)
    o.print(dir(obj))
    o.outdent()


