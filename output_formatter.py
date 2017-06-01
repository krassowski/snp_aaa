level = 0


def formatter(obj, skip_tab=False, name=None):
    global level
    level += 1
    text = ''

    if name:
        text += name + ':' + '\n'

    if not skip_tab:
        text += '\t' * level

    if type(obj) is dict:
        text += '\n'
        for key in sorted(obj.keys()):
            value = obj[key]
            text += '\t' * level + repr(key) + ':' + formatter(value, True) + '\n'
    elif type(obj) in (list, tuple, set):
        text += '\n'
        text += '\t' * level + obj.__class__.__name__ + ': '
        text += '\n'
        for value in tuple(obj):
            text += '\t' * level + formatter(value, True) + '\n'
    else:
        text += '\t' + repr(obj)

    level -= 1
    return text


