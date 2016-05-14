

class DataStore(object):
    """When inheriting from DataStore, make sure you are implementing
    'attributes' list in your __init__ function.
    """
    def __init__(self, attributes, data, iterator=None, single_rows=True):
        self.attributes = attributes
        self.data = data
        self.prev = ''
        self.single_rows = single_rows
        if iterator:
            self.iterator = iterator
        else:
            self.iterator = self.data.__iter__

    def pos(self, attribute):
        try:
            return self.attributes.index(attribute)
        except ValueError:
            raise KeyError

    def parse_single(self, entry):
        return DataRow(self, entry)

    def parse_multi(self, entry):
        return entry

    def parse_entry(self, entry):
        if self.single_rows:
            return self.parse_single(entry)
        else:
            return self.parse_multi(entry)

    def __iter__(self):
        return self

    def is_ready_to_flush(self, data):
        return True

    def next(self):

        if not self.prev:
            self.prev = self.iterator.next()
            prev = ''

        rows = []

        is_end = False

        try:
            while not self.is_ready_to_flush(rows):
                prev = self.iterator.next()
                rows += [prev]

        except StopIteration:
            is_end = True

        entry = self.prev + '\n' + '\n'.join(rows[:-1])

        if is_end:
            entry += prev
            prev = None

        parsed = self.parse_entry(entry)
        self.prev = prev
        return parsed


class DataRow(object):
    """A class for storing csv/tsv row which gives a good trade off between
    easy keyword-based access and memmory-efficient store. The keywords are
    stored in parent class which should be inherited from DataStore class or
    which implements 'pos' function and 'attributes' list, with the interface
    the same as in DataStore.
    """

    def __init__(self, parent, line):
        if not isinstance(line, list) and not isinstance(line, tuple):
            line = line.decode('utf-8').split('\t')
        self._data = line
        self._parent = parent

    def __getattr__(self, k):
        if k.startswith('__') and k.endswith('__'):
            return super(DataRow, self).__getattr__(k)
        try:
            return self._data[self._parent.pos(k)]
        except KeyError:
            raise AttributeError

    def __setattr__(self, k, v):
        if k.startswith('_'):
            super(DataRow, self).__setattr__(k, v)
        else:
            try:
                self._data[self.get_index(k)] = v
            except KeyError:
                self._parent.attributes.append(k)
                self.__setattr__(k, v)
            except IndexError:
                while len(self._data) <= len(self._parent.attributes):
                    self._data.append(None)
                self.__setattr__(k, v)

    def __repr__(self):
        buff = 'DataRow:\n\t'
        data = [self.get_key(i) + ': ' + str(v) for i, v in enumerate(self._data)]
        buff += '\n\t'.join(data)

        return buff

    def get_key(self, index):
        return self._parent.attributes[index]

    def get_index(self, key):
        return self._parent.pos(key)

    def get(self, k):
        return self.__getattr__(k)

    def set(self, k, v):
        return self.__setattr__(k, v)

