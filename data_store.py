class DataStore(object):
    """When inheriting from DataStore, make sure you are implementing
    'attributes' list in your __init__ function.
    """
    def __init__(self, data_raw, data, iterator=None):
        self.data = data
        self.data_raw = data_raw

        if iterator:
            self.iterator = iterator
        else:
            self.iterator = self.data.__iter__

    def __iter__(self):
        return self

    def next(self):

        entry = self.iterator.next()
        parsed = self.data_raw(entry)

        return parsed
