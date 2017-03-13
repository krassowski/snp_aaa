import os
import bsddb


class SetWithCallback(set):
    """A set implementation that triggers callbacks on `add` or `update`.
    It has an important use in BerkleyHashSet database implementation:
    it allows a user to modify sets like native Python's structures while all
    the changes are forwarded to the database, without additional user's action.
    """
    _modifying_methods = {'update', 'add'}

    def __init__(self, items, callback):
        super(SetWithCallback, self).__init__(items)
        self.callback = callback
        for method_name in self._modifying_methods:
            method = getattr(self, method_name)
            setattr(self, method_name, self._wrap_method(method))

    def _wrap_method(self, method):
        def new_method_with_callback(*args, **kwargs):
            result = method(*args, **kwargs)
            self.callback(self)
            return result
        return new_method_with_callback


class ListWithCallback(list):
    """A set implementation that triggers callbacks on `add` or `update`.
    It has an important use in BerkleyHashSet database implementation:
    it allows a user to modify sets like native Python's structures while all
    the changes are forwarded to the database, without additional user's action.
    """
    _modifying_methods = {'append', 'extend'}

    def __init__(self, items, callback):
        super(ListWithCallback, self).__init__(items)
        self.callback = callback
        for method_name in self._modifying_methods:
            method = getattr(self, method_name)
            setattr(self, method_name, self._wrap_method(method))

    def _wrap_method(self, method):
        def new_method_with_callback(*args, **kwargs):
            result = method(*args, **kwargs)
            self.callback(self)
            return result
        return new_method_with_callback


class BerkleyHash(object):

    def __init__(self, name):
        self.name = name
        self.path = self.create_path()
        self.open()

    def create_path(self):
        """Returns path to a file containing the database.
        The file is not guranted to exist, although the 'databases' directory
        will be created (if it does not exist).
        """
        base_dir = os.path.abspath(os.path.dirname(__file__))
        databases_dir = os.path.join(base_dir, 'databases')
        if not os.path.exists(databases_dir):
            os.makedirs(databases_dir)
        return os.path.join(databases_dir, self.name)

    def open(self, mode='c'):
        """Open hash database in a given mode.
        By default it opens a database in read-write mode and in case
        if a database of given name does not exists it creates one.
        """
        self.db = bsddb.hashopen(self.path, mode)

    def __iter__(self):
        return iter(self.db)

    def __len__(self):
        return len(self.db)

    def reset(self):
        """Reset database completely by its removal and recreation."""
        os.remove(self.path)
        self.open()


class BerkleyHashList(BerkleyHash):

    def __getitem__(self, key):
        """key: has to be str"""
        key = bytes(key)
        try:
            items = list(
                filter(
                    bool,
                    self.db.get(key).decode('utf-8').split('|')
                )
            )
        except (KeyError, AttributeError):
            items = []

        return ListWithCallback(
            items,
            lambda new_set: self.__setitem__(key, new_set)
        )

    def __setitem__(self, key, items):
        """key: might be a str or bytes"""
        assert '|' not in items
        if not isinstance(key, bytes):
            key = bytes(key)
        self.db[key] = bytes('|'.join(items))


class BerkleyHashSet(BerkleyHash):
    """A hash-indexed database where values are equivalent to Python's sets.
    It uses Berkley database for storage and accesses it through bsddb3 module.
    """

    def __getitem__(self, key):
        """key: has to be str"""
        key = bytes(key)
        try:
            items = list(
                filter(
                    bool,
                    self.db.get(key).decode('utf-8').split('|')
                )
            )
        except (KeyError, AttributeError):
            items = []

        return SetWithCallback(
            items,
            lambda new_set: self.__setitem__(key, new_set)
        )

    def __setitem__(self, key, items):
        """key: might be a str or bytes"""
        assert '|' not in items
        if not isinstance(key, bytes):
            key = bytes(key)
        self.db[key] = bytes('|'.join(items))

    def items(self):
        for key, value in self.db.iteritems():
            yield (
                key,
                SetWithCallback(
                    value.decode('utf-8').split('|'),
                    lambda new_set: self.__setitem__(key, new_set)
                )
            )
