import os
from abc import abstractmethod, abstractproperty, ABC

import bsddb3 as bsddb
from pathlib import Path


class CollectionWithCallback(ABC):
    """A collection implementation that triggers callbacks on actions
    specified by self._modifying_methods property.

    It has an important use in BerkleyHashSet database implementation:
    it allows a user to modify sets like native Python's structures while all
    the changes are forwarded to the database, without additional user's action.
    """

    @abstractproperty
    def _modifying_methods(self):
        pass

    def __init__(self, items, callback):
        super().__init__(items)
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


class SetWithCallback(CollectionWithCallback, set):
    _modifying_methods = {'update', 'add'}


class ListWithCallback(CollectionWithCallback, list):
    _modifying_methods = {'append', 'extend'}


class BerkleyHash:

    def __init__(self, name):
        self.name = name
        self.path = self.create_path()
        self.open()

    def create_path(self):
        """Returns path to a file containing the database.
        The file is not guaranteed to exist, although the 'databases' directory
        will be created (if it does not exist).
        """
        base_dir = Path.cwd()
        databases_dir = base_dir / 'databases'

        databases_dir.mkdir(exist_ok=True)

        return (databases_dir / self.name).as_posix()

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


class BerkleyHashCollection(ABC, BerkleyHash):

    @abstractproperty
    def python_type(self):
        pass

    @abstractmethod
    def callback_class(self, items, callback):
        pass

    def __setitem__(self, key, items):
        """
        Args:
            key: might be a str or bytes
            items: any iterable, will be casted to desired type
        """

        if not isinstance(key, bytes):
            key = bytes(key, 'utf-8')
        self.db[key] = bytes('|'.join(self.python_type(items)), 'utf-8')

    def __getitem__(self, key):
        """key: has to be str"""
        key = bytes(key, 'utf-8')
        try:
            items = list(
                filter(
                    bool,
                    self.db.get(key).decode('utf-8').split('|')
                )
            )
        except (KeyError, AttributeError):
            items = []

        return self.callback_class(
            items,
            lambda new_set: self.__setitem__(key, new_set)
        )

    def items(self):
        for key, value in self.db.items():
            yield (key, self[key])


class BerkleyHashList(BerkleyHashCollection):

    python_type = list
    callback_class = ListWithCallback


class BerkleyHashSet(BerkleyHashCollection):
    """A hash-indexed database where values are equivalent to Python's sets.
    It uses Berkley database for storage and accesses it through bsddb3 module.
    """

    python_type = set
    callback_class = SetWithCallback
