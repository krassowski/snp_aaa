from __future__ import print_function
import cPickle as pickle
import os
import sys


class CachingError(Exception):
    pass


class cacheable(object):
    """Makes caching a reasonable pleasure. Use as a decorator."""

    def __init__(self, func):
        self.func = func
        self.display_name = func.__name__.replace('_', ' ').title()
        self.cache_name = '.cache:' + func.__name__

    def save(self, *args, **kwargs):
        result = self.create(*args, **kwargs)
        print('"' + self.display_name + '" saved to cache')

        with open(self.cache_name, 'wb') as f:
            pickle.dump(result, f, protocol=pickle.HIGHEST_PROTOCOL)

        return result

    def create(self, *args, **kwargs):
        print('Generating "' + self.display_name + '"...')
        result = self.func(*args, **kwargs)
        return result

    def load_or_create(self, *args, **kwargs):
        """Load or create results of cacheable function.

        If unable to find a file with saved cache, the decorated
        function will be executed and results will be saved.

        If args or kwargs are specified, those will be used
        during decorated function execution.
        """
        if os.path.exists(self.cache_name):
            return self.load()
        else:
            print(
                'Cache data for %s do not exists - running %s' %
                (self.display_name, self.func.__name__)
            )
            try:
                return self.save(*args, **kwargs)
            except TypeError as e:
                print(
                    'Whops! %s probably needs to be loaded manually as there are few arguments required' %
                    self.display_name
                )
                raise e

    def load(self):
        print('Loading %s' % self.display_name)
        with open(self.cache_name, 'rb') as f:
            try:
                result = pickle.load(f)
            except AttributeError as e:
                raise CachingError(
                    'Unable to load cache for %s. '
                    'Please, make sure that the cache files correspond'
                    ' to this version of application.\n'
                    'Original exception: %s.' %
                    (self.display_name, e.message)
                )
            return result
