from __future__ import print_function
import pickle
import os
from pathlib import Path


class CachingError(Exception):
    pass


class cacheable(object):
    """Makes caching a reasonable pleasure. Use as a decorator."""

    cache_dir = Path('.cache')

    def __init__(self, func):
        self.func = func
        self.display_name = func.__name__.replace('_', ' ').title()
        self.cache_name = self.cache_dir / func.__name__
        self.cache_dir.mkdir(exist_ok=True)

    def file_name(self, *args, **kwargs):
        return self.cache_name

    def save(self, *args, **kwargs):
        result = self.create(*args, **kwargs)
        print('"' + self.display_name + '" saved to cache')
        file_name = self.file_name(*args, **kwargs)

        with file_name.open('wb') as f:
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
        file_name = self.file_name(*args, **kwargs)
        if file_name.exists():
            return self.load(*args, **kwargs)
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

    def load(self, *args, **kwargs):
        print('Loading %s' % self.display_name)
        file_name = self.file_name(*args, **kwargs)
        with file_name.open('rb') as f:
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


class args_aware_cacheable(cacheable):

    @staticmethod
    def safe_name(name):
        import string
        safe_chars = "_. %s%s" % (string.ascii_letters, string.digits)
        return ''.join([(c if c in safe_chars else '_')for c in str(name)])

    def file_name(self, *args, **kwargs):
        args_repr = ', '.join(map(self.safe_name, args))
        kwargs_repr = ', '.join('%s=%s' % (self.safe_name(k), self.safe_name(v)) for k, v in kwargs.items())
        if args_repr:
            kwargs_repr = ', ' + kwargs_repr
        return self.cache_dir / ('%s(%s%s)' % (self.func.__name__, args_repr, kwargs_repr))

    def __call__(self, *args, **kwargs):
        return self.load_or_create(*args, **kwargs)
