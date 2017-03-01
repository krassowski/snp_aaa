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
        print('Generating "' + self.display_name + '"...')
        result = self.func(*args, **kwargs)
        print('"' + self.display_name + '" saved to cache')

        with open(self.cache_name, 'wb') as f:
            pickle.dump(result, f, protocol=pickle.HIGHEST_PROTOCOL)

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
            except TypeError:
                print(
                    'Whops! %s needs to be loaded manually as there are few arguments required' %
                    self.display_name
                )
                sys.exit(1)

    def load(self):
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


def cached(action='load'):
    """Makes caching a pleasure. Use as a decorator."""
    def decorator(generating_function):
        def we_all_like_nested_closures(*args, **kwargs):
            display_name = generating_function.__name__.replace('_', ' ').title()
            cache_name = '.cache:' + generating_function.__name__

            if action == 'load':
                with open(cache_name, 'rb') as f:
                    try:
                        variable = pickle.load(f)
                    except AttributeError as e:
                        raise CachingError(
                            'Unable to load cache for %s. '
                            'Please, make sure that the cache files correspond'
                            ' to this version of application.\n'
                            'Original exception: %s.' %
                            (display_name, e.message)
                        )
                # print('"' + display_name + '" data loaded from cache')
            else:
                print('Generating "' + display_name + '"...')
                variable = generating_function(*args, **kwargs)
                if action == 'save':
                    print('"' + display_name + '" saved to cache')
                    with open(cache_name, 'wb') as f:
                        pickle.dump(variable, f, protocol=pickle.HIGHEST_PROTOCOL)
            return variable
        return we_all_like_nested_closures
    return decorator
