from __future__ import print_function
import cPickle as pickle


class CachingError(Exception):
    pass


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
                            'Oryginal exception: %s.' %
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
