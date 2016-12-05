from __future__ import print_function
import cPickle as pickle


def cached(action='load'):
    """Makes caching a pleasure. Use as a decorator."""
    def decorator(generating_function):
        def we_all_like_nested_closures(*args, **kwargs):
            display_name = generating_function.__name__.replace('_', ' ').title()
            cache_name = '.cache:' + generating_function.__name__

            if action == 'load':
                with open(cache_name, 'rb') as f:
                    variable = pickle.load(f)
                print('"' + display_name + '" data loaded from cache')
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
