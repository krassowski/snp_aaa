import datetime
import time


def decorator_maker(storage, decorator_name):

    def decorator(func):
        name = func.__name__

        def informative_reporter(*args, **kwargs):
            print('Executing %s: %s...' % (decorator_name, name))
            with execution_time() as time:
                result = func(*args, **kwargs)
            print(
                'Execution of %s %s finished after %s'
                %
                (name, decorator_name, time.elapsed)
            )
            return result

        storage[name] = informative_reporter
        return informative_reporter

    return decorator


def select_poly_a_related_variants(variants):
    """Return list of variants occurring in a poly(A) track and variants which will create poly(A) track."""

    return [
        variant
        for variant in variants
        if any([
            data.has or data.will_have
            for affected_transcript in variant.affected_transcripts
            for data in affected_transcript.poly_aaa.values()
        ])
    ]


def all_poly_a_variants(variants):
    poly_a_related_variants = select_poly_a_related_variants(variants.values())

    for variant in poly_a_related_variants:
        yield variant


class ExecutionTime:

    def __init__(self, start):
        self.start = start
        self.frozen = False

    @property
    def elapsed(self):
        if self.frozen:
            return self.frozen
        now = time.time()
        return datetime.timedelta(seconds=now - self.start)

    def freeze(self):
        self.frozen = self.elapsed


class execution_time:

    def __enter__(self):
        start = time.time()
        self.time = ExecutionTime(start)
        return self.time

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.time.freeze()
