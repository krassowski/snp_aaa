import gc
import itertools
import subprocess
from contextlib import contextmanager
from time import sleep

from tqdm import tqdm

from multiprocessing import Process, Manager

from cache import args_aware_cacheable


@contextmanager
def fast_gzip_read(file_name, single_thread=False):
    command = 'zcat %s' if single_thread else 'unpigz -p 4 -c %s'
    print(command % file_name)
    p = subprocess.Popen(
        (command % file_name).split(' '),
        #shell=True,
        #bufsize=500000,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT
    )
    yield p.stdout


@args_aware_cacheable
def count_lines(file_name, single_thread=False, grep=None):

    command = 'zcat %s' if single_thread else 'unpigz -p 6 -c %s'
    if grep:
        command += '| grep %s' % grep
    command += ' | wc -l'

    out = subprocess.Popen(
        command % file_name,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT
     ).communicate()[0]

    return int(out.partition(b' ')[0])


def parser(progress_updates, in_queue):
    while True:
        lines = in_queue.get()

        if lines is None:
            in_queue.task_done()
            return

        yield lines

        progress_updates.put(len(lines))


def parse_gz_file(filename, target, static_args=None, shared_args=None, chunk_size=100000, workers=6):

    def progress_bar_listener(queue):
        bar = tqdm(total=count_lines(filename))
        for step in iter(queue.get, None):
            if step is None:
                return
            bar.update(step)
            sleep(0.05)

    progress_updates = manager.Queue()

    progress_bar_process = Process(target=progress_bar_listener, args=[progress_updates])
    progress_bar_process.start()

    in_queue = manager.Queue(workers)

    if not static_args:
        static_args = []

    if not shared_args:
        shared_args = []

    pool = []
    for _ in xrange(workers):
        process = Process(
            target=target,
            args=[progress_updates, in_queue] + static_args + shared_args,
        )
        process.start()
        pool.append(process)

    print('Distributing tasks...')

    with fast_gzip_read(filename) as f:
        iterator = grouper(f, chunk_size)
        for chunk in iterator:
            in_queue.put(chunk)

    print('Work distributed')

    for _ in xrange(workers):
        in_queue.put(None)

    for p in pool:
        p.join()

    progress_updates.put(None)
    progress_bar_process.join()

    print('Joined all')
    gc.collect()

    return shared_args


def grouper(iterable, chunk_size, fill_value=None):
    args = [iter(iterable)] * chunk_size
    return itertools.izip_longest(fillvalue=fill_value, *args)


manager = Manager()