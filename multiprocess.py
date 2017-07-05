import gc
import itertools
import subprocess
from contextlib import contextmanager
from time import sleep

from tqdm import tqdm

from cache import args_aware_cacheable

from multiprocessing import get_context

ctx = get_context('spawn')

@contextmanager
def fast_gzip_read(file_name, processes=4):
    command = 'zcat %s' if processes == 1 else 'unpigz -p ' + str(processes) + ' -c %s'
    p = subprocess.Popen(
        (command % file_name).split(' '),
        #shell=True,
        #bufsize=500000,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        universal_newlines=True
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


def get_manager(manager_store=[]):
    if not manager_store:
        manager = ctx.Manager()
        manager_store.append(manager)
    return manager_store[0]


def progress_bar_listener(filename, queue):
    bar = tqdm(total=count_lines(filename))
    for step in iter(queue.get, None):
        if step is None:
            return
        bar.update(step)
        sleep(0.05)


def parse_gz_file(filename, target, static_args=None, shared_args=None, chunk_size=100000, workers=6):

    manager = get_manager()

    progress_updates = manager.Queue()

    # clean up before fork
    gc.collect()

    progress_bar_process = ctx.Process(target=progress_bar_listener, args=[filename, progress_updates])
    progress_bar_process.start()

    in_queue = manager.Queue(workers)

    if not static_args:
        static_args = []

    if not shared_args:
        shared_args = []

    pool = []
    for _ in range(workers):
        process = ctx.Process(
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

    for _ in range(workers):
        in_queue.put(None)

    for p in pool:
        p.join()

    progress_updates.put(None)
    progress_bar_process.join()

    print('Joined all')

    # clean up after work
    gc.collect()

    return shared_args


def grouper(iterable, chunk_size, fill_value=None):
    args = [iter(iterable)] * chunk_size
    return itertools.zip_longest(fillvalue=fill_value, *args)


