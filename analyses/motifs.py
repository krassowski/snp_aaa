import subprocess
from contextlib import contextmanager

from pathlib import Path

from pyfaidx import Fasta
from tqdm import tqdm

from multiprocess import fast_gzip_read
from helpers import take_transcript_id_without_version

refseq_sequences_fasta = 'ucsc/sequences.fasta'


def get_cds_positions(transcripts):
    cds_positions = {}
    with fast_gzip_read('ucsc/ref_gene.tsv.gz') as f:
        header = next(f)
        # assert header == '#bin    name    chrom   strand  txStart txEnd   cdsStart        cdsEnd  exonCount       exonStarts      exonEnds        score   name2   cdsStartStat cdsEndStat       exonFrames'
        for line in f:
            data = line.split('\t')
            refseq = data[1]
            if refseq not in transcripts:
                continue
            start, end = map(int, data[6:7+1])
            cds_positions[refseq] = (start, end)
    return cds_positions


def dump_sequences(handles, nearby, control):
    handles['variant_flanked'].write('\n'.join(nearby))
    handles['control_unique'].write('\n'.join(list(set(control))))
    handles['control_redundant'].write('\n'.join(control))


class PartitionedFiles:

    def __init__(self, paths, write_callback, slice_by, data_buffers):
        self.part = 0
        self.paths = paths
        self.handles = {}
        self.write_callback = write_callback
        self.slice_by = slice_by
        self.data_buffers = data_buffers

        self.open_handles()

    def dump(self):
        self.write_callback(self.handles, *self.data_buffers)

        # clean buffers
        for buffer in self.data_buffers:
            buffer.reset()

        # move on
        self.part += 1

        # open new handles
        self.close_handles()
        self.open_handles()

    def dump_if_needed(self, i):
        if self.slice_by and (i + 1) % self.slice_by == 0:
            self.dump()

    def open_handles(self):
        self.handles = {}

        for name, path in self.paths.items():
            file_name = '%s_part_%s' % (path.as_posix(), self.part)
            self.handles[name] = open(file_name, 'w')

    def close_handles(self):
        for handle in self.handles:
            handle.close()


@contextmanager
def partitioned_files(paths, slice_by, write, *args):
    files = PartitionedFiles(paths, write, slice_by, *args)
    yield files
    files.dump()
    files.close_handles()


def prepare_files_with_motifs(variants, dir_name, control_sequences, slice_by=None, max_sequence_length=None):
    """Write flanked sequences of variants and relevant control sequences to Fasta files.

    If slice_by or max_sequence_length is given (it might be desired e.g.
    in order to upload sequences to online meme service which has multiple
    restriction on input size), sequences will be partitioned into several
    files according to:
        slice_by: the maximal count of sequence in a single "part" file
        max_sequence_length: combined sequence length to adjust the number of
            sequences to be saved in a "part" file

    Returns: a dictionary pointing to locations of created files
    """

    # determine how many sequences can be put in file
    if max_sequence_length:
        some_sequence = variants[0].sequence
        slice_by = min(slice_by or 0, max_sequence_length // len(some_sequence) - 1)

    location = Path('motifs_discovery') / dir_name

    location.mkdir(exist_ok=True)

    sequences_paths = {
        'variant_flanked': location / 'nearby.fa',
        'control_unique': location / 'control.fa',
        'control_redundant': location / 'control_not_unique.fa'
    }

    nearby = []
    control = []

    with partitioned_files(sequences_paths, slice_by, dump_sequences, [nearby, control]) as files:
        for i, variant in tqdm(enumerate(variants), total=len(variants)):
            files.dump_if_needed(i)

            control.append('>%s\n%s' % (
                variant.refseq_transcript,
                control_sequences[variant.refseq_transcript]    # this is full sequence of given gene
            ))

            nearby.append('>%s_%s_%s\n%s' % (
                variant.chr_name,
                variant.chr_start,
                variant.alts[0],
                variant.sequence    # this is sequence in place of mutation +/- offset
            ))

    return sequences_paths


def find_motifs(motifs_files, min_motif_length, max_motif_length, description='', program='dreme', use_control=True):
    """Find motifs using either dreme or meme package.

    By default uses provided control sequences set if 'dreme' is selected.
    Control sequences cannot be used with meme as offline version does not support that.

    Arguments:
        motifs_files: a dict with keys pointing to variant sequences and control sequences files
        min_motif_length: mink/minw
        max_motif_length: maxk/maxw
        program: either 'meme' or 'dreme'
        use_control: should provided control sequences be used or should dreme use shuffled sequences?
    """
    if program == 'meme' and use_control:
        raise Exception('Cannot use control sequences with meme. Change use_control to False.')

    variants_path = motifs_files['variant_flanked']
    control_path = motifs_files['control_unique']

    output_path = variants_path.parent / 'dreme_out'

    if program == 'dreme':
        args = (
            [
                'dreme-py3',
                '-p', variants_path
            ]
            + (['-n', control_path] if use_control else []) +
            [
                '-desc', description,
                '-mink', min_motif_length,
                '-maxk', max_motif_length,
                '-oc', output_path,
                '-norc'     # "Search only the given primary sequences for motifs":
                #  excludes searching by reverse complement of sequences
            ]
        )

    elif program == 'meme':
        args = [
            'meme',
            variants_path,
            '-desc', description,
            '-minw', min_motif_length,
            '-maxw', max_motif_length,
            '-nmotifs', 15,
            '-oc', output_path,
            # '-revcomp' is not deactivated by default (complementary to norc)
        ]
    else:
        raise ValueError('Unsupported program: %s' % program)

    args = list(map(str, args))
    print('Executing', ' '.join(args))
    result = subprocess.check_output(args)
    print(result)

    return result


def load_refseq_sequences():
    return Fasta(refseq_sequences_fasta, key_function=take_transcript_id_without_version)
