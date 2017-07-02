import os
import subprocess
from collections import defaultdict
from os.path import dirname

from tqdm import tqdm


def get_cds_positions(transcripts):
    cds_positions = {}
    with open('ucsc/hgTables') as f:
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
    handles['variant_flanked_file'].write('\n'.join(nearby))
    handles['control_unique_file'].write('\n'.join(list(set(control))))
    handles['control_redundant_file'].write('\n'.join(control))

    for handle in handles.values():
        handle.close()


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

    location = 'motifs_discovery/' + dir_name

    if not os.path.exists(location):
        os.makedirs(location)

    location += '/'

    sequences_paths = {
        'variant_flanked_file': location + 'nearby.fa',
        'control_unique_file': location + 'control.fa',
        'control_redundant_file': location + 'control_not_unique.fa'
    }

    nearby = []
    control = []

    part = None

    if slice_by or max_sequence_length:
        part = 0

    def create_handles(part):
        return {
            name: open('%s_part_%s' % (file_name, part) if part is not None else file_name, 'w')
            for name, file_name in sequences_paths.items()
        }

    handles = create_handles(part)

    for i, variant in tqdm(enumerate(variants), total=len(variants)):
        if slice_by and (i + 1) % slice_by == 0:
            dump_sequences(handles, nearby, control)
            part += 1
            handles = create_handles(part)

            nearby = []
            control = []

        full_gene_sequence = control_sequences[variant.refseq_transcript]
        control.append('>%s\n%s' % (variant.refseq_transcript, full_gene_sequence))

        sequence_nearby_variant = variant.sequence
        nearby.append('>%s_%s_%s\n%s' % (variant.chr_name, variant.chr_start, variant.alts[0], sequence_nearby_variant))

    dump_sequences(handles, nearby, control)

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

    variants_path = motifs_files['variant_flanked_file']
    control_path = motifs_files['control_unique_file']

    files_path = dirname(variants_path)
    output_path = files_path + '/' + 'dreme_out'

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

    args = list(map(str, args))
    print('Executing', ' '.join(args))
    result = subprocess.check_output(args)
    print(result)

    return result


def load_refseq_sequences(transcripts):
    db = defaultdict(str)
    with open('ucsc/sequences.fasta') as f:
        refseq = None
        for line in f:
            if line.startswith('>'):
                refseq = line[1:].split('.')[0].rstrip()
            else:
                db[refseq] += line.rstrip()
    return db