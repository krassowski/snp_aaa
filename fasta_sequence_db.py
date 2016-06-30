#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function


def complement(seq):
    basic = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    try:
        return ''.join([basic[n] for n in seq])
    except KeyError:
        # http://arep.med.harvard.edu/labgc/adnan/projects/Utilities/revcomp.html
        IUPAC = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'U': 'A', 'Y': 'R',
                 'R': 'Y', 'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K', 'B': 'V',
                 'V': 'B', 'D': 'H', 'H': 'D', 'N': 'N'}
        return ''.join([IUPAC[n] for n in seq])


class BasicSequenceDB(object):

    directory = 'ensembl'
    assembly = 'GRCh38'
    species = 'Homo_sapiens'

    def __init__(self, sequence_type, id_type, path=None):

        self.sequence_type = sequence_type
        self.id_type = id_type

        if path is None:
            data = [self.species, self.assembly, self.sequence_type, self.id_type]
            path = self.directory + '/' + '.'.join(data) + '.fa'

        self.path = path

    def parse_coordinates(self, s, e):
        start = int(s)
        end = int(e)
        if self.sequence_type in ['cds', 'cdna']:
            start -= 1
        elif self.sequence_type == 'dna':
            end += 1
        else:
            raise Exception('Unknown seq_type')
        return start, end


class SequenceDB(BasicSequenceDB):

    def __init__(self, index_by='transcript', sequence_type='cds', id_type='all', path=None, restrict_to=None):

        super(SequenceDB, self).__init__(sequence_type, id_type, path)

        self.restrict_to = restrict_to
        print('[loading database: ' + self.path + ']')
        self.load(index_by, self.path)
        print('[database loaded]')

    def load(self, index_by, path):
        with open(path, 'r') as f:
            if index_by == 'transcript':
                self.load_by_transcript(f)
            else:
                raise Exception('Unknown index_by')

    def load_by_transcript(self, f):
        self.db = {}
        transcript_id = None
        skip = bool(self.restrict_to)
        for line in f:
            if line.startswith('>'):
                # skip '>' character, get first space separated string
                transcript_id_versioned = line[1:].split(' ')[0]
                # get rid of version number
                transcript_id = transcript_id_versioned.split('.')[0]
                if self.restrict_to and transcript_id in self.restrict_to:
                    skip = False
                else:
                    skip = bool(self.restrict_to)
                if not skip:
                    assert transcript_id not in self.db
                    self.db[transcript_id] = ''
            elif transcript_id and not skip:
                    self.db[transcript_id] += line.rstrip()
            elif not skip:
                print('Warning: record empty at line: \n' + line)

    def has(self, seq_id):
        return seq_id in self.db

    def get(self, name):
        return self.db.get(name, '')

    def fetch(self, name, strand, raw_start, raw_end, offset):
        if not raw_start or not raw_end:
            return None
        seq = None
        if raw_start and raw_end and self.has(name):
            start, end = self.parse_coordinates(raw_start, raw_end)
            whole_seq = self.get(name)

            if strand == -1:
                whole_seq = complement(whole_seq)[::-1]
                start, end = len(whole_seq) - end, len(whole_seq) - start
            cut_from = start - offset
            cut_to = end + offset
            seq = whole_seq[cut_from if cut_from >= 0 else 0:cut_to if cut_to >=0 else 0]
            # if we are on the edge of sequence
            if cut_from < 0:
                seq = '-' * (-cut_from) + seq
            if cut_to > len(whole_seq):
                seq += '-' * (cut_to - len(whole_seq))

            assert len(seq) == offset * 2 + int(raw_end) - int(raw_start) + 1

        return seq


class FastSequenceDB(BasicSequenceDB):

    """
    This operates on chromosome files with all sequence continuous
    without breakes or headers except the one in the first line.
    """

    def __init__(self, id_type, sequence_type='dna', length=60):

        super(FastSequenceDB, self).__init__(sequence_type, id_type)
        self.length = length

    def fetch(self, start, end, offset=0):

        start -= 1

        start -= offset
        end += offset

        if start < 0:
            start = 0

        result = ''

        length = self.length

        with open(self.path, 'r') as f:
            header = f.readline()

            first_line = start // length

            position = first_line * length

            # (length + 1) stands for length of line with new line character
            pos_in_file = first_line * (length + 1) + len(header)

            f.seek(pos_in_file)

            for line in f:

                position += length

                if start <= position:

                    pos_from = start - position + length
                    pos_to = end - position + length

                    if pos_from < 0:
                        pos_from = 0

                    if pos_to > length:
                        pos_to = length

                    result += line.rstrip()[pos_from:pos_to]

                    if pos_to < length:
                        break

        return result

