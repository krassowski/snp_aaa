#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
import gzip
from pyfaidx import Fasta


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

    directory = 'ensembl/v'
    species = 'Homo_sapiens'

    def __init__(self, sequence_type, id_type, version, assembly, path=None):

        self.sequence_type = sequence_type
        self.id_type = id_type
        self.version = version
        self.assembly = assembly

        # filenames of older assemblies have version included
        if int(version) <= 75:
            self.assembly += '.' + version

        if path is None:
            data = [self.species, self.assembly, self.sequence_type, self.id_type]
            path = self.directory + self.version + '/' + '.'.join(data) + '.fa.gz'

        self.path = path


class TranscriptSequenceDB(BasicSequenceDB):
    """Used for cds and cdna"""

    def __init__(self, version, assembly, sequence_type='cds', id_type='all', path=None, restrict_to=None):

        super(TranscriptSequenceDB, self).__init__(sequence_type, id_type, version, assembly, path)

        self.restrict_to = restrict_to
        print('[loading database: ' + self.path + ']')
        self.db = {}
        self.load(self.path)
        print('[database loaded]')

    def load(self, path):
        self.db = {}
        from variant_sources.ensembl import fast_gzip_read
        with fast_gzip_read(path) as f:

            transcript_id = None
            skip = False

            for line in f:
                if line.startswith('>'):
                    # skip '>' character, get first space separated string
                    transcript_id_versioned = line[1:].split(' ')[0]
                    # get rid of version number
                    transcript_id = transcript_id_versioned.split('.')[0]
                    if self.restrict_to and transcript_id in self.restrict_to:
                        skip = False
                    else:
                        skip = self.restrict_to
                    if not skip:
                        assert transcript_id not in self.db
                        self.db[transcript_id] = ''
                elif transcript_id and not skip:
                        self.db[transcript_id] += line.rstrip()
                elif not skip:
                    print('Warning: record empty at line: \n' + line)

    def fetch(self, name, strand, raw_start, raw_end, offset):
        if not raw_start or not raw_end:
            return

        seq = None

        whole_seq = self.db.get(name, None)

        if whole_seq:
            start = raw_start - 1
            end = raw_end

            cut_from = start - offset
            cut_to = end + offset

            seq = whole_seq[max(cut_from, 0):max(cut_to, 0)]

            # if we are on the edge of sequence
            if cut_from < 0:
                seq = '-' * (-cut_from) + seq
            if cut_to > len(whole_seq):
                seq += '-' * (cut_to - len(whole_seq))

            assert len(seq) == offset * 2 + raw_end - raw_start + 1

            #s = self.transcripts[name][start:end]
            #print('compare:')
            #print(seq)
            #print(s)

        return seq


class FastSequenceDB(BasicSequenceDB):
    """
    This operates on chromosome files with all sequence continuous
    without breaks or headers except the one in the first line.
    """

    def __init__(self, version, assembly, id_type, sequence_type='dna', length=60):

        super(FastSequenceDB, self).__init__(sequence_type, id_type, version, assembly)
        self.length = length

    def fetch(self, start, end, strand, offset=0):

        start -= 1

        start -= offset
        end += offset

        if start < 0:
            start = 0

        result = ''

        length = self.length

        with gzip.open(self.path, 'r') as f:
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

        # print('strand: %s' % strand)
        if strand == '-1':
            # print('inverting')
            result = complement(result)[::-1]

        return result
