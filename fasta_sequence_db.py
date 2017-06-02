#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
import gzip


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
