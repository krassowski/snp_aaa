import os
import sys
biomart_fork_path = os.path.realpath(os.path.join(os.curdir, 'biomart'))
sys.path.insert(0, biomart_fork_path)
from biomart import BiomartDataset
from data_store import DataStore


class BiomartSequence(object):

    __slots__ = 'header', 'sequence'

    def __init__(self, entry, attributes):

        entry = entry.split('\n')
        self.header = entry[0].lstrip('>')
        self.sequence = entry[1:]


class BiomartData(DataStore):

    def __init__(self, dataset, attributes, filters, fasta=False):

        assert isinstance(dataset, BiomartDataset)

        self.attributes = attributes
        self.filters = filters
        self.dataset = dataset
        self.fasta = fasta

        # Either I do not understand how to force biomart
        # to give the header data or it does not work
        # response = self.fetch_data(header=1)
        # print(response.text)
        data = self.fetch_data()
        iterator = data.iter_lines()

        super(BiomartData, self).__init__(attributes, data, iterator, single_rows=not fasta)

    def parse_multi(self, entry):
        return BiomartSequence(entry, self.attributes)

    def is_ready_to_flush(self, rows):
        if self.fasta:
            return filter(bool, [row.startswith('>') for row in rows])
        else:
            return True

    def fetch_data(self, header=0):
        return self.dataset.search({
            'filters': self.filters,
            'attributes': [unicode(x) for x in self.attributes]
            }, formatter='FASTA' if self.fasta else 'TSV', header=header)

    def count(self):
        return self.dataset.count()