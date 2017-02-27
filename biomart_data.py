from biomart import BiomartDataset
from data_store import DataStore


class BiomartSequence(object):

    __slots__ = 'header', 'sequence'

    def __init__(self, entry, attributes):

        entry = entry.split('\n')
        self.header = entry[0].lstrip('>')
        self.sequence = entry[1:]


class BiomartData(DataStore):

    def __init__(self, dataset, data_raw, filters):

        assert isinstance(dataset, BiomartDataset)

        self.attributes = data_raw.attributes
        self.filters = filters
        self.dataset = dataset

        # Either I do not understand how to force biomart
        # to give the header data or it does not work
        # response = self.fetch_data(header=1)
        # print(response.text)
        data = self.fetch_data()
        iterator = data.iter_lines()

        super(BiomartData, self).__init__(data_raw, data, iterator)

    def fetch_data(self, header=0):
        return self.dataset.search({
            'filters': self.filters,
            'attributes': [unicode(x) for x in self.attributes]
            }, formatter='TSV', header=header)

    def count(self):
        return self.dataset.count()
