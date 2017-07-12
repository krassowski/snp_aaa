from collections import Counter, defaultdict

from math import log
from networkx import nx, nx_agraph
from recordclass import recordclass
from analyses import report, reporter
from helpers import all_poly_a_variants
import subprocess
import multiprocess
import matplotlib.pyplot as plt


def vep(input_filename, use_ensembl_online=True, overwrite=False):
    input_filename = 'reports/' + input_filename
    output = input_filename + '.vep'
    cmd = (
        'ensembl-vep/vep --appris --biotype --check_existing --fork 8 --port 3337 '
        '--polyphen b --regulatory --sift b --species homo_sapiens --symbol --tsl '
        '--coding_only --tab --assembly GRCh37 --format hgvs --input_file %s -o %s '
    ) % (input_filename, output)
    if overwrite:
        cmd += ' --force_overwrite'
    if use_ensembl_online:
        cmd += ' --database --db_version 88'
    else:
        cmd += ' --cache --offline --fasta Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz'

    if not os.exsists(output):
        p = subprocess.Popen(
            cmd.split(' '),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            universal_newlines=True
        )
        for line in p.stdout:
            print(line)

    return output


def load_sequence_ontology():
    import obonet
    url = 'https://raw.githubusercontent.com/The-Sequence-Ontology/SO-Ontologies/master/so.obo'
    graph = obonet.read_obo(url)
    return graph


@reporter
def summarize_poly_aaa_variants(variants):

    columns = [
        'snp_id', 'gene', 'poly_aaa_increase',
        'poly_aaa_decrease', 'poly_aaa_change',
        'chr', 'start', 'end', 'ref', 'alt',
        'transcript', 'cds_start', 'cds_end'
    ]

    Record = recordclass('RecordPolyA', columns)

    aaa_records = []
    aaa_variants = set()
    up_variants = {}
    down_variants = {}
    mutations_in_cds_hgvs_format = defaultdict(list)
    all_variants_ids = []
    variants_sources = Counter()
    transcripts = set()

    for variant in all_poly_a_variants(variants, preserve_sources=True):

        all_variants_ids.extend(variant.snp_id.split(','))

        for transcript in variant.affected_transcripts:

            if not transcript.poly_aaa:
                continue

            for alt, aaa_data in transcript.poly_aaa.items():

                if aaa_data.increased:
                    category = 'increased'
                elif aaa_data.decreased:
                    category = 'decreased'
                else:
                    category = 'constant'

                mutations_in_cds_hgvs_format[category].append(
                    transcript.as_hgvs(variant.ref, alt)
                )

                record = Record(
                    variant.snp_id,
                    None,#variant.ensembl_gene_stable_id    # TODO
                    aaa_data.increased,
                    aaa_data.decreased,
                    aaa_data.change,
                    variant.chr_name,
                    variant.chr_start,
                    variant.chr_end,
                    variant.ref,
                    alt,
                    transcript.ensembl_id,
                    transcript.cds_start,
                    transcript.cds_end
                )

                if aaa_data.increased:
                    up_variants[variant] = True
                if aaa_data.decreased:
                    down_variants[variant] = True
                transcripts.add(transcript.ensembl_id)

                aaa_records.append(record)

            aaa_variants.add(variant)

        for source in set(variant.source.split(',')):
            variants_sources[source] += 1

    report(
        'poly aaa increase and decrease by variants',
        aaa_records,
        columns
    )
    report(
        'poly aaa sources',
        variants_sources.items(),
        ['source', 'count']
    )

    for category, muts in mutations_in_cds_hgvs_format.items():
        report(
            'Mutations which result in ' + category + ' in cds hgvs formats',
            muts
        )

    report(
        'all ids',
        all_variants_ids
    )

    print('Affected transcripts: %s' % len(transcripts))
    print('Down variants: %s' % len(down_variants))
    print('Up variants: %s' % len(up_variants))
    print('Unique variants: %s' % len(aaa_variants))
    print('Variants identifiers: %s' % sum(v.snp_id.count(',') + 1 for v in aaa_variants))
    print(variants_sources)


def do_work(progress, in_queue, static_args, consequences):

    transcripts_to_check, ids_to_check = static_args

    for lines in multiprocess.parser(progress, in_queue):
        for line in lines:

            if line is None:
                in_queue.task_done()
                return

            data = line.split('\t')

            # ignore 'LRG' transcripts at all
            if data[2].startswith('LRG_'):
                continue

            ensembl_id = data[2]
            if ensembl_id not in transcripts_to_check:
                continue

            variant_id = data[7]
            print(variant_id)

            if variant_id not in ids_to_check:
                continue

            consequence = data[5]

            print(consequence)
            consequences.append((variant_id, ensembl_id, consequence))


@reporter
def poly_aaa_consequences(variants):

    mutations_in_cds_hgvs_format = defaultdict(list)

    for variant in all_poly_a_variants(variants, preserve_sources=True):

        for transcript in variant.affected_transcripts:

            if not transcript.poly_aaa:
                continue

            for alt, aaa_data in transcript.poly_aaa.items():

                if aaa_data.increased:
                    category = 'increased'
                elif aaa_data.decreased:
                    category = 'decreased'
                else:
                    category = 'constant'

                mutations_in_cds_hgvs_format[category].append(
                    transcript.as_hgvs(variant.ref, alt)
                )

    consequences = defaultdict(Counter)
    skipped = Counter()
    for category, muts in mutations_in_cds_hgvs_format.items():
        filename = report(
            'Mutations which result in ' + category + ' in cds hgvs formats',
            muts
        )
        vep_filename = vep(filename)
        for line in open(vep_filename):
            if line.startswith('#'):
                continue
            line = line.split('\t')
            tested_transcript = line[0].split(':')[0]
            vep_transcript = line[4]
            if line[5] != 'Transcript':
                skipped['Not a transcript feature'] += 1
                continue
            if tested_transcript != vep_transcript:
                skipped['Different transcript'] += 1
                continue

            variant_consequences = line[6].split(',')
            for consequence in variant_consequences:
                consequences[category][consequence] += 1

        print(skipped)
        print('Raw consequences')
        print(consequences)

    expanded_consequences = defaultdict(Counter)
    # expand consequences
    graph = load_sequence_ontology()

    def match_by_name(name, negate=False):
        if negate:
            def test(node_tuple):
                node, data = node_tuple
                return data['name'] != name
        else:
            def test(node_tuple):
                node, data = node_tuple
                return data['name'] == name
        return test

    def match_not_in(set_in):
        def test_in(node_tuple):
            node, data = node_tuple
            return node not in set_in
        return test_in

    def make_node_tuple(node):
        return node, graph.node[node]

    stop_on = [
        'coding_sequence_variant',
        'feature_variant',
        'internal_feature_elongation'
    ]

    for category, category_consequences in consequences.items():
        for consequence, count in category_consequences.items():
            visited = set()
            nodes = list(filter(match_by_name(consequence), graph.nodes(data=True)))
            assert len(nodes) == 1
            next_nodes = True
            while next_nodes:
                next_nodes = []
                for node in nodes:
                    node, data = node
                    expanded_consequences[category][data['name']] += count
                    visited.add(node)
                    if data['name'] in stop_on:
                        continue

                    candidate_nodes = [dest for this_node, dest in graph.out_edges(node)]
                    candidate_nodes = [make_node_tuple(node) for node in candidate_nodes]
                    candidate_nodes = list(filter(match_not_in(visited), candidate_nodes))
                    for node_tuple in candidate_nodes:
                        node, data = node_tuple
                        if node not in [n for n, d in next_nodes]:
                            next_nodes.append(node_tuple)
                nodes = next_nodes

    for category, counts in expanded_consequences.items():

        consequences_to_include = ['coding_sequence_variant']
        consequences_to_include.extend(counts.keys())
        g = graph.subgraph([node for node, data in graph.nodes(data=True) if data['name'] in consequences_to_include])
        g = g.reverse()

        max_count = max(counts.values())

        for node, data in g.nodes_iter(data=True):
            name = data['name']
            count = counts[name]
            color = (255 - int(log((count / max_count) + 1) * 255), 255, 255)
            g.node[node]['style'] = 'filled'
            g.node[node]['shape'] = 'box'
            color = '#%02x%02x%02x' % color
            g.node[node]['fillcolor'] = color
            if name not in consequences[category]:
                g.node[node]['style'] = 'dashed,filled'

        g = nx.relabel_nodes(g, {node: data['name'] + ': %s' % counts.get(data['name']) for node, data in g.nodes(data=True)})

        a = nx_agraph.to_agraph(g)

        a.layout('dot', args='-Nfontsize=10 -Nwidth=".2" -Nheight=".2" -Nmargin=.1 -Gfontsize=8 -Earrowsize=.5')
        a.draw('reports/poly_a_consequences_dag_' + category + '.svg')

        plt.show()

    print('Expanded consequences:')
    print(expanded_consequences)
