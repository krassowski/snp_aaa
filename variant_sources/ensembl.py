import gzip
import json
import subprocess
from collections import defaultdict, OrderedDict
import gc

from tqdm import tqdm

from cache import cacheable
from commands import SourceSubparser
from variant import Variant, BiomartVariant, Transcript
from variant_sources import variants_getter
from variant_sources.biomart import gene_names_from_patacsdb_csv


def count_lines(file_object):
    out = subprocess.Popen(
        'zcat ' + file_object.filename + ' | wc -l',
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT
     ).communicate()[0]
    return int(out.partition(b' ')[0])


ensembl_args = SourceSubparser(
    'ensembl',
    help='Arguments for biomart arguments fetching'
)

ensembl_args.add_command(
    '--filters',
    type=json.loads,
    default={},
    help='Additional variants filters to be used when querying biomart'
)

ensembl_args.add_command(
    '--genes_list',
    nargs='+',
    type=str,
    default=gene_names_from_patacsdb_csv(),
    help=(
             'list of human genes from which variants should be extracted for use in '
             'available analyses. By default all human genes from PATACSDB '
             'will be used. To use all human genes specify "all_human_genes".'
         )
)


@variants_getter
def ensembl(args):
    """Download variants from given --biomart, belonging to
    coding sequences of genes specified with --genes_list and
    filtered with --filters.
    """
    raw_variants_by_gene = load_ensembl_variants.load_or_create(
        args.genes_list,
        filters=args.filters
    )
    return raw_variants_by_gene


@cacheable
def load_ensembl_variants(gene_names, filters={}):

    types = OrderedDict(
        (
            ('transcript_variation_id', int),
            ('variation_feature_id', int),
            ('feature_stable_id', str),
            ('allele_string', str),
            ('somatic', bool),
            ('consequence_types', str),
            ('cds_start', int),
            ('cds_end', int),
            ('cdna_start', int),
            ('cdna_end', int),
            ('translation_start', int),
            ('translation_end', int),
            ('distance_to_transcript', int),
            ('codon_allele_string', str),
            ('pep_allele_string', str),
            ('hgvs_genomic', str),
            ('hgvs_transcript', str),
            ('hgvs_protein', str),
            ('polyphen_prediction', str),
            ('polyphen_score', float),
            ('sift_prediction', str),
            ('sift_score', float),
            ('display', bool)
        )
    )
    loc = 'ensembl/v88/GRCh37/variation_database/'
    #loc = 'ensembl/v88/GRCh37/variation_database/head/'

    by_id = {}
    keys = types.keys()
    somatic_pos = keys.index('somatic')
    consequence_pos = keys.index('consequence_types')

    from snp_parser import jit

    @jit
    def int_or_none(x):
        if x == '\\N':
            return None
        return int(x)

    transcript_strand = {}
    with gzip.open(loc + 'transcript.txt.gz') as f:
        for line in tqdm(f, total=count_lines(f)):
            data = line.split('\t')
            transcript_strand[data[14]] = data[6]

    @jit
    def check_and_split(line):
            return data, ref

    with gzip.open(loc + 'transcript_variation.txt.gz') as f:
        for line in tqdm(f, total=count_lines(f)):
            data = line.split('\t')
            if data[somatic_pos] == '1' and 'coding_sequence_variant' in data[consequence_pos].split(','):
                alleles = data[3].split('/')
                ref = alleles[0]
                if len(ref) != 1:
                    continue
                if len(alleles) > 1 and alleles[1] != 'COSMIC_MUTATION':
                    if len(alleles[1]) != 1:
                        continue
                    # should i keep it?
                    if all(allele != 'A' for allele in alleles):
                        continue
                variant_id = int(data[1])

                transcript = Transcript(
                    cds_start=int_or_none(data[6]),
                    cds_end=int_or_none(data[7]),
                    cdna_start=int_or_none(data[8]),
                    cdna_end=int_or_none(data[9]),
                    strand=transcript_strand[data[2]],
                    ensembl_id=data[2]
                )

                if variant_id in by_id:
                    v = by_id[variant_id]

                    if v.allele_1 != ref:
                        print(
                            'Reference does not match: %s vs %s'
                            'for variant %s, with transcripts: %s %s'
                            %
                            (
                                v.allele_1,
                                ref,
                                variant_id,
                                data[2],
                                list(v.affected_transcripts)[0].ensembl_id
                            )
                        )
                        continue

                    v.affected_transcripts.add(
                        transcript
                    )
                    # if v.ensembl_transcript_stable_id == data[2]:
                    #    v
                else:
                    v = Variant(
                        allele_1=ref,
                        affected_transcripts={
                            transcript
                        }
                    )
                    by_id[variant_id] = v
    gc.collect()

    print('Accepted:', len(by_id))

    seq_region = {}
    # seq_region_id, name, cord_system_fk
    with gzip.open(loc + 'seq_region.txt.gz') as f:
        for line in tqdm(f, total=count_lines(f)):
            data = line.split('\t')
            seq_region[int(data[0])] = data[1]

    sources = {}
    # 'source_id', 'name', 'version', 'description', 'url', 'type', 'somatic_status', 'data_types'
    with gzip.open(loc + 'source.txt.gz') as f:
        for line in tqdm(f, total=count_lines(f)):
            data = line.split('\t')
            sources[int(data[0])] = data[1]
    gc.collect()

    feature_table = [
        'variation_feature_id',
        'seq_region_id',
        'seq_region_start',
        'seq_region_end',
        'seq_region_strand',
        'variation_id',
        'allele_string',
        'variation_name',
        'map_weight',
        'flags',
        'source_id',
        'consequence_types',
        'variation_set_id',
        'class_attrib_id',
        'somatic',
        'minor_allele',
        'minor_allele_freq',
        'minor_allele_count',
        'alignment_quality',
        'evidence_attribs',
        'clinical_significance',
        'display',
    ]

    with gzip.open(loc + 'variation_feature.txt.gz') as f:
        for line in tqdm(f, total=count_lines(f)):
            data = line.split('\t')
            if int(data[0]) in by_id:
                v = by_id[int(data[0])]
                v.chr_name = seq_region[int(data[1])]
                v.chrom_start = int(data[2])
                v.chrom_end = int(data[3])
                v.chrom_strand = data[4]
                v.refsnp_id = data[7]
                v.refsnp_source = sources[int(data[10])]
                # if v.refsnp_source != 'COSMIC':
                #    print(v)

    gc.collect()
    by_name = defaultdict(list)
    no_name = 0
    no_u = 0
    for v in tqdm(by_id.itervalues(), total=len(by_id)):
        if not hasattr(v, 'refsnp_id'):
            no_name += 1
            continue
        if v.refsnp_id in by_name:
            by_name[v.refsnp_id].affected_transcripts.update(v.affected_transcripts)
            no_u += 1
            #print('Compare:')
            #print(by_name[v.refsnp_id])
            #print(v)
        else:
            by_name[v.refsnp_id] = v
            # transcript = v.extract_transcript()
    print('Noname: %s' % no_name)
    print('Nouniqe names: %s' % no_u)

    #del by_id
    gc.collect()
    #print('Press enter to continue')
    #x = raw_input()

    return {v_id: [v] for v_id, v in by_name.iteritems()}
    # THERE ARE MANY GENES FOR A SINGLE VARIANT
