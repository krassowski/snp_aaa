import os
import subprocess
import tabix
from collections import defaultdict

from tqdm import tqdm

from analyses import reporter
from analyses.spidex import spidex_get_variant, choose_record
from analyses.spidex import convert_to_strand
from expression_database import ExpressedGenes, iterate_over_expression, Gene
from expression_database import count_all, TISSUES_LIST, import_expressed_genes
from snp_parser import SPIDEX_LOCATION, jit
from variant import Variant
from scipy.stats import pearsonr
import numpy as np
import sys
import math
from IPython.core import ultratb

sys.excepthook = ultratb.FormattedTB(
    mode='Verbose',
    color_scheme='Linux',
    call_pdb=1
)


def create_path_for_genes_db(tissues):
    """Create almost-certainly unique path for a database which will contain
    information about genes having variant-egene pairs in given tissues."""

    from hashlib import sha224

    tissues_hash_code = sha224(','.join(tissues)).hexdigest()

    return 'genes_{hash_code}.db'.format(
        hash_code=tissues_hash_code
    )


@jit
def sign(x):
    return math.copysign(1, x)


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
            print(refseq, start, end)
    return cds_positions


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


def get_sequence(variant, offset):
    location = 'ensembl/v{e_version}/Homo_sapiens.GRCh{version}.{e_version}.dna.chromosome.{chrom}.fa.gz'.format(
        version=37,
        e_version=75,
        chrom=variant.chr_name
    )
    fasta = subprocess.check_output([
        'samtools',
        'faidx',
        location,
        '{chrom}:{start}-{end}'.format(
            chrom=variant.chr_name,
            start=variant.chrom_start - offset,
            end=variant.chrom_end + offset
        )
    ])
    assert fasta.startswith('>')
    return ''.join(fasta.split('\n')[1:])


def find_motifs(variants, name, sequences, min_motif_length, max_motif_length):
    location = 'motifs_discovery/'
    location += name

    if not os.path.exists(location):
        os.makedirs(location)
    location += '/'

    nearby_mutated_name = location + 'nearby.fa'
    control_name = location + 'control.fa'

    with open(nearby_mutated_name, 'w') as fn, open(control_name, 'w') as fc:
        nearby = []
        control = set()

        for variant in variants:
            full_gene_sequence = sequences[variant.refseq_transcript]
            control.add('>%s\n%s' % (variant.refseq_transcript, full_gene_sequence))

            sequence_nearby_variant = variant.sequence
            nearby.append('>%s\n%s' % (variant.refsnp_id, sequence_nearby_variant))

        fn.write('\n'.join(nearby))
        fc.write('\n'.join(list(control)))

    # TODO use meme web version for discriminate analysis with longer motifs
    args = list(map(str, [
        'dreme-py3',
        '-p', nearby_mutated_name,
        '-n', control_name,
        '-desc', name,
        '-mink', min_motif_length,
        '-maxk', max_motif_length,
        '-norc'     # "Search only the given primary sequences for motifs" - excludes searching by reverse complement of sequences
    ]))
    print(args)
    result = subprocess.check_output(args)
    return result


@reporter
def gtex_on_spidex_for_motifs(_):
    print('Gtex_on_spidex_for_motifs')

    cds_offset = 20
    min_motif_length = 10
    max_motif_length = 14

    variants = {}

    class Record(object):
        __slots__ = ('variant', 'zscore', 'valid', 'inconsistent', 'spidex_record')

        def __init__(self, variant, zscore):
            self.zscore = zscore
            self.variant = variant
            self.valid = True
            self.inconsistent = False

    for variant, tissue, slope, spidex_record, g in iterate_gtex_vs_spidex(cds_type='CDS'):

        if spidex_record.location != 'exonic':
            continue

        var_repr = variant.as_hgvs()

        slope = float(slope)
        zscore = float(spidex_record.dpsi_zscore)

        if var_repr not in variants:
            record = Record(variant, zscore)
            record.spidex_record = spidex_record
            variants[var_repr] = record

            if variant.chrom_start - g.start < cds_offset or g.end - variant.chrom_end < cds_offset:
                record.valid = False
                continue
        else:
            record = variants[var_repr]

            if not record.valid:
                continue

        # TODO filter out mutants with zscore/slope too close to zero
        if sign(zscore) != sign(slope):
            record.inconsistent = True

    variants_down = []
    variants_up = []
    variants_inconsistent = []

    transcripts = [r.variant.refseq_transcript for r in variants.values()]
    cds_positions = get_cds_positions(transcripts)

    no_cds_cnt = 0
    invalid_seq_cnt = 0
    too_close_cnt = 0

    for r in variants.values():
        ensembl_gene = r.variant.ensembl_gene_stable_id.split('.')[0]
        variant = r.variant
        # check that we are one 'offset' away from UTRs
        try:
            cds = cds_positions[variant.refseq_transcript]
        except KeyError:
            no_cds_cnt += 1
            r.valid = False
            continue

        if cds[0] > variant.chrom_start - cds_offset or cds[1] < variant.chrom_end + cds_offset:
            too_close_cnt += 1
            r.valid = False
            continue

        # Get sequence
        sequence = get_sequence(r.variant, cds_offset)
        if len(sequence) == 2 * cds_offset + 1:
            r.variant.sequence = sequence
            r.valid = True
        else:
            invalid_seq_cnt += 1
            # print('Invalid sequence:', sequence)
            r.valid = False

    print('Skipped %s variants: no cds data)' % no_cds_cnt)
    print('Skipped %s variants: sequence retrieved does not mach requested offsets)' % invalid_seq_cnt)
    print('Skipped %s variants: too close CDS edge)' % too_close_cnt)

    for r in variants.values():
        variant = r.variant
        if r.valid:
            if not r.inconsistent:
                if r.zscore > 0:
                    variants_up.append(variant)
                else:
                    variants_down.append(variant)
            else:
                variants_inconsistent.append(variant)

    sequences = load_refseq_sequences(transcripts)

    groups = {
        'consistent_up': variants_up,
        'consistent_down': variants_down,
        'inconsistent': variants_inconsistent
    }

    try:
        for group_name, variants_list in groups.items():
            out = find_motifs(variants_list, group_name, sequences, min_motif_length, max_motif_length)
            print(out)
    except Exception as e:
        import traceback
        traceback.print_exc(e)

    print('Finished, entering interactive mode')
    from IPython import embed
    embed()


@reporter
def gtex_on_spidex(_):
    print('Gtex_on_spidex')

    effect_sizes = []
    z_scores = []

    with open('gtex_vs_spidex.csv', 'w') as f:
        f.write(
            'chrom,pos,ref,alt,gene_name,gene_start,gene_end,gene_strand,'
            'gtex_slope,gtex_tissue,spidex_dpsi_zscore,spidex_dpsi_max_tissue\n'
        )
        for variant, tissue, slope, record, g in iterate_gtex_vs_spidex():
            f.write(','.join(map(str, [
                variant.chr_name, variant.chrom_start, variant.ref,
                variant.alts[0], g.name, g.start, g.end, g.strand,
                slope, tissue, record.dpsi_zscore, record.dpsi_max_tissue
            ])) + '\n')

            effect_sizes.append(float(slope))
            #z_scores.append(float(record.dpsi_max_tissue))
            z_scores.append(float(record.dpsi_zscore))

    print('Found %s pairs gtex-spidex' % len(effect_sizes))

    with open('effect_sizes.txt', 'w') as f:
        f.writelines([str(x) + '\n' for x in effect_sizes])

    with open('z_scores.txt', 'w') as f:
        f.writelines([str(x) + '\n' for x in z_scores])

    effect_sizes = np.array(effect_sizes)
    z_scores = np.array(z_scores)

    # "The null hypothesis is that the two variables are uncorrelated."
    pearson_coef, p_value = pearsonr(effect_sizes, z_scores)

    print('Pearson\'s r; p-value')
    print(pearson_coef, p_value)


def iterate_gtex_vs_spidex(**kwargs):
    """
    In spidex there are only SNPs (single!)

    Definitions for GTEx (from http://www.gtexportal.org/home/documentationPage):
        The effect size of the eQTLs is defined as the slope of the linear regression,
        and is computed as the effect of the alternative allele (ALT) relative to the
        reference allele (REF) in the human genome reference GRCh37/hg19
        (i.e., the eQTL effect allele is the ALT allele).


    Definitions for SPIDEX (more in spidex/README):
        dpsi_max_tissue: The delta PSI. This is the predicted change in
                         percent-inclusion due to the variant, reported
                         as the maximum across tissues (in percent).
        dpsi_zscore: This is the z-score of dpsi_max_tissue relative to the
                     distribution of dPSI that are due to common SNP.

        ref_allele: The reference allele at the variant position (forward-strand)
        mut_allele: The mutant allele at the variant position (forward-strand)

    """

    tissues_list = TISSUES_LIST
    # Use "Brain cortex" for basic tests - it's very small
    # tissues_list = ['Brain_Cortex']
    # Use "Adipose Subcutaneous" for larger tests
    # tissues_list = ['Adipose_Subcutaneous']

    path = create_path_for_genes_db(tissues_list)
    genes = ExpressedGenes(path)
    genes.reset()

    import_expressed_genes(
        genes,
        tissues_list=tissues_list
    )

    tb = tabix.open(SPIDEX_LOCATION)

    count = count_all(tissues_list)

    for mutation_code, tissue, slope, gene in tqdm(iterate_over_expression(tissues_list), total=count):
        chrom, pos, ref, alt, _ = mutation_code.split('_')

        # In spidex there are only SNPs (single!)
        if len(ref) != 1 or len(alt) != 1:
            continue

        pos = int(pos)

        g = Gene(*genes[gene])

        if not g:
            print('gene %s not present in data' % gene)
            continue

        # only genes overlapping with given mutation
        if not (g.start <= pos <= g.end):
            continue

        variant = Variant(
            chr_name=chrom,
            chrom_start=pos,
            chrom_end=pos,
            chrom_strand=g.strand,
            refsnp_id='-',
            ref=convert_to_strand(ref, g.strand),
            alts=(convert_to_strand(alt, g.strand),),
            ensembl_gene_stable_id=gene,
        )

        records = spidex_get_variant(tb, variant)
        records = [
            record
            for record in records
            if all(
                getattr(record, key) == value
                for key, value in kwargs.items()
            )
        ]
        record = choose_record(records, variant, variant.alts[0], strict=False)

        if record:
            if g.name != record.gene:
                print('Gene name mismatch %s %s!' % (g.name, record.gene))
                continue

            variant.refseq_transcript = record.transcript

            yield variant, tissue, slope, record, g

