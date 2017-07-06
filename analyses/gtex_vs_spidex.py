import math
import tabix
from collections import Counter

import numpy as np
import subprocess
from scipy.stats import pearsonr, spearmanr
from tqdm import tqdm

import settings
from analyses import reporter, report
from analyses.motifs import get_cds_positions, prepare_files_with_motifs, find_motifs, load_refseq_sequences
from analyses.spidex import convert_to_strand
from exceptions import StrandMismatch, Intronic
from analyses.spidex import spidex_get_variant, choose_record
from cache import cacheable
from expression_database import ExpressedGenes, iterate_over_expression
from expression_database import count_all, GTEX_TISSUES, import_expressed_genes
from jit import jit
from settings import SPIDEX_LOCATION
from models import SingleAltVariant, Gene


def create_path_for_genes_db(tissues):
    """Create almost-certainly unique path for a database which will contain
    information about genes having variant-egene pairs in given tissues."""

    from hashlib import sha224

    tissues_serialized = ','.join(tissues).encode('utf-8')
    tissues_hash_code = sha224(tissues_serialized).hexdigest()

    return 'genes_{hash_code}.db'.format(
        hash_code=tissues_hash_code
    )


@jit
def sign(x):
    return math.copysign(1, x)


def get_sequence(variant, offset):
    location = 'ensembl/v{e_version}/Homo_sapiens.{version}.{e_version}.dna.chromosome.{chrom}.fa.gz.bgz'.format(
        version=settings.GRCH_VERSION,
        e_version=settings.ENSEMBL_VERSION,
        chrom=variant.chr_name
    )
    fasta = subprocess.check_output([
        'samtools',
        'faidx',
        location,
        '{chrom}:{start}-{end}'.format(
            chrom=variant.chr_name,
            start=variant.chr_start - offset,
            end=variant.chr_end + offset
        )
    ])
    assert fasta.startswith('>')
    return ''.join(fasta.split('\n')[1:])


@cacheable
def get_muts_groups_and_seqs(cds_offset):
    """
    When searing for motifs using unchanged sequence, I assume that there
    are more mutations destroying already existing motifs.

    Conversely, if I use mutated sequences, I expect that more motifs
    has been created and the new motifs occurrence has an important effect
    on expression regulation.

    The second seems to be more difficult to detect as the entropy-derived
    probabilities of such events are unfavourable.
    """
    variants = {}

    class Record:
        __slots__ = ('variant', 'zscore', 'g', 'valid', 'inconsistent', 'spidex_record')

        def __init__(self, variant, zscore, gene):
            self.g = gene
            self.zscore = zscore
            self.variant = variant
            self.valid = True
            self.inconsistent = False

    for variant, tissue, slope, spidex_record, g in iterate_gtex_vs_spidex(cds_type='CDS'):

        if spidex_record.location != 'exonic':
            continue

        key = (variant.chr_name, variant.chr_start, variant.alt, g.name)

        slope = float(slope)
        zscore = float(spidex_record.dpsi_zscore)

        if key not in variants:
            record = Record(variant, zscore, g)
            record.spidex_record = spidex_record
            variants[key] = record
            if variant.chr_start - g.start < cds_offset or g.end - variant.chr_end < cds_offset:
                record.valid = False
                continue
        else:
            record = variants[key]
            # most of the time, if there is a key, the validity of the variant
            # will be already known. But in if executing with multithreading it
            # can be not yet computed; therefore after 'if not valid' condition
            # there is a second check: is it really valid?
            if not record.valid:
                continue

            if variant.chr_start - g.start < cds_offset or g.end - variant.chr_end < cds_offset:
                record.valid = False
                continue

        # TODO filter out mutants with zscore/slope too close to zero
        if sign(zscore) != sign(slope):
            record.inconsistent = True

    print('Grouping spidex-gtex pairs by effect')

    variants_down = []
    variants_up = []
    variants_inconsistent = []

    transcripts = [r.variant.refseq_transcript for r in variants.values()]
    cds_positions = get_cds_positions(transcripts)

    no_cds_cnt = 0
    invalid_seq_cnt = 0
    too_close_cnt = 0

    for r in variants.values():
        ensembl_gene = r.variant.gene.split('.')[0]
        variant = r.variant
        # check that we are one 'offset' away from UTRs
        try:
            cds = cds_positions[variant.refseq_transcript]
        except KeyError:
            no_cds_cnt += 1
            r.valid = False
            continue

        if cds[0] > variant.chr_start - cds_offset or cds[1] < variant.chr_end + cds_offset:
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
    return groups, sequences


@reporter
def gtex_on_spidex_motifs_dreme(_):
    cds_offset = 50
    min_motif_length = 8
    max_motif_length = 14

    groups, sequences = get_muts_groups_and_seqs.load_or_create(cds_offset)

    for group_name, variants_list in groups.items():
        motifs_files = prepare_files_with_motifs(variants_list, group_name, sequences)

        out = find_motifs(
            motifs_files,
            min_motif_length,
            max_motif_length,
            description=group_name,
            program='dreme',
            use_control=False
        )
        print(out)


@reporter
def gtex_on_spidex_motifs_dreme_with_control(_):
    cds_offset = 50
    min_motif_length = 8
    max_motif_length = 14

    groups, sequences = get_muts_groups_and_seqs.load_or_create(cds_offset)

    for group_name, variants_list in groups.items():
        motifs_files = prepare_files_with_motifs(variants_list, group_name, sequences)

        out = find_motifs(
            motifs_files,
            min_motif_length,
            max_motif_length,
            description=group_name,
            program='dreme',
            use_control=True
        )
        print(out)


@reporter
def gtex_on_spidex_motifs_meme_online(_):
    cds_offset = 50

    groups, sequences = get_muts_groups_and_seqs.load_or_create(cds_offset)

    for group_name, variants_list in groups.items():
        motifs_files = prepare_files_with_motifs(
            variants_list,
            group_name,
            sequences,
            slice_by=1000,
            max_sequence_length=60000
        )
        print('Motifs files to upload online written to part files:')
        print(motifs_files)


@reporter
def gtex_on_spidex(_):
    effect_sizes = []
    z_scores = []
    report_data = []

    for variant, tissue, slope, record, g in iterate_gtex_vs_spidex():
        report_data.append((
            variant.chr_name, variant.chr_start, variant.ref,
            variant.alt, g.name, g.start, g.end, g.strand,
            slope, tissue, record.dpsi_zscore, record.dpsi_max_tissue
        ))

        effect_sizes.append(float(slope))
        #z_scores.append(float(record.dpsi_max_tissue))
        z_scores.append(float(record.dpsi_zscore))  # this is based on dpsi_max_tissues!

    report(
        'gtex_vs_spidex',
        report_data,
        'chrom,pos,ref,alt,gene_name,gene_start,gene_end,gene_strand,'
        'gtex_slope,gtex_tissue,spidex_dpsi_zscore,spidex_dpsi_max_tissue'.split(',')
    )

    print('Found %s pairs gtex-spidex' % len(effect_sizes))

    report('effect_sizes', effect_sizes)
    report('z_scores', z_scores)

    effect_sizes = np.array(effect_sizes)
    z_scores = np.array(z_scores)

    # "The null hypothesis is that the two variables are uncorrelated."
    pearson_coef, p_value = pearsonr(effect_sizes, z_scores)

    print('Pearson\'s r; p-value')
    print(pearson_coef, p_value)


@reporter
def same_tissues_gtex_on_spidex(_):
    effect_sizes = {}
    z_scores = {}

    spidex_tissues = [
        'Adipose_Subcutaneous',
        'Adipose_Visceral_Omentum',
        'Adrenal_Gland',
        'Brain_Anterior_cingulate_cortex_BA24',
        'Brain_Caudate_basal_ganglia',
        'Brain_Cerebellar_Hemisphere',
        'Brain_Cerebellum',
        'Brain_Cortex',
        'Brain_Frontal_Cortex_BA9',
        'Brain_Hippocampus',
        'Brain_Hypothalamus',
        'Brain_Nucleus_accumbens_basal_ganglia',
        'Brain_Putamen_basal_ganglia',
        'Breast_Mammary_Tissue',
        'Cells_EBV-transformed_lymphocytes', # not sure if it's better to include or not
        'Colon_Sigmoid',
        'Colon_Transverse',
        'Heart_Atrial_Appendage',
        'Heart_Left_Ventricle',
        'Liver',
        'Lung',
        'Muscle_Skeletal',
        'Ovary',
        'Prostate',
        'Testis',
        'Thyroid',

    ]

    #for variant, tissue, slope, spidex_record, g in iterate_gtex_vs_spidex(tissues=spidex_tissues, strict=True, only_the_same_gene=True):
    for variant, tissue, slope, spidex_record, g in iterate_gtex_vs_spidex(tissues=spidex_tissues, strict=True):
        slope = float(slope)
        z_score = float(spidex_record.dpsi_zscore)

        if variant not in effect_sizes:
            effect_sizes[variant] = slope
            z_scores[variant] = z_score

        elif abs(slope) > abs(effect_sizes[variant]):
            effect_sizes[variant] = slope
            assert z_score == z_scores[variant]

    print('Mutations mapped:', len(effect_sizes))
    print('Spidex records:')
    print('GTEx records from selected tissues:', count_all(spidex_tissues))
    coef, p_value = pearsonr(list(effect_sizes.values()), list(z_scores.values()))
    print('Pearson\'s r', coef)
    print('P-value', p_value)
    coef, p_value = spearmanr(list(effect_sizes.values()), list(z_scores.values()))
    print('Spearman\'s r', coef)
    print('P-value', p_value)


def iterate_gtex_vs_spidex(strict=False, tissues=GTEX_TISSUES, location=None, filters=None, only_the_same_gene=False):
    """
    Yield records representing GTEx-SPIDEX pairs which are matching
    (the same position, reference and alternative alleles).

    Args:
        strict: should critical data discrepancies raise errors or be collected as statistics?
        tissues: list of tissues to be used
        only_the_same_gene: accept only mutations acting on their gene (strong cis- acting)
        location: 'intronic' or 'exonic' - filters SPIDEX records
        filters: dict - filtering criteria for GTEx records.
    
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

    # Use "Brain cortex" for basic tests - it's very small
    # tissues = ['Brain_Cortex']
    # Use "Adipose Subcutaneous" for larger tests
    # tissues = ['Adipose_Subcutaneous']

    path = create_path_for_genes_db(tissues)
    genes = ExpressedGenes(path)
    genes.reset()

    import_expressed_genes(
        genes,
        tissues=tissues
    )

    tb = tabix.open(SPIDEX_LOCATION)

    count = count_all(tissues)

    counter = Counter()

    for mutation_code, tissue, slope, ensembl_gene_id in tqdm(iterate_over_expression(tissues), total=count):
        chrom, pos, ref, alt, _ = mutation_code.split('_')

        # In spidex there are only SNPs (single!)
        if len(ref) != 1 or len(alt) != 1:
            counter['not_single'] += 1
            continue

        pos = int(pos)

        gene = Gene(*genes[ensembl_gene_id])

        if not gene:
            print('gene %s not present in data' % ensembl_gene_id)
            continue

        if only_the_same_gene:
            if not (gene.start <= pos <= gene.end):
                counter['not_within_gene'] += 1
                continue

        variant = SingleAltVariant(
            chr_name=chrom,
            chr_start=pos,
            chr_end=pos,
            chr_strand=gene.strand,
            snp_id='-',
            ref=convert_to_strand(ref, gene.strand),
            alt=convert_to_strand(alt, gene.strand),
            gene=ensembl_gene_id
        )

        records = spidex_get_variant(tb, variant)
        if filters:
            records = [
                record
                for record in records
                if all(
                    getattr(record, key) == value
                    for key, value in filters.items()
                )
            ]

        record = None

        try:
            # if genes are the same there is no need to test strands, but its better to double check
            record = choose_record(records, variant, variant.alt, strict=strict, test_strand=True, location=location)
        except StrandMismatch:
            counter['strand_mismatch'] += 1
        except Intronic:
            counter['intronic'] += 1

        if record:
            if gene.name != record.gene:
                counter['gene_name_mismatch'] += 1
                continue

            variant.refseq_transcript = record.transcript

            yield variant, tissue, slope, record, gene
        else:
            counter['Not found in SPIDEX'] += 1

    if strict:
        print(counter)
