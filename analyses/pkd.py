import requests
import sys
from collections import defaultdict
from analyses import variants_getter

from biomart_data import BiomartDataset
from cache import cacheable
from parse_variants import decode_hgvs_code
from variant import Variant
from recordclass import recordclass


PKDB_URL = 'http://pkdb.mayo.edu/cgi-bin/v2_display_mutations.cgi?apkd_mode=PROD&username='

Gene = recordclass(
    'Gene',
    [
        'ensembl_transcript_stable_id',
        'chrom', 'strand',
        'tx_start', 'tx_end',
        'cds_start', 'cds_end',
    ]
)


"""
In PDKB they used following genes:
PKD1: NM_00296.2; PKD2: NM_000297

To quickly map mutations to genome in GRCh 37 I used http://genome.ucsc.edu/cgi-bin/hgTables
#bin	name	chrom	strand	txStart	txEnd	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	score	name2	cdsStartStat	cdsEndStat	exonFrames

PKD1:
601	NM_000296	chr16	-	2138710	2185899	2139727	2185690	46	2138710,2140285,2140674,2140884,2141423,2141781,2142047,2142480,2142954,2143544,2143811,2144092,2147148,2147319,2147728,2147868,2149644,2149861,2150166,2150396,2152061,2152381,2152814,2153266,2154498,2155322,2155865,2156091,2156398,2156805,2157883,2158252,2162340,2162788,2163161,2164170,2165378,2165992,2166529,2166833,2167489,2167791,2168676,2169114,2169307,2185475,	2140195,2140591,2140809,2141175,2141598,2141907,2142189,2142593,2143094,2143739,2144014,2144211,2147242,2147504,2147778,2147985,2149771,2150072,2150310,2150567,2152257,2152634,2152971,2153896,2154643,2155475,2156025,2156305,2156678,2156949,2158033,2161872,2162474,2162964,2163293,2164926,2165626,2166119,2166645,2167054,2167673,2168463,2168846,2169186,2169379,2185899,	0	PKD1	cmpl	cmpl	0,0,0,0,2,2,1,2,0,0,1,2,1,2,0,0,2,1,1,1,0,2,1,1,0,0,2,1,0,0,0,1,2,0,0,0,1,0,1,2,1,1,2,2,2,0,

PKD2:
157	NM_000297	chr4	+	88928798	88998931	88928885	88996846	15	88928798,88940609,88957371,88959402,88964384,88967793,88973142,88977237,88979134,88983057,88986525,88986913,88989049,88995963,88996609,	88929480,88940723,88957505,88959653,88964609,88968022,88973310,88977419,88979255,88983156,88986647,88987031,88989213,88996111,88998931,	0	PKD2	cmpl	cmpl	0,1,1,0,2,2,0,0,2,0,0,2,0,2,0,

(i.e., c.1232_1234delATC in lieu of g.44003delATC or p.LEU440del).

I converted refseq ids to ensembl using biomart: http://grch37.ensembl.org/biomart/

<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
    <Dataset name = "hsapiens_gene_ensembl" interface = "default" >
        <Filter name = "refseq_mrna" value = "NM_000296"/>
        <Attribute name = "ensembl_gene_id" />
        <Attribute name = "ensembl_transcript_id" />
    </Dataset>
</Query>

"""

GENES = {
    'PDK1': Gene(
        ensembl_transcript_stable_id='ENST00000423118',
        chrom='chr16', strand='-',
        tx_start=2138710, tx_end=2185899,
        cds_start=2139727, cds_end=2185690
    )
}


@cacheable
def get_html_page(url):

    r = requests.post(
        url,
        data={
            'GENE': 'PKD1',
            'GERM': 'germline',  # or somatic
            'CLINICAL': 'All',
            'apkd_mode': 'PROD',
            'username': ''
        }
    )

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    return r.text


def get_raw_table():
    from bs4 import BeautifulSoup

    html = get_html_page.load_or_create(PKDB_URL)
    soup = BeautifulSoup(html, 'html.parser')

    variants_table = soup.find('table', {'id': 'main_table'})
    rows = variants_table.find_all('tr')

    raw_rows = []

    for row in rows:
        cells = [
            cell.get_text().strip()
            for cell in row.find_all('td')
            ]
        raw_rows.append(cells)

    return raw_rows


@variants_getter
def polycystic_kidney_disease_variants(exonic_only=True):
    """Fetch variants associated with PKD1 from PKDB."""
    variants_by_genes = defaultdict(list)

    table = get_raw_table()
    rows = iter(table)

    headers = next(rows)

    headers = [
        header.lower().replace(' ', '_').replace('#', 'count').replace('%', 'percent')
        for header in headers
    ]

    VariantRecord = recordclass('VariantRecord', headers)

    for row in rows:
        record = VariantRecord(*row)

        if exonic_only:
            if not record.region.startswith('EX'):
                continue

        # currently hardcoded
        hgvs_code = 'PDK1:c.' + record.cdnachange

        try:
            gene_name, pos, ref, alt = decode_hgvs_code(hgvs_code)
        except ValueError as orginal_e:
            # handle manually mutation code which do not comply with HGVS

            # range specified with '-' instead of underscore:
            if '-' in hgvs_code:
                print('Assuming - means range in %s.' % hgvs_code)
                hgvs_code = hgvs_code.replace('-', '_')

                try:
                    gene_name, pos, ref, alt = decode_hgvs_code(hgvs_code)
                except ValueError as e:
                    print(e.message)
                    continue
            else:
                print(orginal_e.message)
                continue

        gene = GENES[gene_name]

        variant = Variant(
            None, ref=ref, alts=[alt],
            chr_name=gene.chrom,
            chrom_start=gene.cds_start + pos,
            ensembl_transcript_stable_id=gene.ensembl_transcript_stable_id,
            refsnp_id=hgvs_code
        )
        # TODO: tried to use
        # http://myvariant.info/ to map variants, get variants data.
        # as_hgvs won work until ref is given.
        # back in the same point
        print(variant.as_hgvs())

        variants_by_genes[gene_name].append(variant)

    return variants_by_genes
