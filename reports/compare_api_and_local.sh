cut -f1 expression_table_for_variants_with_tissues_\(based_on_data_from_gtex\).tsv | sort | uniq > gtex_local_ids
cut -f1 API_expression_table_for_variants_with_tissues_\(based_on_data_from_gtex\).tsv | sort | uniq > gtex_api_ids

diff gtex_api_ids gtex_local_ids
