#!/usr/bin/env python3

import sys
import argparse
import json
import hail as hl
from utils.hail_utils import join_tables, add_percentiles
from analyses_for_manuscript.scores_to_use import RAW_SCORES
from ukbb_common.utils.generic import create_broadcast_dict



vsm_f = sys.argv[1]

with open(vsm_f, 'r') as f: 
    data = json.load(f)

score_list = []
for vsm in data['vsm_info']: 
    score_list.append(vsm['score_col'])

print(1)

final_linker = 'gs://missense-scoring/mutation/linker_with_uniprot_name.ht'
dir_to_save = data['dir_to_save']
join_all_ht = dir_to_save + 'join_all.ht'
percentile_all_ht = dir_to_save +  'percentile_all.ht'
median_by_gene_ht = dir_to_save + 'median_by_gene.ht'

print(2)

def join_all():
    ht = hl.read_table(final_linker)
    for vsm in data['vsm_details']:
      #  print(score_table)
        vsm_ht = hl.read_table(vsm['hail_path'])
        score_col = vsm['score_col']
        ht = join_tables(ht, vsm_ht, score_col, filter=False)
        print(3)

        if vsm['higher_is_less_deleterious'] == 'True': 
            ht = ht.annotate(
                **{f'{score_col}_neg': 1-ht[score_col]}
            )
        print(4)

    print('writing file')
    ht.write('gs://genetics-gym-not-public/Trisha/join_all.ht', True)
    print('done')

def percentile_all():
    ht = hl.read_table('gs://missense-scoring/all_models_scores.ht')
    ht = ht.filter(
        hl.all([hl.is_defined(ht[x]) for x in score_list])
    )
    ht = add_percentiles(ht, score_list)

    for vsm in data['vsm_details']:
        if vsm['higher_is_less_deleterious'] == 'True':
            score_col = vsm['score_col']
            ht = ht.annotate(
                **{f'{score_col}_neg_percentile': 1-ht[f'{score_col}_percentile']}
            )

    print('writing file')
    ht.write('gs://genetics-gym-not-public/Trisha/percentiles.ht', True)
    print('done')

def median_by_gene():
    ht = hl.read_table('gs://missense-scoring/all_models_scores_intersection_percentiles.ht')
    ht.group_by(ht.uniprot_id).aggregate(
        **{f'{s}_gene_median': hl.agg.approx_median(ht[s], k=1000)
           for s in RAW_SCORES},
        **{f'{s}_gene_mean': hl.agg.mean(ht[s])
           for s in RAW_SCORES}
    ).write('gs://genetics-gym-not-public/Trisha/temp_agg.ht', True)

    agg_ht = hl.read_table('gs://missense-scoring/temp_agg.ht')
    ht.annotate(
        **create_broadcast_dict(agg_ht.uniprot_id)[ht.uniprot_id]
    ).write('gs://genetics-gym-not-public/Trisha/temp2_agg.ht', True)

    ht = hl.read_table('gs://missense-scoring/temp2_agg.ht')
    median_score_list = [f'{s}_gene_median' for s in RAW_SCORES]
    mean_score_list = [f'{s}_gene_mean' for s in RAW_SCORES]
    ht = add_percentiles(ht, median_score_list)
    ht = add_percentiles(ht, mean_score_list)

    ht.write('gs://genetics-gym-not-public/Trisha/all_models_scores_intersection_with_aggs_percentile.ht', True)

if __name__ == '__main__':
    print('hi')
    join_all()
    # percentile_all()
    #median_by_gene()

