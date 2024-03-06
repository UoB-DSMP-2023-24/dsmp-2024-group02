# Import packages
import pandas as pd
import numpy as np
from tcrdist.repertoire import TCRrep

# # read the vdjdb.txt as dataframe
# tcell_df = pd.read_csv("C:/Users/86137/Desktop/mini_project/vdjdb-2023-06-01/vdjdb.txt", sep='\t')
# # print(tcell_df)
#
# url = 'https://raw.githubusercontent.com/kmayerb/tcrdist2/API2/tcrdist/test_files_compact/dash.csv'
# data = pd.read_csv(url)
# # print(data.head())
#
# """
# If you just want a 'tcrdistances' using pre-set default setting.
#
#     You can access distance matrices:
#         tr.pw_alpha     - alpha chain pairwise distance matrix
#         tr.pw_beta      - alpha chain pairwise distance matrix
#         tr.pw_cdr3_a_aa - cdr3 alpha chain distance matrix
#         tr.pw_cdr3_b_aa - cdr3 beta chain distance matrix
# """
# tr = TCRrep(cell_df = data,
#             organism = 'mouse',
#             chains = ['alpha','beta'],
#             db_file = 'alphabeta_gammadelta_db.tsv')
#
# # print(tr.pw_alpha)
# #tr.pw_beta
# #tr.pw_cdr3_a_aa
# #tr.pw_cdr3_b_aa


def main():
    # Read the vdjdb.txt from local
    tcell_df = pd.read_csv("C:/Users/86137/Desktop/mini_project/vdjdb-2023-06-01/vdjdb.txt", sep='\t')
    # Initialize dataframe, just select some columns
    t_cell_df = tcell_df[['complex.id', 'gene', 'cdr3', 'v.segm', 'j.segm', 'species']]
    # clean the data, remove any v.segm or j.segm is NaN
    clean_t_cell_df = t_cell_df[t_cell_df['v.segm'].notnull() & t_cell_df['j.segm'].notnull()]

    # --------------------------------------------------------------------------------------------------------------
    # select gene is TRA and species is HomoSapiens(human)
    human_alpha_df = clean_t_cell_df[(clean_t_cell_df['gene'] == 'TRA') & (clean_t_cell_df['species'] == 'HomoSapiens')]
    # Rename the columns
    human_alpha_df.rename(columns={'cdr3': 'cdr3_a_aa', 'v.segm': 'v_a_gene', 'j.segm': 'j_a_gene'}, inplace=True)
    # Only needs cdr3_a_aa,v_a_gene and j_a_gene columns
    human_alpha_df = human_alpha_df[['cdr3_a_aa', 'v_a_gene', 'j_a_gene']]
    # Reset index
    human_alpha_df.reset_index(drop=True)
    # print(human_alpha_df)

    # Human alpha chain
    tr_alpha_human = TCRrep(cell_df=human_alpha_df,
                            organism='human',
                            chains=['alpha'],
                            compute_distances=False)

    tr_alpha_human.cpus = 2
    # Original matrix is too big, sparse matrix is used here
    tr_alpha_human.compute_sparse_rect_distances(radius=50, chunk_size=100)

    # --------------------------------------------------------------------------------------------------------------
    # Select gene is TRB and species is HomoSapiens(human)
    human_beta_df = clean_t_cell_df[(clean_t_cell_df['gene'] == 'TRB') & (clean_t_cell_df['species'] == 'HomoSapiens')]
    # Rename the columns
    human_beta_df.rename(columns={'cdr3': 'cdr3_b_aa', 'v.segm': 'v_b_gene', 'j.segm': 'j_b_gene'}, inplace=True)
    # Only needs cdr3_b_aa,v_b_gene and j_b_gene columns
    human_beta_df = human_beta_df[['cdr3_b_aa', 'v_b_gene', 'j_b_gene']]
    # Reset index
    human_beta_df.reset_index(drop=True)
    # print(human_beta_df)

    # Human beta chain
    tr_beta_human = TCRrep(cell_df=human_beta_df,
                           organism='human',
                           chains=['beta'],
                           compute_distances=False)

    tr_beta_human.cpus = 2
    tr_beta_human.compute_sparse_rect_distances(radius=50, chunk_size=100)

    # --------------------------------------------------------------------------------------------------------------
    # select gene is TRA and species is MusMusculus(mouse)
    MusMusculus_alpha_df = clean_t_cell_df[(clean_t_cell_df['gene'] == 'TRA') & (clean_t_cell_df['species'] == 'MusMusculus')]
    # Rename the columns
    MusMusculus_alpha_df.rename(columns={'cdr3': 'cdr3_a_aa', 'v.segm': 'v_a_gene', 'j.segm': 'j_a_gene'}, inplace=True)
    # Only needs cdr3_a_aa,v_a_gene and j_a_gene columns
    MusMusculus_alpha_df = MusMusculus_alpha_df[['cdr3_a_aa', 'v_a_gene', 'j_a_gene']]
    # Reset index
    MusMusculus_alpha_df.reset_index(drop=True)

    # print(MusMusculus_alpha_df)

    # Mouse alpha chain
    tr_alpha_mouse = TCRrep(cell_df=MusMusculus_alpha_df,
                            organism='mouse',
                            chains=['alpha'],
                            compute_distances=False)

    tr_alpha_mouse.cpus = 2
    # Original matrix is too big, sparse matrix is used here
    tr_alpha_mouse.compute_sparse_rect_distances(radius=50, chunk_size=100)

    # --------------------------------------------------------------------------------------------------------------
    # select gene is TRB and species is MusMusculus(mouse)
    MusMusculus_beta_df = clean_t_cell_df[(clean_t_cell_df['gene'] == 'TRA') & (clean_t_cell_df['species'] == 'MusMusculus')]
    # Rename the columns
    MusMusculus_beta_df.rename(columns={'cdr3': 'cdr3_a_aa', 'v.segm': 'v_a_gene', 'j.segm': 'j_a_gene'}, inplace=True)
    # Only needs cdr3_a_aa,v_a_gene and j_a_gene columns
    MusMusculus_beta_df = MusMusculus_beta_df[['cdr3_a_aa', 'v_a_gene', 'j_a_gene']]
    # Reset index
    MusMusculus_beta_df.reset_index(drop=True)

    # print(MusMusculus_alpha_df)

    # Mouse beta chain
    tr_beta_mouse = TCRrep(cell_df=MusMusculus_beta_df,
                            organism='mouse',
                            chains=['beta'],
                            compute_distances=False)

    tr_beta_mouse.cpus = 2
    # Original matrix is too big, sparse matrix format is used here
    tr_beta_mouse.compute_sparse_rect_distances(radius=50, chunk_size=100)



    print(tr_beta_human.rw_beta)

if __name__ == '__main__':
    # This is the "entry point" guard for multiprocessing on Windows.
    main()

