# Mouse, alpha, antigen.species figure

# Import packages
import pandas as pd
import numpy as np
from tcrdist.repertoire import TCRrep
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE

# Read dataset
tcell_df = pd.read_csv("vdjdb.txt", sep='\t')

# Clean the data, remove any v.segm or j.segm is NaN
t_cell_df = tcell_df[tcell_df['v.segm'].notnull() & tcell_df['j.segm'].notnull()]

t_cell_df.reset_index(drop=True, inplace=True)
# Initialize dataframe, just select some columns
clean_t_cell_df = t_cell_df[
    ['complex.id', 'gene', 'cdr3', 'v.segm', 'j.segm', 'species', 'antigen.species', 'antigen.epitope']]

# --------------------------------------------------------------------------------------------------------------
# Select gene is TRA and species is MusMusculus(mouse)
mouse_alpha_df = clean_t_cell_df[(clean_t_cell_df['gene'] == 'TRA') & (clean_t_cell_df['species'] == 'MusMusculus')]
# Rename the columns
mouse_alpha_df.rename(columns={'cdr3': 'cdr3_a_aa', 'v.segm': 'v_a_gene', 'j.segm': 'j_a_gene'}, inplace=True)
# Only needs cdr3_b_aa,v_b_gene and j_b_gene columns
mouse_alpha_df_as = mouse_alpha_df[['cdr3_a_aa', 'v_a_gene', 'j_a_gene', 'antigen.species']]
mouse_alpha_df_ae = mouse_alpha_df[['cdr3_a_aa', 'v_a_gene', 'j_a_gene', 'antigen.epitope']]
# Reset index
mouse_alpha_df_as.reset_index(drop=True, inplace=True)
mouse_alpha_df_ae.reset_index(drop=True, inplace=True)

# Select gene is TRB and species is MusMusculus(mouse)
mouse_beta_df = clean_t_cell_df[(clean_t_cell_df['gene'] == 'TRB') & (clean_t_cell_df['species'] == 'MusMusculus')]
# Rename the columns
mouse_beta_df.rename(columns={'cdr3': 'cdr3_b_aa', 'v.segm': 'v_b_gene', 'j.segm': 'j_b_gene'}, inplace=True)
# Only needs cdr3_b_aa,v_b_gene and j_b_gene columns
mouse_beta_df_as = mouse_beta_df[['cdr3_b_aa', 'v_b_gene', 'j_b_gene', 'antigen.species']]
mouse_beta_df_ae = mouse_beta_df[['cdr3_b_aa', 'v_b_gene', 'j_b_gene', 'antigen.epitope']]
# Reset index
mouse_beta_df_as.reset_index(drop=True, inplace=True)
mouse_beta_df_ae.reset_index(drop=True, inplace=True)


# --------------------------------------------------------------------------------------------------------------

# This function is for mouse only
def calculate_tr_alpha(input_df):
    tr_alpha_mouse = TCRrep(cell_df=input_df,
                            organism='mouse',
                            chains=['alpha'],
                            db_file='alphabeta_gammadelta_db.tsv')
    mouse_alpha_matrix = tr_alpha_mouse.pw_alpha

    df_alpha = tr_alpha_mouse.clone_df
    # return the distance matrix and corresponding dataframe
    return mouse_alpha_matrix, df_alpha


# print(mouse_alpha_matrix)
# print(mouse_alpha_df)
# print(df_alpha)
# print(tr_alpha_mouse.clone_df.columns)


def calculate_tr_beta(input_df):
    tr_beta_mouse = TCRrep(cell_df=input_df,
                           organism='mouse',
                           chains=['beta'],
                           db_file='alphabeta_gammadelta_db.tsv')
    mouse_beta_matrix = tr_beta_mouse.pw_beta
    df_beta = tr_beta_mouse.clone_df
    return mouse_beta_matrix, df_beta


# Define plot function, focus on TCR and antigen species relation
def plot2d_antigen_species_mouse(chain_name, input_matrix, clone_df, random_state):
    # Here specificity is assumed as antigen.species
    specificity = clone_df['antigen.species']
    # Use t-SNE to reduce dimension of distance matrix
    tsne_reducer = TSNE(n_components=2, random_state=random_state, learning_rate='auto', init='pca')
    tsne_results = tsne_reducer.fit_transform(input_matrix)

    # Get unique antigen species
    unique_antigen_species = np.unique(specificity)
    # Generate colors [R,G,B, transparency]
    colors = plt.cm.jet(np.linspace(0, 1, len(unique_antigen_species)))
    # Create a dictionary, contains each antigen species and its color
    # Mapping color
    color_map = dict(zip(unique_antigen_species, colors))

    # -------------------------------------------------------------------------
    # Visualization
    plt.figure(figsize=(12, 20))

    # For loop obtains corresponding antigen.species and color
    for antigen_species, color in color_map.items():
        # specificity is a series contains all antigen species
        true_indiex = (specificity == antigen_species)
        # Plot scatter
        plt.scatter(tsne_results[true_indiex, 0], tsne_results[true_indiex, 1], color=color, label=antigen_species,
                    alpha=0.7)

    # Add title, xy label and legend
    plt.title('t-SNE 2-dimensional plot of TCRs based on antigen.species ({})'.format(chain_name), fontweight='bold', fontsize=15)
    plt.xlabel('t-SNE x-axis')
    plt.ylabel('t-SNE y-axis')
    plt.legend(title='antigen.species', bbox_to_anchor=(1, 1), loc='upper left', fontsize=8)
    plt.show()


# Define plot functon focus on TCR and the antigen epitope
def plot2d_antigen_epitope_mouse(chain_name, input_matrix, clone_df, random_state):
    # Here specificity is assumed as antigen.epitope
    specificity = clone_df['antigen.epitope']
    # Use t-SNE to reduce dimension of distance matrix
    tsne_reducer = TSNE(n_components=2, random_state=random_state, learning_rate='auto', init='pca')
    tsne_results = tsne_reducer.fit_transform(input_matrix)

    # Get unique antigen epitope
    unique_antigen_epitope = np.unique(specificity)
    # Generate colors [R,G,B, transparency]
    colors = plt.cm.jet(np.linspace(0, 1, len(unique_antigen_epitope)))
    # Create a dictionary, contains each antigen epitope and its color
    # Mapping color
    color_map = dict(zip(unique_antigen_epitope, colors))

    # -------------------------------------------------------------------------
    # Visualization
    plt.figure(figsize=(12, 20))

    # For loop obtains corresponding antigen.epitope and color
    for antigen_epitope, color in color_map.items():
        # specificity is a series contains all antigen epitope
        true_indiex = (specificity == antigen_epitope)
        # Plot scatter
        plt.scatter(tsne_results[true_indiex, 0], tsne_results[true_indiex, 1], color=color, label=antigen_epitope,
                    alpha=0.7)

    # Add title, xy label and legend
    plt.title('t-SNE 2-dimensional plot of TCRs based on antigen.epitope ({})'.format(chain_name), fontweight='bold', fontsize=15)
    plt.xlabel('t-SNE x-axis')
    plt.ylabel('t-SNE y-axis')
    plt.legend(title='antigen.epitope', bbox_to_anchor=(1, 1), loc='upper left', fontsize=6)
    plt.show()


# Mouse alpha chain(for antigen species and epitope)
mouse_alpha_matrix_as, df_alpha_as = calculate_tr_alpha(mouse_alpha_df_as)
mouse_alpha_matrix_ae, df_alpha_ae = calculate_tr_alpha(mouse_alpha_df_ae)
# Mouse beta chain(for antigen species and epitope)
mouse_beta_matrix_as, df_beta_as = calculate_tr_beta(mouse_beta_df_as)
mouse_beta_matrix_ae, df_beta_ae = calculate_tr_beta(mouse_beta_df_ae)

# Plot 2-dimensional figures
plot2d_antigen_species_mouse('alpha', mouse_alpha_matrix_as, df_alpha_as, 46)
plot2d_antigen_species_mouse('beta', mouse_beta_matrix_as, df_beta_as, 45)
plot2d_antigen_epitope_mouse('alpha', mouse_alpha_matrix_ae, df_alpha_ae, 47)
plot2d_antigen_epitope_mouse('beta', mouse_beta_matrix_ae, df_beta_ae, 48)
