import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
import seaborn as sns
from scipy.spatial.distance import pdist, squareform
from scipy.stats import chi2
from dtaidistance import dtw
from dtaidistance import dtw_visualisation as dtwvis


def squared_euclidean_distance(series1, series2):
    """Calculates the squared Euclidean distance between two series."""
    return ((series1 - series2) ** 2).sum()

def get_intra_distance_for_WT(gene_expression_data: pd.DataFrame) -> pd.DataFrame:
    """Calculates pairwise distances within groups defined by 'light_regime' and 'plate' for WT mutant.

    Args:
        gene_expression_data: DataFrame containing gene expression data, with columns
            prefixed by 'y2_' and 'ynpq_' and a 'mutant_ID' column.

    Returns:
        DataFrame with pairwise distances, 'light_regime', 'plate', mean intra-group distances,
        and total WT cells per plate and light regime.
    """
    # Filter data for mutant_ID == 'WT'
    wt_data = gene_expression_data[gene_expression_data['mutant_ID'] == 'WT']

    # Count total WT cells per plate and light regime
    wt_counts = wt_data.groupby(['plate', 'light_regime']).size().reset_index(name='total_wt_cells')

    # Group by 'light_regime' and 'plate'
    grouped = wt_data.groupby(['light_regime', 'plate'])

    # Initialize lists to store distances and group information
    all_distances_y2 = []
    all_distances_ynpq = []
    all_light_regimes = []
    all_plates = []
    all_well_ids = []

    # Iterate over each group and calculate pairwise distances
    for (light_regime, plate), group in grouped:
        group = group.sort_values('well_id')  # ensure stable order
        # Extract and validate data for y2 and ynpq
        y2_data = group.filter(like='y2_').dropna(axis=1).values
        ynpq_data = group.filter(like='ynpq_').dropna(axis=1).values
        well_ids = group['well_id'].values

        if y2_data.shape[0] > 1 and y2_data.shape[1] > 0:
            distances_y2 = pdist(y2_data, metric=squared_euclidean_distance)
            all_distances_y2.extend(distances_y2)
            # Generate well ID pairs for y2
            well_id_pairs_y2 = [(well_ids[i], well_ids[j]) for i in range(len(well_ids)) for j in range(i + 1, len(well_ids))]
            all_well_ids.extend(well_id_pairs_y2)
            
        if ynpq_data.shape[0] > 1 and ynpq_data.shape[1] > 0:
            distances_ynpq = pdist(ynpq_data, metric=squared_euclidean_distance)
            all_distances_ynpq.extend(distances_ynpq)

        # Store group information for each pairwise comparison
        n_comparisons_y2 = len(distances_y2) if y2_data.shape[0] > 1 and y2_data.shape[1] > 0 else 0
        n_comparisons_ynpq = len(distances_ynpq) if ynpq_data.shape[0] > 1 and ynpq_data.shape[1] > 0 else 0
        n_comparisons = max(n_comparisons_y2, n_comparisons_ynpq)

        all_light_regimes.extend([light_regime] * n_comparisons)
        all_plates.extend([plate] * n_comparisons)

    # Create the result DataFrame
    result_df = pd.DataFrame({
        'light_regime': all_light_regimes,
        'plate': all_plates,
        'well_id_pair': all_well_ids,  # Add well ID pairs
        'pairwise_distances_y2': all_distances_y2,
        'pairwise_distances_ynpq': all_distances_ynpq,
    })

    # Calculate mean intra-group distances for each plate and light regime
    mean_distances = result_df.groupby(['plate', 'light_regime'])['pairwise_distances_y2'].transform('mean')
    result_df['mean_intra_distance_y2'] = mean_distances

    # Merge with the total WT cell counts
    result_df = result_df.merge(wt_counts, on=['plate', 'light_regime'], how='left')

    return result_df