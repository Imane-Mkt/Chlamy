import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np



def plot_raw_y(data, ts = 'y2', which_plate = 'all', mutant = 'all'):
    """
    This function plots overlapping time series for each light regime for a specified plate or WT or all
    Valid values:
    - ts: ['y2', 'ynpq']  type of time series to plot ('y2', 'ynpq', etc.)
    - which_plate: 'all' or a valid plate prefix (e.g., '33', '98').
    - mutant: ['all', 'WT', 'nonWT']
    """
    # Identifier les colonnes y2 et elapsed_hours
    y2_columns = [col for col in data.columns if col.startswith(ts) and col.split('_')[-1].isdigit() and 'std' not in col]
    elapsed_hours_columns = [col for col in data.columns if col.startswith('elapsed_hours_')]

    # Régimes lumineux
    light_regimes = ['20h_HL', '20h_ML', '2h-2h', '10min-10min', '5min-5min', '1min-5min', '1min-1min', '30s-30s']

    # Créer une grille 4x4 pour les graphiques
    fig, axes = plt.subplots(4, 4, figsize=(16, 12))  # 4x4 grille
    axes = axes.flatten()  # Aplatir la grille pour un accès facile

    # Tracer les graphiques pour chaque régime lumineux
    for i, light in enumerate(light_regimes):
        ax = axes[i]  # Sélectionner l'axe correspondant
        data_light = data[data['light_regime'] == light]
        if which_plate != 'all':
            if which_plate == '99':
                data_light[data_light['plate']=='99']
            data_light[data_light['plate'].astype(str).str.startswith(which_plate)]  # Filtrer les données pour le régime lumineux actuel
        
        labeled_plates = set()

        # Tracer les courbes pour chaque plaque
        for plate in data_light['plate'].unique():
            data_plate = data_light[(data_light['plate'] == plate)]
            if mutant == 'WT':
                data_plate[data_plate['mutant_ID'].isin(['WT', 'WT CLiP'])]
            if mutant == 'nonWT':
                data_plate[data_plate['mutant_ID'].isnotin(['WT', 'WT CLiP'])]
                
            for _, row in data_plate.iterrows():
                y2_values = row[y2_columns].values  # Extraire les valeurs y2
                time_values = row[elapsed_hours_columns].astype(float).values  # Extraire les heures écoulées
                # Filtrer les valeurs valides (non-NaN)
                valid_indices = ~np.isnan(time_values)
                time_values = time_values[valid_indices]

                y2_values = y2_values[:len(time_values)]  # Assurez-vous que y2_values a la même longueur que time_values
                time_values = time_values[:len(y2_values)]  # Assurez-vous que time_values a la même longueur que y2_values
                # Ajouter un label uniquement si la plaque n'a pas encore été étiquetée
                if plate not in labeled_plates:
                    ax.plot(time_values, y2_values, alpha=0.7, label=f'Plate {plate}')
                    labeled_plates.add(plate)  # Marquer la plaque comme étiquetée
                else:
                    ax.plot(time_values, y2_values, alpha=0.7)
        
        # Ajouter des détails au graphique
        ax.set_title(f'Régime lumineux : {light}')
        ax.set_xlabel('Elapsed Hours')
        ax.set_ylabel(f'{ts} for {mutant} in {which_plate}')
        ax.grid(True)

    # Supprimer les axes inutilisés si moins de 16 régimes lumineux
    for j in range(len(light_regimes), len(axes)):
        fig.delaxes(axes[j])

    # Ajuster l'espacement entre les graphiques
    plt.tight_layout()
    plt.show()
    

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load your dataset (wide format: one slope column per regime)
def plot_slope_correlation(data, regime_order = ["20h_HL", "2h-2h","10min-10min", "5min-5min", "1min-5min", "1min-1min", "30s-30s"], phase=2):
    labels = {
    '20h_HL': 'HL',
    '2h-2h': '2h',
    '10min-10min': '10min',
    '5min-5min': '5min',
    '1min-1min': '1min',
    '30s-30s': '30s'
}

    # Pivot to wide format: one row per mutant_ID, one column per light regime
    pivoted = data.pivot_table(
        index="mutant_ID",
        columns="light_regime",
        values="y2_slope",
        aggfunc="mean"
    )
    # Define desired light regime order

    # Reorder columns after pivot
    pivoted = pivoted[regime_order]

    # Drop rows with too many missing values (optional)
    pivoted = pivoted.dropna(thresh=3)  # keep rows with at least 3 regimes

    # Compute correlation
    corr_matrix = pivoted.corr()

    # Plot heatmap
    #plt.figure(figsize=(10, 8))
    # sns.heatmap(corr_matrix, cmap="magma", annot=True, fmt=".2f",
    #             vmin=0, vmax=1, linewidths=0.5)
    # plt.title("Correlation of average y2_slope between light regimes")
    # plt.xticks(rotation=45, ha="right")
    # plt.tight_layout()
    # plt.show()
    ax = sns.heatmap(corr_matrix, cmap="magma", annot=True, fmt=".2f",
                 vmin=0, vmax=1, linewidths=0.5)

    plt.title(f"Correlation between slopes in Phase{phase}")
    ax.set_xlabel("")
    ax.set_ylabel("")

    # Change xtick labels
    ax.set_xticklabels([labels.get(t.get_text(), t.get_text()) for t in ax.get_xticklabels()], rotation=45, ha="right")

    # Change ytick labels
    ax.set_yticklabels([labels.get(t.get_text(), t.get_text()) for t in ax.get_yticklabels()], rotation=0)

    plt.tight_layout()
    plt.show()


def plot_slope_lights_trend_log(data, genes):
    import numpy as np
    import pandas as pd
    from numpy.polynomial.polynomial import Polynomial
    import matplotlib.cm as cm
    import matplotlib.colors as mcolors

    # Load data

    # Light regime order and frequency mapping
    # regime_order = ['20h_HL', '20h_ML', '2h-2h', '10min-10min', '5min-5min', '1min-1min', '30s-30s']
    # regime_to_freq = {regime: i / (len(regime_order) - 1) for i, regime in enumerate(regime_order)}

    #  frequencies in 1/min
    regime_to_freq = {
        '20h_HL': 10000,                     
        '2h-2h': 120,                
        '10min-10min': 10,          
        '5min-5min': 5,
        '1min-1min': 1,
        '30s-30s': 0.5           
    }
    regime_order = list(regime_to_freq.keys())

    for gene in genes:
        # if gene not in data['mutated_genes'].unique():
        #     print(f'{gene} not present in dataset')
        #     continue
        gene_data = data[
            (data['mutated_genes'].astype(str).str.contains(gene)) &
            (data['light_regime'].isin(regime_order)) &
            (~data['y2_slope'].isna())
        ].copy()
        gene_data['freq'] = gene_data['light_regime'].map(regime_to_freq)

        plates = gene_data['plate'].unique()
        wt_data = data[
            (data['mutated_genes'] == 'WT') &
            (data['plate'].isin(plates)) &
            (data['light_regime'].isin(regime_order)) &
            (~data['y2_slope'].isna())
        ].copy()
        wt_data['freq'] = wt_data['light_regime'].map(regime_to_freq)

        mutant_means = gene_data.groupby(['mutant_ID', 'light_regime']).agg(
            mean_y2_slope=('y2_slope', 'mean')
        ).reset_index()

        gene_agg = mutant_means.groupby('light_regime').agg(
            mean_y2_slope=('mean_y2_slope', 'mean')
        ).reset_index()

        gene_agg['freq'] = gene_agg['light_regime'].map(regime_to_freq)

        wt_agg = wt_data.groupby('light_regime').agg(
            mean_y2_slope=('y2_slope', 'mean')
        ).reset_index()
        wt_agg['freq'] = wt_agg['light_regime'].map(regime_to_freq)

        x_gene = gene_agg['freq'].values
        y_gene = gene_agg['mean_y2_slope'].values

        x_wt = wt_agg['freq'].values
        y_wt = wt_agg['mean_y2_slope'].values
        gene_agg = gene_agg.sort_values("freq")
        wt_agg = wt_agg.sort_values("freq")


        plt.figure(figsize=(10, 6))
        mutants = gene_data['mutant_ID'].unique()
        colors = cm.get_cmap('tab10', len(mutants)) 
        for i, mutant in enumerate(mutants):
            mutant_df = gene_data[gene_data['mutant_ID'] == mutant]
            plt.scatter(
                mutant_df['freq'],
                mutant_df['y2_slope'],
                label=mutant,
                color=colors(i),
                alpha=0.7
            )

        plt.xscale('log')

        plt.scatter(x_gene, y_gene, color='blue', label="Gene mean", zorder=3)
        plt.plot(gene_agg['freq'], gene_agg['mean_y2_slope'], color='orange', label="Gene mean trend", zorder=2)

        plt.scatter(x_wt, y_wt, color='black', label="WT mean", zorder=3, marker='x')
        plt.plot(wt_agg['freq'], wt_agg['mean_y2_slope'], color='black', linestyle='--', label="WT mean trend", zorder=1)

        # equidistant
        # plt.xticks(ticks=list(regime_to_freq.values()), labels=regime_order, rotation=45)
        # plt.xlabel("Light Regime (low to high frequency)")
        
        # freq scale
        plt.xticks(
            ticks=[regime_to_freq[r] for r in regime_order],
            labels=regime_order,
            rotation=45
        )
        plt.xlabel("Light switching frequency (1/min)")

        plt.ylabel("y2_slope")
        plt.title(f"Mutant vs WT Slope Trend – Gene: {gene}")
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.show()
        
