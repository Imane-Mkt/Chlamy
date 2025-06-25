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
def plot_slope_correlation(data):
    # Pivot to wide format: one row per mutant_ID, one column per light regime
    pivoted = data.pivot_table(
        index="mutant_ID",
        columns="light_regime",
        values="y2_slope",
        aggfunc="mean"
    )
    # Define desired light regime order
    regime_order = [
        "20h_HL", "2h-2h",
        "10min-10min", "5min-5min", "1min-5min", "1min-1min", "30s-30s"
    ]

    # Reorder columns after pivot
    pivoted = pivoted[regime_order]

    # Drop rows with too many missing values (optional)
    pivoted = pivoted.dropna(thresh=3)  # keep rows with at least 3 regimes

    # Compute correlation
    corr_matrix = pivoted.corr()

    # Plot heatmap
    plt.figure(figsize=(10, 8))
    sns.heatmap(corr_matrix, cmap="magma", annot=True, fmt=".2f",
                vmin=0, vmax=1, linewidths=0.5)
    plt.title("Correlation of average y2_slope between light regimes")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.show()
