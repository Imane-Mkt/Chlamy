import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

data_path = pd.read_csv('/home/imokhtatif/.vscode-server/Chlamy_Project_v2-main/Data/database_2025-03-25.csv')

def sanity_checks(data):
    """
    This function performs sanity checks on the data.
    It checks for missing values, duplicates, data types, dates to distinguish between phases.
    """
    data['start_date'] = pd.to_datetime(data['start_date'], errors='coerce')
    earliest_date = data['start_date'].min()
    most_recent_date = data['start_date'].max()
    print(f"Earliest start date : {earliest_date}")
    print(f"Most recent start date : {most_recent_date}")
    
    light_regimes = data['light_regime'].unique()
    print(f"Light regimes : {light_regimes}")
    
    print(f"Number of unique plates : {data['plate'].nunique()}")
    print(f"Unique plates : {data['plate'].unique()}")
    
    # Light regimes per plate
    from collections import Counter
    plate_light_regime_counts = {}
    for plate_number in data['plate'].unique():
        plate_data = data[data['plate'] == plate_number]
        if 'light_regime' in plate_data.columns:
            unique_light_regimes = plate_data['light_regime'].unique()
            plate_light_regime_counts[plate_number] = len(unique_light_regimes)
    regime_count_distribution = Counter(plate_light_regime_counts.values())
    for regime_count, plate_count in regime_count_distribution.items():
        print(f"{plate_count} plates have {regime_count} unique light regimes.")

def get_format_data(data):
    """
    This function formats the data for plotting.
    Drop NaN lines and missing data lines
    Filter out Nan y2 values > max non NaN index (the time series have != lengths)
    Format mutated_genes into a list and put != genes in different rows
    Add elapsed_hours (measurement_time) to plot the time series
    """    
     # Delete Nan Columns
    #print(data.shape,"Colonnes après 2 :", data.columns)
    
    # Set mutated_genes to NaN for all WT or WT CLiP rows
    data.loc[data['mutant_ID'].isin(['WT', 'WT_CLiP']), 'mutated_genes'] = 'WT'
    # 2. Drop WT rows with well_id starting with 'P' on plate starting with '31'
    condition_drop_p_wt = (
        (data['plate'].astype(str).str.startswith('31')) &
        (data['mutant_ID'].isin(['WT', 'WT_CLiP'])) &
        (data['well_id'].astype(str).str.startswith('P'))
    )
    data = data[~condition_drop_p_wt]
    
    # Create a new column 'fv_fm_end' and initialize it with NaN
    data['fv_fm_end'] = None
    data['ynpqend'] = None
    ynpq_columns = [col for col in data.columns if col.startswith('ynpq_')]
    y2_columns = [col for col in data.columns if col.startswith('y2_') and 'std' not in col]
    time_columns = [col for col in data.columns if col.startswith('measurement_time_')]
    data['fv_fm_end'] = data[y2_columns].apply(lambda row: row.dropna().iloc[-1] if not row.dropna().empty else np.nan, axis=1)
    data['ynpqend'] = data[ynpq_columns].apply(lambda row: row.dropna().iloc[-1] if not row.dropna().empty else np.nan, axis=1)
    # for index, row in data.iterrows():
    #     valid_time_indices = row[time_columns].notna()
    #     valid_y2_indices = row[y2_columns].notna()
    #     valid_ynpq_indices = row[ynpq_columns].notna()
    #     common_valid_indices = valid_time_indices & valid_y2_indices & valid_ynpq_indices
    #     for col in time_columns:         # Mettre à NaN les valeurs en dehors des indices communs
    #         if not common_valid_indices[col]:
    #             data.at[index, col] = np.nan
    #     for col in y2_columns:
    #         if not common_valid_indices[col]:
    #             data.at[index, col] = np.nan
    #     for col in ynpq_columns:
    #         if not common_valid_indices[col]:
    #             data.at[index, col] = np.nan
            
    # Supprimer les colonnes sans nom
    if None in data.columns:
        #print("Colonne 'None' détectée, suppression...")
        data = data.loc[:, data.columns.notnull()]

    #print(data.shape,"Colonnes après 0 :", data.columns)

    # Drop y2 Nan lines
    y2_columns = [col for col in data.columns if col.startswith('y2_') and col.split('_')[-1].isdigit() and 'std' not in col]
    rows_no_y2_with_mutant_id = data[data[y2_columns].isna().all(axis=1)]
    data = data.drop(rows_no_y2_with_mutant_id.index)
    #print(data.shape,"Colonnes après 3 :", data.columns)
    # Drop missing mutant ID lines
    rows_to_delete = data[data[y2_columns].notna().any(axis=1) & data['mutant_ID'].isna()]
    data = data.drop(rows_to_delete.index)
    #print(data.shape,"Colonnes après 4 :", data.columns)

    # Add elapsed_hours (measurement_time) to plot the time series
    time_columns = [col for col in data.columns if col.startswith('measurement_time_')]
    #print("Colonnes de temps disponibles :", time_columns)
    for col in time_columns:
        #print(f"Colonne {col} - Valeurs manquantes :", data[col].isnull().sum())
        data[col] = pd.to_datetime(data[col], format= 'mixed')
    elapsed_hours = data[time_columns].apply(lambda row: (row - row.min()).dt.total_seconds() / 3600, axis=1)
    elapsed_hours.columns = [f'elapsed_hours_{i}' for i in range(len(elapsed_hours.columns))]
    # Find the last index of a measurement_time_ column
    time_col_indices = [data.columns.get_loc(col) for col in time_columns]
    insert_at = max(time_col_indices) + 1
    # Insert each elapsed_hours column at the right position
    for i, col in enumerate(elapsed_hours.columns):
        data.insert(insert_at + i, col, elapsed_hours.iloc[:, i])
    
    # Remove suspicious plates
    data = data[data['plate']!='20']
    # Merge 99 and '99'
    data['plate'] = data['plate'].apply(lambda x: '99' if x == 99 else x)

    # Put genes of the same mutant in differnt rows
    # rows = []
    # for _, row in data.iterrows():
    #     genes = [g.strip() for g in str(row['mutated_genes']).split(',')]#replace('&', ',').split(',')]
    #     features = [f.strip() for f in str(row['feature']).split(',')]#replace('&', ',').split(',')]
    #     if len(genes) != len(features):
    #         for gene in genes:
    #             new_row = row.copy()
    #             new_row['mutated_genes'] = gene
    #             rows.append(new_row)
    #     else:
    #         for gene, feature in zip(genes, features):
    #             new_row = row.copy()
    #             new_row['mutated_genes'] = gene
    #             new_row['feature'] = feature
    #             rows.append(new_row)
    # data = pd.DataFrame(rows).reset_index(drop=True)
    # if not rows:
    #     print("La liste 'rows' est vide. Aucune donnée n'a été ajoutée.")
    
    # Delete the High pulse end points in continuous light
    # Supprimer le dernier point valide (non-NaN) dans chaque time series y2_
    for idx, row in data.iterrows():
        valid_cols = [col for col in y2_columns if pd.notna(row[col])]
        if valid_cols:
            last_col = valid_cols[-1]
            data.at[idx, last_col] = np.nan
    ynpq_columns = [col for col in data.columns if col.startswith('ynpq_') and col.split('_')[-1].isdigit()]
    for idx, row in data.iterrows():
        valid_cols = [col for col in ynpq_columns if pd.notna(row[col])]
        if valid_cols:
            last_col = valid_cols[-1]
            data.at[idx, last_col] = np.nan
    # Supprimer le dernier point valide dans chaque time series measurement_time_ et elapsed_hours_
    measurement_time_columns = [col for col in data.columns if col.startswith('measurement_time_') and col.split('_')[-1].isdigit()]
    elapsed_hours_columns = [col for col in data.columns if col.startswith('elapsed_hours_') and col.split('_')[-1].isdigit()]

    for idx, row in data.iterrows():
        # Pour measurement_time_
        valid_cols = [col for col in measurement_time_columns if pd.notna(row[col])]
        if valid_cols:
            last_col = valid_cols[-1]
            data.at[idx, last_col] = np.nan

        # Pour elapsed_hours_
        valid_cols = [col for col in elapsed_hours_columns if pd.notna(row[col])]
        if valid_cols:
            last_col = valid_cols[-1]
            data.at[idx, last_col] = np.nan
            
        # Keep the non Nan y2 values and shave the rest
    y2_columns = [col for col in data.columns if col.startswith('y2_') and col.split('_')[-1].isdigit() and 'std' not in col]
    light_regimes = data['light_regime'].unique()
    max_timesteps_per_light = {}
    for light in light_regimes:
        data_light = data[data['light_regime'] == light]
        max_timestep = 0
        for col in y2_columns:
            if data_light[col].notna().any():  # non-NaN
                max_timestep = max(max_timestep, int(col.split('_')[-1]))
        max_timesteps_per_light[light] = max_timestep
        max_timesteps_df = pd.DataFrame(list(max_timesteps_per_light.items()), columns=['light_regime', 'max_timestep'])
    #print(data.shape, max_timesteps_df)
    columns_to_drop = [
    col for col in data.columns if (
        (col.startswith('y2_') or col.startswith('y2_std_') or col.startswith('ynpq_') or col.startswith('measurement_time_') or col.startswith('elapsed_hours_'))
        and col.split('_')[-1].isdigit()
        and int(col.split('_')[-1]) > max(max_timesteps_df['max_timestep']) 
    )]
    #print(data.shape,"Colonnes après 1 :", data.columns)
    data = data.drop(columns=columns_to_drop)
    #print  (data.shape, "Colonnes de données après formatage :", data.columns)
    formatted_data = data.copy()
    return formatted_data

def remove_oldbad_plates(data):
    """
    Removes duplicate entries with the same plate, light_regime, mutant_ID, and well_id,
    keeping only the one with the most recent start_date.
    """
    print('myea')
    data['start_date'] = pd.to_datetime(data['start_date'])
    deduped = data.sort_values('start_date').drop_duplicates(
        subset=['plate', 'light_regime', 'mutant_ID', 'mutated_genes', 'well_id'],
        keep='last'
    )
    return deduped