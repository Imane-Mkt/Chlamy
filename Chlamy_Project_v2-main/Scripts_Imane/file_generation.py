import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.cm as cm
import statsmodels.api as sm
from statsmodels.formula.api import ols

def pdf_slopetrend(data, genes, chunk_size=60, out_prefix="Gene_Plots_Report"):
    regime_to_freq = {
        '20h_HL': 10000,                     
        '2h-2h': 120,                
        '10min-10min': 10,          
        '5min-5min': 5,
        '1min-1min': 1,
        '30s-30s': 0.5           
    }
    regime_order = list(regime_to_freq.keys())

    # loop through chunks of the gene list
    for i in range(0, len(genes), chunk_size):
        chunk_genes = genes[i:i+chunk_size]
        pdf_name = f"{out_prefix}_part{i//chunk_size+1}.pdf"
        with PdfPages(pdf_name) as pdf:
            for gene in chunk_genes:
                gene_data = data[
                    (data['mutated_genes'].astype(str).str.contains(gene)) &
                    (data['light_regime'].isin(regime_order)) &
                    (~data['y2_slope'].isna())
                ].copy()
                if gene_data.empty:
                    continue

                gene_data['freq'] = gene_data['light_regime'].map(regime_to_freq)
                plates = gene_data['plate'].unique()
                wt_data = data[
                    (data['mutated_genes'] == 'WT') &
                    (data['plate'].isin(plates)) &
                    (data['light_regime'].isin(regime_order)) &
                    (~data['y2_slope'].isna())
                ].copy()
                wt_data['freq'] = wt_data['light_regime'].map(regime_to_freq)

                # Aggregates
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

                gene_agg = gene_agg.sort_values("freq")
                wt_agg = wt_agg.sort_values("freq")

                # --- Plot ---
                plt.figure(figsize=(10, 6))
                mutants = gene_data['mutant_ID'].unique()
                colors = cm.get_cmap('tab10', len(mutants)) 
                for j, mutant in enumerate(mutants):
                    mutant_df = gene_data[gene_data['mutant_ID'] == mutant]
                    plt.scatter(
                        mutant_df['freq'],
                        mutant_df['y2_slope'],
                        label=mutant,
                        color=colors(j),
                        alpha=0.7
                    )

                plt.xscale('log')
                plt.scatter(gene_agg['freq'], gene_agg['mean_y2_slope'], color='blue', label="Gene mean")
                plt.plot(gene_agg['freq'], gene_agg['mean_y2_slope'], color='orange', label="Gene mean trend")
                plt.scatter(wt_agg['freq'], wt_agg['mean_y2_slope'], color='black', marker='x', label="WT mean")
                plt.plot(wt_agg['freq'], wt_agg['mean_y2_slope'], color='black', linestyle='--', label="WT trend")

                plt.xticks(
                    ticks=[regime_to_freq[r] for r in regime_order],
                    labels=regime_order,
                    rotation=45
                )
                plt.xlabel("Light switching frequency (1/min)")
                plt.ylabel("y2_slope")
                plt.title(f"Mutant vs WT – Gene: {gene}")
                plt.grid(True)
                plt.legend()
                plt.tight_layout()

                pdf.savefig()
                plt.close()

                # # --- ANCOVA results page ---
                # gene_data['is_mutant'] = 1
                # wt_data['is_mutant'] = 0
                # combined_data = pd.concat([gene_data, wt_data], ignore_index=True)

                # if not combined_data.empty:
                #     model = ols('y2_slope ~ freq * is_mutant', data=combined_data).fit()
                #     anova_table = sm.stats.anova_lm(model, typ=2)

                #     fig = plt.figure(figsize=(8, 4))
                #     plt.axis('off')
                #     text = f"ANCOVA results for {gene}:\n\n{anova_table.to_string()}"
                #     fig.text(0.01, 0.9, text, fontsize=10, va='top')
                #     pdf.savefig(fig)
                #     plt.close(fig)

        print(f"Saved {pdf_name} with {len(chunk_genes)} genes")
