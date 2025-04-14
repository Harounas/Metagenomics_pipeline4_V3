import pandas as pd
import random
from collections import defaultdict
import plotly.express as px
import plotly.io as pio
import os
from .trimmomatic import run_trimmomatic
from .metaspades import run_spades
from .bowtie2 import run_bowtie2
from .kraken2 import run_kraken2
import distinctipy
import numpy as np
import matplotlib.pyplot as plt
import logging

def process_sample(forward, reverse, base_name, bowtie2_index, kraken_db, output_dir, threads, run_bowtie, use_precomputed_reports, use_assembly, skip_preprocessing=False):
    try:
        if use_precomputed_reports:
            kraken_report = os.path.join(output_dir, f"{base_name}_kraken_report.txt")
            if not os.path.exists(kraken_report):
                raise FileNotFoundError(f"Precomputed Kraken2 report not found: {kraken_report}")
            return kraken_report

        if skip_preprocessing:
            # Skip trimming, host depletion, and initial Kraken2 classification
            logging.info(f"Skipping preprocessing for sample {base_name}")
            contigs_file = os.path.join(output_dir, f"{base_name}_contigs.fasta")
            if not os.path.exists(contigs_file):
                logging.info(f"Performing de novo assembly for sample {base_name}")
                contigs_file = run_spades(forward, reverse, base_name, output_dir, threads)
            kraken_input = contigs_file
        else:
            # Existing preprocessing steps
            # Step 1: Run Trimmomatic
            trimmed_forward = os.path.join(output_dir, f"{base_name}_1_trimmed_paired.fq.gz")
            trimmed_reverse = os.path.join(output_dir, f"{base_name}_2_trimmed_paired.fq.gz")
            if not (os.path.exists(trimmed_forward) and os.path.exists(trimmed_reverse)):
                logging.info(f"Running Trimmomatic for sample {base_name}")
                trimmed_forward, trimmed_reverse = run_trimmomatic(
                    forward, reverse, base_name, output_dir, threads
                )
            else:
                logging.info(f"Using existing trimmed files for sample {base_name}")

            # Step 2: Optional host depletion
            if run_bowtie:
                bowtie_unmapped_r1 = os.path.join(output_dir, f"{base_name}_1_unmapped.fq.gz")
                bowtie_unmapped_r2 = os.path.join(output_dir, f"{base_name}_2_unmapped.fq.gz")
                if not (os.path.exists(bowtie_unmapped_r1) and os.path.exists(bowtie_unmapped_r2)):
                    logging.info(f"Running Bowtie2 host depletion for sample {base_name}")
                    unmapped_r1, unmapped_r2 = run_bowtie2(
                        trimmed_forward, trimmed_reverse, base_name, bowtie2_index, output_dir, threads
                    )
                else:
                    logging.info(f"Using existing Bowtie2 output files for sample {base_name}")
                    unmapped_r1, unmapped_r2 = bowtie_unmapped_r1, bowtie_unmapped_r2
            else:
                unmapped_r1, unmapped_r2 = trimmed_forward, trimmed_reverse

            # Step 3: Perform assembly if requested
            if use_assembly:
                contigs_file = os.path.join(output_dir, f"{base_name}_contigs.fasta")
                if not os.path.exists(contigs_file):
                    logging.info(f"Performing de novo assembly for sample {base_name}")
                    contigs_file = run_spades(unmapped_r1, unmapped_r2, base_name, output_dir, threads)
                kraken_input = contigs_file
            else:
                kraken_input = (unmapped_r1, unmapped_r2)

        # Step 4: Run Kraken2
        kraken_report = os.path.join(output_dir, f"{base_name}_kraken_report.txt")
        if not os.path.exists(kraken_report):
            if isinstance(kraken_input, tuple):
                input_r1, input_r2 = kraken_input
            else:
                input_r1, input_r2 = kraken_input, None
            kraken_report = run_kraken2(
                input_r1, input_r2, base_name, kraken_db, output_dir, threads
            )
        else:
            logging.info(f"Using existing Kraken2 report for sample {base_name}")

        return kraken_report

    except Exception as e:
        logging.error(f"Error processing sample {base_name}: {e}")
        return None

def generate_sample_ids_csv(kraken_dir):
    """
    Generates a CSV file containing sample IDs extracted from Kraken report filenames.

    Parameters:
    - kraken_dir (str): Path to the directory containing Kraken report files.

    Returns:
    - str: Path to the generated sample_ids.csv file.
    """
    try:
        sample_ids = []
        for fname in os.listdir(kraken_dir):
            if fname.endswith('_kraken_report.txt'):
                sample_id = fname.replace('_kraken_report.txt', '')
                sample_ids.append(sample_id)

        sample_ids_df = pd.DataFrame({'Sample_ID': sample_ids})
        csv_path = os.path.join(kraken_dir, 'sample_ids.csv')
        sample_ids_df.to_csv(csv_path, index=False)
        logging.info(f"Sample IDs written to {csv_path}")
        return csv_path

    except Exception as e:
        logging.error(f"Error generating sample IDs CSV: {e}")
        return None

def aggregate_kraken_results(kraken_dir, metadata_file=None, sample_id_df=None, read_count=1,max_read_count=1000000000000000000000000000):
    """
    Aggregates Kraken results, merging metadata or using sample IDs if metadata is not provided.

    Parameters:
    - kraken_dir (str): Path to the directory containing Kraken report files.
    - metadata_file (str, optional): Path to the metadata CSV file. Defaults to None.
    - sample_id_df (DataFrame, optional): DataFrame of sample IDs. Used if metadata_file is not provided.
    - read_count (int): Minimum read count threshold for filtering results.
    - read_count (int): Maximum read count threshold for filtering results.
    
    Returns:
    - str: Path to the generated merged TSV file.
    """
    try:
        # Load metadata from file first, fall back to sample_id_df if not provided
        if metadata_file:
            metadata = pd.read_csv(metadata_file, sep=",")
            print("Using metadata from the provided metadata file.")
        elif sample_id_df is not None:
            metadata = sample_id_df
            print("Using sample IDs as metadata.")
        else:
            raise ValueError("Either metadata_file or sample_id_df must be provided.")

        sample_id_col = metadata.columns[0]  # Assume the first column is the sample ID

        # Dictionary to store aggregated results
        aggregated_results = {}

        # Iterate over each Kraken report file
        for file_name in os.listdir(kraken_dir):
            if file_name.endswith("_report.txt"):
                with open(os.path.join(kraken_dir, file_name), 'r') as f:
                    for line in f:
                        fields = line.strip().split('\t')
                        perc_frag_cover = fields[0]
                        nr_frag_cover = fields[1]
                        nr_frag_direct_at_taxon = int(fields[2])
                        rank_code = fields[3]
                        ncbi_ID = fields[4]
                        scientific_name = fields[5]
                        parts = file_name.split('_')
                        extracted_part = '_'.join(parts[:-2])
                        sampleandtaxonid = extracted_part + str(ncbi_ID)

                        # Check if rank code is species-level and meets the read count threshold
                        if (rank_code == 'S' or rank_code == 'S1' or rank_code == 'S2' or rank_code == 'S3') and (nr_frag_direct_at_taxon >= read_count and nr_frag_direct_at_taxon <=max_read_count):
                            if extracted_part in metadata[sample_id_col].unique():
                                sample_metadata = metadata.loc[metadata[sample_id_col] == extracted_part].iloc[0].to_dict()
                                aggregated_results[sampleandtaxonid] = {
                                    'Perc_frag_cover': perc_frag_cover,
                                    'Nr_frag_cover': nr_frag_cover,
                                    'Nr_frag_direct_at_taxon': nr_frag_direct_at_taxon,
                                    'Rank_code': rank_code,
                                    'NCBI_ID': ncbi_ID,
                                    'Scientific_name': scientific_name,
                                    'SampleID': extracted_part,
                                    **sample_metadata
                                }

        # Output aggregated results to a TSV file
        merged_tsv_path = os.path.join(kraken_dir, "merged_kraken.tsv")
        with open(merged_tsv_path, 'w') as f:
            # Write headers dynamically
            headers = ['Perc_frag_cover', 'Nr_frag_cover', 'Nr_frag_direct_at_taxon', 'Rank_code', 'NCBI_ID', 'Scientific_name', 'SampleID'] + metadata.columns[1:].tolist()
            f.write("\t".join(headers) + "\n")
            for sampleandtaxonid, data in aggregated_results.items():
                f.write("\t".join(str(data[col]) for col in headers) + "\n")

        return merged_tsv_path

    except Exception as e:
        print(f"Error aggregating Kraken results: {e}")
        return None



def generate_abundance_plots(merged_tsv_path, top_N,col_filter,pat_to_keep):
    try:
        df = pd.read_csv(merged_tsv_path, sep="\t")
        df.columns = df.columns.str.replace('/', '_').str.replace(' ', '_')
        df = df.apply(lambda col: col.map(lambda x: x.strip() if isinstance(x, str) else x))
        df = df[df['Scientific_name'] != 'Homo sapiens']  # Remove human reads
        if col_filter:
            df=df[~df['Scientific_name'].isin(col_filter)] 
        if pat_to_keep:
            df=df[df['Scientific_name'].isin(pat_to_keep)] 
        #if max_read_count:
            #df=df[df['Nr_frag_direct_at_taxon']<=max_read_count]
            
            
        # Generate both viral and bacterial abundance plots
        for focus, filter_str, plot_title in [
            ('Virus_Type', 'Virus', 'Viral'),
            ('Bacteria_Type', 'Virus', 'Bacterial')
        ]:
            if focus == 'Bacteria_Type':
                df_focus = df[~df['Scientific_name'].str.contains(filter_str, case=False, na=False)]
                #summary_csv_path = os.path.join("Bacteria_summary.csv")
                #df_focus.to_csv(summary_csv_path , index=False)
            else:
                df_focus = df[df['Scientific_name'].str.contains(filter_str, case=False, na=False)]
                #summary_csv_path = os.path.join("Virus_summary.csv")
                #df_focus.to_csv(summary_csv_path , index=False)
            df_focus = df_focus.rename(columns={'Scientific_name': focus})

            if top_N:
                top_N_categories = df_focus[focus].value_counts().head(top_N).index
                df_focus = df_focus[df_focus[focus].isin(top_N_categories)]

            categorical_cols = df_focus.select_dtypes(include=['object']).columns.tolist()
            categorical_cols.remove(focus)

            for col in categorical_cols:
                grouped_sum = df_focus.groupby([focus, col])['Nr_frag_direct_at_taxon'].mean().reset_index()
                # Create a color mapping based on unique values in the 'focus' column
                #colordict = dict(zip(grouped_sum[focus].unique(), distinctipy.get_colors(len(grouped_sum[focus].unique()))))
                colordict = defaultdict(int)
                random_colors0 = ["#{:06X}".format(random.randint(0, 0xFFFFFF)) for _ in range(len(grouped_sum[focus].unique()))]
                #random_colors = ['#{:06x}'.format(random.randint(000000, 0xffff00)) for _ in range(len(grouped_sum[focus].unique()))]
                random_colors1 =['#000000','#FF0000','#556B2F','#ADD8E6','#6495ED','#00FF00','#0000FF','#FFFF00','#00FFFF','#FF00FF','#C0C0C0','#808080','#800000','#808000','#008000',
'#008080','#000080','#CD5C5C','#DAA520','#FFA500','#F0E68C','#ADFF2F','#2F4F4F','#E0FFFF','#4169E1','#8A2BE2','#4B0082','#EE82EE','#D2691E','#BC8F8F','#800080','#DDA0DD','#FF1493','#8B4513','#A0522D','#708090','#B0C4DE','#FFFFF0','#DCDCDC','#FFEFD5','#F5DEB3','#7FFFD4','#FFC0CB','#A52A2A','#040720','#34282C','#3B3131','#3A3B3C','#52595D','#FFFFFF','#FFFFF4','#FFF9E3']              #for target, color in zip(grouped_sum[focus].unique(), random_colors):
                if (len(grouped_sum[focus].unique())<=len(random_colors1)):
                  for target, color in zip(grouped_sum[focus].unique(), random_colors1[:len(grouped_sum[focus].unique())]):
                    colordict[target] = color
                else :
                  for target, color in zip(grouped_sum[focus].unique(), random_colors0):
                    colordict[target] = color
               
                #colordict=distinctipy.get_colors(len(grouped_sum[col].unique()))
                #colordict = dict(zip(grouped_sum[focus].unique(), distinctipy.get_colors(len(grouped_sum[focus].unique()))))
                # Generate a unique color for each unique item in the 'focus' column
                #random_colors = distinctipy.get_colors(len(grouped_sum[col].unique()))
                #colordict = {focus: color for category, color in zip(grouped_sum[col].unique(), random_colors)}
                def color_distance(c1, c2):
                  return np.sqrt(sum((a - b) ** 2 for a, b in zip(c1, c2)))

                # Function to generate colors with a minimum distance
               # Function to generate unique colors with a minimum distance
                def generate_distant_colors(num_colors, min_distance=30):
                   colors = []
                   seen_colors = set()  # To track already generated colors in tuple format
                   while len(colors) < num_colors:
                      new_color = (random.randint(0, 255), random.randint(0, 255), random.randint(0, 255))
        
                       # Check if the new color is sufficiently distant and not already seen
                      if all(color_distance(new_color, existing) >= min_distance for existing in colors) and new_color not in seen_colors:
                         colors.append(new_color)
                         seen_colors.add(new_color)  # Add the new color to the set
                   # Convert RGB colors to hex format
                   hex_colors = [f"#{r:02x}{g:02x}{b:02x}" for r, g, b in colors]
                   return hex_colors
                #num_categories = len(grouped_sum[focus].unique())
                # Generate distinct colors for each category
               # colors = generate_distant_colors(num_categories, min_distance=90)
                #colors = base_colors[:len(unique_targets)]
                #colors = plt.get_cmap('tab20').colors 
                #colordict = dict(zip(grouped_sum[focus].unique(), colors))
                # Generate colors using 'tab20' colormap
                #colors = plt.get_cmap('tab20').colors

                  # Ensure there are enough colors for unique items in grouped_sum[focus]
                #unique_targets = grouped_sum[focus].unique()
                #if len(unique_targets) > len(base_colors):
                  #raise ValueError("Not enough unique colors in 'tab20' colormap for the number of unique targets.")
                  #colors = plt.get_cmap('tab20')(np.linspace(0, 1, len(unique_targets)))
                #else:
                # Use base 'tab20' colors if they are sufficient
                 # colors = base_colors[:len(unique_targets)]
                 # Map each unique target to a color from the colormap
                #olordict = dict(zip(grouped_sum[focus].unique(), colors[:len(unique_targets)]))
                plot_width = 1100 + 5 * len(grouped_sum[col].unique())
                plot_height = 800 + 5 * len(grouped_sum[col].unique())
                font_size = max(10, 14 - len(grouped_sum[col].unique()) // 10)

                fig = px.bar(
                    grouped_sum,
                    x=col,
                    y='Nr_frag_direct_at_taxon',
                    color=focus,
                    color_discrete_map=colordict,
                    title=f"{plot_title} Abundance by {col}"
                )
                summary_csv_path = os.path.join(f"{plot_title}_summary.csv")
                grouped_sum.to_csv(summary_csv_path , index=False)
                fig.update_layout(
                    xaxis=dict(tickfont=dict(size=font_size), tickangle=45),
                    yaxis=dict(tickfont=dict(size=font_size)),
                    title=dict(text=f'Average {plot_title} Abundance by {col}', x=0.5, font=dict(size=16)),
                    bargap=0.5,
                    legend=dict(
                        font=dict(size=font_size),
                        x=1,
                        y=1,
                        traceorder='normal',
                        orientation='v',
                        itemwidth=30,
                        itemsizing='constant',
                        itemclick='toggleothers',
                        itemdoubleclick='toggle'
                    ),
                    width=plot_width,
                    height=plot_height
                )

                fig.write_image(f"{plot_title}_Abundance_by_{col}.png", format='png', scale=3,width=1920, height=1080)
           

    except Exception as e:
        print(f"Error generating abundance plots: {e}")

