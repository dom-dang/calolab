#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import re
import numpy as np

# Map contigs to chromosome names
contig_to_chr = {
    "CM000663.2": "chr1", "CM000664.2": "chr2", "CM000665.2": "chr3",
    "CM000666.2": "chr4", "CM000667.2": "chr5", "CM000668.2": "chr6",
    "CM000669.2": "chr7", "CM000670.2": "chr8", "CM000671.2": "chr9",
    "CM000672.2": "chr10", "CM000673.2": "chr11", "CM000674.2": "chr12",
    "CM000675.2": "chr13", "CM000676.2": "chr14", "CM000677.2": "chr15",
    "CM000678.2": "chr16", "CM000679.2": "chr17", "CM000680.2": "chr18",
    "CM000681.2": "chr19", "CM000682.2": "chr20", "CM000683.2": "chr21",
    "CM000684.2": "chr22", "CM000685.2": "chrX", "CM000686.2": "chrY"
}

# Create a reverse mapping for lookup
chr_to_contig = {v: k for k, v in contig_to_chr.items()}

def load_trna_data_from_csv():
    """
    Load tRNA data from the provided CSV file
    """
    try:
        # Load tRNA gene data from the CSV file
        trna_data = pd.read_csv("human_tRNA_data.csv")
        
        # Count tRNA genes by chromosome
        trna_counts = trna_data['chromosome'].value_counts().reset_index()
        trna_counts.columns = ['chromosome', 'trna_gene_count']
        
        # Make sure we have entries for all standard chromosomes
        all_chromosomes = set(contig_to_chr.values())
        existing_chromosomes = set(trna_counts['chromosome'])
        
        # Add missing chromosomes with count 0
        missing_chromosomes = all_chromosomes - existing_chromosomes
        for chrom in missing_chromosomes:
            trna_counts = pd.concat([trna_counts, pd.DataFrame({'chromosome': [chrom], 'trna_gene_count': [0]})], ignore_index=True)
        
        # Sort chromosomes in natural order
        def chr_sort_key(chr_name):
            if chr_name == 'chrM':
                return 100 
            elif chr_name == 'chrX':
                return 23
            elif chr_name == 'chrY':
                return 24
            else:
                num = re.sub(r'chr', '', chr_name)
                try:
                    return int(num)
                except ValueError:
                    return 99  # For any other non-standard chromosome
        
        trna_counts['sort_key'] = trna_counts['chromosome'].apply(chr_sort_key)
        trna_counts = trna_counts.sort_values('sort_key').drop('sort_key', axis=1)
        
        return trna_counts
    
    except Exception as e:
        print(f"Error loading tRNA data from CSV: {e}")
        # Return an empty DataFrame with the correct columns
        return pd.DataFrame(columns=['chromosome', 'trna_gene_count'])

def load_read_mapping_data():
    """
    Load read mapping data 
    """
    try:
        # Process the first dataset: R24-5134
        df1 = pd.read_csv("R24-5134_bwa_v1.tsv", sep='\t')
        df1.columns = ['reference', 'count']
        
        # Process the second dataset: r24-5135
        df2 = pd.read_csv("r24-5135_bwa_v1.tsv", sep='\t')  
        df2.columns = ['reference', 'count']
        
        # Combine the two datasets
        mapping_df = pd.DataFrame({
            'reference': df1['reference'],
            'count_5134': df1['count'],
            'count_5135': df2['count']
        })
        
        # Calculate total counts
        mapping_df['total_count'] = mapping_df['count_5134'] + mapping_df['count_5135']
        
        # Filter to only keep the main chromosomes
        main_refs = [ref for ref in contig_to_chr.keys()]
        main_mapping_df = mapping_df[mapping_df['reference'].isin(main_refs)]
        
        # Map contig IDs to chromosome names
        main_mapping_df['chromosome'] = main_mapping_df['reference'].map(contig_to_chr)
        
        # Group by chromosome and sum counts (in case there are multiple entries per chromosome)
        result_df = main_mapping_df.groupby('chromosome').agg({
            'count_5134': 'sum',
            'count_5135': 'sum',
            'total_count': 'sum'
        }).reset_index()
        
        # Rename to 'mapped_reads' for consistency with rest of code
        result_df = result_df.rename(columns={'total_count': 'mapped_reads'})
        
        return result_df
    
    except Exception as e:
        print(f"Error loading read mapping data: {e}")
        return pd.DataFrame(columns=['chromosome', 'count_5134', 'count_5135', 'mapped_reads'])

def calculate_and_plot_comparisons():
    """
    Calculate the fraction of mapped reads to tRNA gene counts
    and visualize the results
    """
    # Get tRNA gene counts from CSV
    trna_df = load_trna_data_from_csv()
    
    if trna_df.empty:
        print("No tRNA gene data available")
        return
    
    # Get mapped reads
    mapped_df = load_read_mapping_data()
    
    if mapped_df.empty:
        print("No mapped read data available")
        return
    
    # Merge datasets
    merged_df = pd.merge(trna_df, mapped_df, on='chromosome', how='outer')
    
    # Handle any missing values
    merged_df['trna_gene_count'] = merged_df['trna_gene_count'].fillna(0)
    merged_df['mapped_reads'] = merged_df['mapped_reads'].fillna(0)
    merged_df['count_5134'] = merged_df['count_5134'].fillna(0)
    merged_df['count_5135'] = merged_df['count_5135'].fillna(0)
    
    # Calculate the fraction (reads per gene)
    # Avoid division by zero by adding a small epsilon where gene count is 0
    merged_df['reads_per_gene'] = merged_df.apply(
        lambda row: row['mapped_reads'] / row['trna_gene_count'] if row['trna_gene_count'] > 0 else 0, 
        axis=1
    )
    
    # Add a sorting key to ensure chromosomes are shown in the correct order
    def chr_sort_key(chr_name):
        if chr_name == 'chrM':
            return 100  
        elif chr_name == 'chrX':
            return 23
        elif chr_name == 'chrY':
            return 24
        else:
            num = re.sub(r'chr', '', chr_name)
            try:
                return int(num)
            except ValueError:
                return 99 
    
    merged_df['sort_key'] = merged_df['chromosome'].apply(chr_sort_key)
    merged_df = merged_df.sort_values('sort_key')
    
    # Create ordered chromosome list for the x-axis
    ordered_chromosomes = merged_df['chromosome'].tolist()
    
    # Set up the figure with subplots - add two more for individual samples (now 7 plots total)
    fig, axs = plt.subplots(7, 1, figsize=(14, 42))
    
    # Plot 1: tRNA gene counts by chromosome
    sns.barplot(x='chromosome', y='trna_gene_count', data=merged_df, ax=axs[0], color='skyblue', order=ordered_chromosomes)
    axs[0].set_title('tRNA Gene Counts by Chromosome', fontsize=14)
    axs[0].set_xlabel('Chromosome', fontsize=12)
    axs[0].set_ylabel('Number of tRNA Genes', fontsize=12)
    axs[0].tick_params(axis='x', rotation=45)
    
    # Plot 2: Mapped reads by chromosome
    sns.barplot(x='chromosome', y='mapped_reads', data=merged_df, ax=axs[1], color='salmon', order=ordered_chromosomes)
    axs[1].set_title('Mapped Reads by Chromosome', fontsize=14)
    axs[1].set_xlabel('Chromosome', fontsize=12)
    axs[1].set_ylabel('Number of Mapped Reads', fontsize=12)
    axs[1].tick_params(axis='x', rotation=45)
    
    # Plot 3: Reads per gene by chromosome
    sns.barplot(x='chromosome', y='reads_per_gene', data=merged_df, ax=axs[2], color='lightgreen', order=ordered_chromosomes)
    axs[2].set_title('Reads per tRNA Gene by Chromosome', fontsize=14)
    axs[2].set_xlabel('Chromosome', fontsize=12)
    axs[2].set_ylabel('Reads per Gene', fontsize=12)
    axs[2].tick_params(axis='x', rotation=45)
    
    # Add a horizontal line at the mean in the third plot
    mean_reads_per_gene = merged_df['reads_per_gene'].mean()
    axs[2].axhline(y=mean_reads_per_gene, color='red', linestyle='--', 
                  label=f'Mean: {mean_reads_per_gene:.2f}')
    axs[2].legend()
    
    # Plot 4: Sample R24-5134 reads by chromosome
    sns.barplot(x='chromosome', y='count_5134', data=merged_df, ax=axs[3], color='purple', order=ordered_chromosomes)
    axs[3].set_title('Sample R24-5134 Reads by Chromosome', fontsize=14)
    axs[3].set_xlabel('Chromosome', fontsize=12)
    axs[3].set_ylabel('Number of Reads', fontsize=12)
    axs[3].tick_params(axis='x', rotation=45)
    
    # Plot 5: Sample r24-5135 reads by chromosome
    sns.barplot(x='chromosome', y='count_5135', data=merged_df, ax=axs[4], color='teal', order=ordered_chromosomes)
    axs[4].set_title('Sample r24-5135 Reads by Chromosome', fontsize=14)
    axs[4].set_xlabel('Chromosome', fontsize=12)
    axs[4].set_ylabel('Number of Reads', fontsize=12)
    axs[4].tick_params(axis='x', rotation=45)
    
    # Plot 6: Comparison between samples (grouped bar plot)
    sample_melted = pd.melt(merged_df, 
                           id_vars=['chromosome'], 
                           value_vars=['count_5134', 'count_5135'],
                           var_name='sample', 
                           value_name='value')
    
    sns.barplot(x='chromosome', y='value', hue='sample', 
               data=sample_melted,
               palette={'count_5134': 'purple', 'count_5135': 'teal'},
               order=ordered_chromosomes,
               ax=axs[5])
    
    axs[5].set_title('Sample Comparison: R24-5134 vs r24-5135 Reads by Chromosome', fontsize=14)
    axs[5].set_xlabel('Chromosome', fontsize=12)
    axs[5].set_ylabel('Number of Reads', fontsize=12)
    axs[5].tick_params(axis='x', rotation=45)
    axs[5].legend(title='Sample', labels=['R24-5134', 'r24-5135'])
    

    # Plot 7: Combined tRNA gene counts with both sample reads using dual Y axes
    # Plot the two read count datasets on the primary y-axis
    read_data = merged_df[['chromosome', 'count_5134', 'count_5135']].copy()
    read_data = pd.melt(read_data,
                      id_vars=['chromosome'],
                      value_vars=['count_5134', 'count_5135'],
                      var_name='sample',
                      value_name='read_count')
    
    # Create a mapping for better labels
    read_data['sample'] = read_data['sample'].replace({
        'count_5134': 'Sample R24-5134',
        'count_5135': 'Sample r24-5135'
    })
    
    # Plot read counts on the primary y-axis
    sns.barplot(x='chromosome', y='read_count', hue='sample',
               data=read_data,
               palette={'Sample R24-5134': 'purple', 'Sample r24-5135': 'teal'},
               order=ordered_chromosomes,
               ax=axs[6])
    
    # Create a twin axis for the tRNA gene counts
    ax_twin = axs[6].twinx()
    
    # Plot tRNA gene counts on the secondary y-axis
    sns.barplot(x='chromosome', y='trna_gene_count', 
               data=merged_df,
               color='skyblue',
               order=ordered_chromosomes,
               ax=ax_twin,
               alpha=0.7)
    
    # Set labels and titles
    axs[6].set_title('tRNA Gene Counts vs Sample Reads by Chromosome (Dual Y-axis)', fontsize=14)
    axs[6].set_xlabel('Chromosome', fontsize=12)
    axs[6].set_ylabel('Number of Reads', fontsize=12, color='darkblue')
    ax_twin.set_ylabel('Number of tRNA Genes', fontsize=12, color='skyblue')
    
    # Color the tick labels to match the respective datasets
    ax_twin.tick_params(axis='y', colors='skyblue')
    axs[6].tick_params(axis='y', colors='darkblue')
    axs[6].tick_params(axis='x', rotation=45)
    
    # Create a combined legend
    lines1, labels1 = axs[6].get_legend_handles_labels()
    # Create a custom handle for the tRNA gene counts
    from matplotlib.lines import Line2D
    custom_line = Line2D([0], [0], color='skyblue', lw=4)
    lines1.append(custom_line)
    labels1.append('tRNA Genes')
    
    # Remove the automatically created legend and create a new combined one
    axs[6].get_legend().remove()
    axs[6].legend(lines1, labels1, title='Data Type', loc='upper right')
    
    # Adjust layout
    plt.tight_layout()
    
    # Save the figure
    plt.savefig('trna_analysis_results.png', dpi=300)
    
    # Print summary statistics
    print("\nSummary Statistics:")
    print(f"Total tRNA genes: {merged_df['trna_gene_count'].sum()}")
    print(f"Total mapped reads: {merged_df['mapped_reads'].sum()}")
    print(f"Total reads in sample R24-5134: {merged_df['count_5134'].sum()}")
    print(f"Total reads in sample r24-5135: {merged_df['count_5135'].sum()}")
    
    print("\nChromosomes with highest tRNA gene density:")
    top_gene_density = merged_df.sort_values('trna_gene_count', ascending=False).head(5)
    print(top_gene_density[['chromosome', 'trna_gene_count']])
    
    print("\nChromosomes with highest reads per gene:")
    # Filter out chromosomes with zero tRNA genes to avoid infinity values
    valid_reads_df = merged_df[merged_df['trna_gene_count'] > 0]
    top_reads_per_gene = valid_reads_df.sort_values('reads_per_gene', ascending=False).head(5)
    print(top_reads_per_gene[['chromosome', 'reads_per_gene']])
    
    # Identify chromosomes with significantly different read ratios
    mean_ratio = valid_reads_df['reads_per_gene'].mean()
    std_ratio = valid_reads_df['reads_per_gene'].std()
    
    # Define threshold for "significant" difference (e.g., 1.5 standard deviations)
    threshold = 1.5 * std_ratio
    
    high_ratio_chrs = valid_reads_df[valid_reads_df['reads_per_gene'] > mean_ratio + threshold]
    low_ratio_chrs = valid_reads_df[valid_reads_df['reads_per_gene'] < mean_ratio - threshold]
    
    print("\nChromosomes with significantly higher reads per gene:")
    print(high_ratio_chrs[['chromosome', 'reads_per_gene', 'trna_gene_count', 'mapped_reads']])
    
    print("\nChromosomes with significantly lower reads per gene:")
    print(low_ratio_chrs[['chromosome', 'reads_per_gene', 'trna_gene_count', 'mapped_reads']])
    
    # Add sample-specific analysis
    print("\nSample comparison - R24-5134 vs r24-5135:")
    merged_df['ratio_5135_to_5134'] = merged_df.apply(
        lambda row: row['count_5135'] / row['count_5134'] if row['count_5134'] > 0 else float('inf'), 
        axis=1
    )
    
    # Replace infinity values with NaN for better reporting
    merged_df['ratio_5135_to_5134'].replace([float('inf'), float('-inf')], np.nan, inplace=True)
    
    # Find chromosomes with large differences between samples
    valid_ratio_df = merged_df.dropna(subset=['ratio_5135_to_5134'])
    high_diff_samples = valid_ratio_df.sort_values('ratio_5135_to_5134', ascending=False).head(5)
    print("\nChromosomes with higher r24-5135 to R24-5134 ratio:")
    print(high_diff_samples[['chromosome', 'count_5134', 'count_5135', 'ratio_5135_to_5134']])
    
    # Return the merged DataFrame for further analysis
    return merged_df

def main():
    """
    Main function to execute the workflow
    """
    print("Starting tRNA gene analysis...")
    
    # Run the analysis
    result_df = calculate_and_plot_comparisons()
    
    if result_df is not None:
        # Remove the sort_key column before saving
        if 'sort_key' in result_df.columns:
            result_df = result_df.drop('sort_key', axis=1)
            
        # Save the results to CSV
        result_df.to_csv('trna_analysis_results.csv', index=False)
        
        # Also save the individual sample results
        sample_df = result_df[['chromosome', 'count_5134', 'count_5135', 'mapped_reads', 'trna_gene_count', 'reads_per_gene']]
        sample_df.to_csv('trna_sample_comparison.csv', index=False)
        
        print("\nAnalysis complete! Results saved to:")
        print("- 'trna_analysis_results.csv'")
        print("- 'trna_sample_comparison.csv'")
        print("- 'trna_analysis_results.png' (main figure with all plots)")
        print("- 'sample_comparison_plot.png' (Plot 6: sample comparison)")
        print("- 'trna_genes_and_samples_plot.png' (Plot 7: tRNA genes + both samples)")
    else:
        print("Analysis failed.")

if __name__ == "__main__":
    main()