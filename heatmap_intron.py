import matplotlib.pyplot as plt
import numpy as np
import re
import seaborn as sns
import argparse
from matplotlib.colors import LinearSegmentedColormap
from Bio import SeqIO
from Bio.Seq import Seq

def parse_input_file(file_path):
    """
    Parse the input file containing sequences and their descriptions.
    Returns a list of dictionaries with sequence information.
    """
    sequences = []
    current_seq = None
    
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            
            # Skip empty lines
            if not line:
                continue
                
            # Handle FASTA header lines
            if line.startswith('>'):
                continue
                
            # Parse sequence lines with labels
            if line.startswith('Input:'):
                if current_seq:
                    sequences.append(current_seq)
                current_seq = {'name': 'Input', 'sequence': line.split('Input:')[1].strip()}
            elif line.startswith('RC Input:'):
                if current_seq:
                    current_seq['rc_sequence'] = line.split('RC Input:')[1].strip()
            elif line.startswith('Gene Segment:'):
                if current_seq:
                    current_seq['gene_segment'] = line.split('Gene Segment:')[1].strip()
            elif line.startswith('Gene Segment Relative to Input:'):
                if current_seq:
                    current_seq['gene_segment_rel'] = line.split('Gene Segment Relative to Input:')[1].strip()
            # If it's not a labeled line and we have a current sequence, it might be a continuation
            elif current_seq and not line.startswith('>'):
                # Handle continuation lines if needed
                pass
    
    # Add the final sequence if there is one
    if current_seq:
        sequences.append(current_seq)
        
    return sequences

def find_all_occurrences(main_str, sub_str):
    """
    Find all positions where sub_str appears in main_str (case-insensitive).
    """
    main_str_lower = main_str.lower()
    sub_str_lower = sub_str.lower()
    positions = []
    pos = main_str_lower.find(sub_str_lower)
    
    while pos != -1:
        positions.append(pos)
        pos = main_str_lower.find(sub_str_lower, pos + 1)
    
    return positions

def create_alignment_matrix(intron, sequences):
    """
    Create a matrix representing all alignments.
    1 indicates a match between the sequence and the intron at that position.
    """
    # Initialize matrix with zeros
    matrix = np.zeros((len(sequences) * 4, len(intron)))  # 4 types per sequence (Input, RC, Gene, Gene Rel)
    labels = []
    
    row_index = 0
    for seq_info in sequences:
        # Process each type of sequence for this entry
        for seq_type in ['sequence', 'rc_sequence', 'gene_segment', 'gene_segment_rel']:
            if seq_type in seq_info:
                seq = seq_info[seq_type]
                # Find all occurrences of the sequence in the intron
                positions = find_all_occurrences(intron, seq)
                
                # If we found matches, mark them in the matrix
                if positions:
                    for pos in positions:
                        for i in range(len(seq)):
                            if pos + i < len(intron):
                                matrix[row_index, pos + i] = 1
                
                # Create label for this row
                type_name = {
                    'sequence': 'Input',
                    'rc_sequence': 'RC Input',
                    'gene_segment': 'Gene Segment',
                    'gene_segment_rel': 'Gene Segment Rel'
                }[seq_type]
                
                labels.append(type_name + seq)
                row_index += 1
    
    # Remove empty rows
    non_empty_rows = np.any(matrix, axis=1)
    matrix = matrix[non_empty_rows]
    labels = [label for i, label in enumerate(labels) if non_empty_rows[i]]
    
    return matrix, labels

def plot_heatmap(intron, matrix, labels, output_file):
    """
    Create and save a heatmap visualization of the alignments.
    """
    plt.figure(figsize=(20, max(8, len(labels) * 0.4)))
    
    colors = [(1, 1, 1), (0, 0.5, 1)]  # White to blue
    cmap = LinearSegmentedColormap.from_list('custom_cmap', colors, N=2)
    
    # Plot the heatmap
    ax = sns.heatmap(matrix, cmap=cmap, cbar=False, linewidths=0.5)
    
    plt.yticks(np.arange(len(labels)) + 0.5, labels, fontsize=8)
    
    # Set x-axis labels (the intron sequence)
    plt.xticks(np.arange(len(intron)) + 0.5, list(intron.upper()), fontsize=8)
    
    # Add position numbers below the intron sequence
    # ax2 = ax.twiny()
    # ax2.set_xticks(np.arange(len(intron)) + 0.5)
    # ax2.set_xticklabels([str(i+1) for i in range(len(intron))], fontsize=6)
    # ax2.tick_params(axis='x', which='both', length=0)
    
    plt.title('Intron Alignment Heatmap', fontsize=16)
    plt.tight_layout()
    
    # Save the figure
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print("Heatmap saved to" + output_file)
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='Create alignment heatmap for intron sequences')
    parser.add_argument('intron', help='Intron sequence')
    parser.add_argument('input_file', help='Path to the input file containing sequences')
    parser.add_argument('--output', '-o', default='intron_alignment_heatmap.png', help='Output file name')
    
    args = parser.parse_args()
    
    sequences = parse_input_file(args.input_file)

    matrix, labels = create_alignment_matrix(args.intron, sequences)
    
    if matrix.size == 0:
        print("No alignments found between the sequences and the intron.")
        return

    plot_heatmap(args.intron, matrix, labels, args.output)
    
    print("Found {} alignments across {} input sequences.".format(len(labels), len(sequences)))
    print("Intron length: {} bases".format(len(args.intron)))


if __name__ == "__main__":
    main()