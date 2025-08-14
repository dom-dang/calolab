# Copyright (C) 2024 The Ohio State University

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

# Searches a genome for complementarity to an input sequence(s). Can allow/specify mismatches.
# Searches for perfect complementarity first, then, if selected, for mismatched complementarity. After each input sequence, complementary genes are listed. Genes that fit both perfect and mismatched criteria are listed only once. As the program runs, data for each region of complementarity is listed.
# Sequence(s) can be input manually or the program can be adjusted to allow multiple sequences to run consecutively from a text file [‘Input Sequences’, ‘Sequence Names’]. Can be used to measure complementarity in a list of sequences generated using ‘Random Sequence Generator.’
# Needed for use: properly formatted genome (see: genome formatter) [‘Formatted_Genome’].
# Needed for use properly formatted gene names (see: gene name formatter) [‘Formatted_Gene_Names’].

import multiprocessing
import time


def reverse_complement(inputSequence):
    reverseSeq = inputSequence[::-1].upper().replace('U','T')
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reverseComplement = ''.join(complement[base] for base in reverseSeq)
    return reverseComplement


def substring_match(string1, string2, segment_length, mismatches, gene_index, geneNames):
    # function used for mismatched complementarity
    substrings1 = []
    for i in range(len(string1) - (segment_length - 1)):
        substrings1.append(string1[i:i + segment_length])
    substrings2 = []
    for j in range(len(string2) - (segment_length - 1)):
        substrings2.append(string2[j:j + segment_length])
    # splits the genome and input sequence (rc) into all possible substrings of a specified length.
    for substring1 in substrings1:
        for substring2 in substrings2:
            sum1 = 0
            for x in range(segment_length):
                if substring1[x] == substring2[x]:
                    sum1 += 1
                else:
                    sum1 += 0
                # counts the overlapping nucleotides between each combination of substrings
                if sum1 >= (segment_length - mismatches):
                    a = string2.find(substring2) - string1.find(substring1)
                    print(f'Input: {reverse_complement(substring1)}')
                    print(f'RC Input: {substring1}')
                    print(f'Gene Segment: {substring2}')
                    print(f'Gene Segment Relative to Input: {string2[a:a+len(string1)]}')
                    print(geneNames[gene_index], end='')
                    print()
                    return gene_index
                    # if the overlapping nucleotides exceeds the minimum required, returns the gene index and prints info.

if __name__ == "__main__":
    genome = open('Formatted_Genome', 'r').readlines()
    geneNames = open('Formatted_Gene_Names', 'r').readlines()
    seq = input('Input Sequence: ')
    allInputSeq = [seq]
    # can also use text file input for multiple sequences at once by including the following lines.
    # allInputSeq = open('Input Sequences', 'r').readlines()
    inputSeqNames = ['Sequence 1']
    # inputSeqNames = open('Sequence Names', 'r').readlines()



    segmentLength = int(input('Segment Length? (nts of perfect complementarity): '))
    YNmismatch = input('Allow mismatches Y/N: ').upper()
    if YNmismatch == 'Y':
        segmentLengthWithMismatches = int(input('Segment Length? (with mismatches): '))
        allowedMismatches = int(input('Allowed Mismatches: '))
    geneIndexListGlobal = []
    sequenceComp = []

    seqNumber = 0
    for sequence in allInputSeq:
        print(inputSeqNames[seqNumber])
        seqAdjusted = reverse_complement(sequence)
        # creates the reverse complement of the input sequence.
        # to search for matching sequences, include the following line:
        # seqAdjusted = sequence.upper().replace('U','T')


        geneIndexList = []
        foundGeneNames = []

        n1 = 0
        n2 = segmentLength
        segmentList = []

        while n2 <= len(seqAdjusted):
            segmentList.append(seqAdjusted[n1:n2])
            n1 += 1
            n2 += 1
            # splits input sequence (rc) into all possible segments (substrings) of a given length
        geneIndexPerfect = 0
        for gene in genome:
            exist = False
            for segment in segmentList:
                if segment in gene and exist is False:
                    a = gene.find(segment) - seqAdjusted.find(segment)
                    print(f'Input: {reverse_complement(segment)} *perfect')
                    print(f'RC Input: {segment}')
                    print(f'Gene Segment Relative to Input: {gene[a:a + len(seqAdjusted)]}')
                    exist = True
                    # searches each gene in the genome for each substring of the input sequence (rc).
                    # if a match is found in a gene, prints information about the gene and the paired region.
                    # used for perfect complementarity
                    # note: if there are multiple positive regions in the same gene, this will only detect the first region
            if exist is True:
                geneIndexList.append(geneIndexPerfect)
                # appends indexes of genes with complementarity fitting criteria.
                print(geneNames[geneIndexPerfect])
            geneIndexPerfect += 1

        if YNmismatch == 'Y':
            # ran only if mismatches are selected (computationally intensive).
            print('*perfect complementarity search complete, transitioning to mismatched complementarity')
            print()
            pool = multiprocessing.Pool()
            start_time = time.perf_counter()
            processes = [pool.apply_async(substring_match, args=(seqAdjusted, gene, segmentLengthWithMismatches, allowedMismatches, geneIndex, geneNames)) for gene, geneIndex in zip(genome, range(len(genome)))]
            # runs the mismatched complementarity function in parallel for each gene.
            for z in processes:
                data = z.get()
                if data is not None:
                    geneIndexList.append(data)
                    # appends list positions for genes with complementarity fitting criteria.

        geneIndexList = set(geneIndexList)
        # keeps only one copy of duplicates fitting both perfect and mismatched criteria.

        for gene in geneIndexList:
            foundGeneNames.append(geneNames[gene])
            geneIndexListGlobal.append(geneNames[gene])

        foundGeneNames = set(foundGeneNames)
        # per input sequence, list of genes fitting criteria.
        
        print('List of hits: ')
        for entry in foundGeneNames:
            print(entry, end='')
        print(f'Number of hits: {str(len(foundGeneNames))}')
        print()
        # prints number of ORFs/genes/UTRs/etc. fitting criteria.
        sequenceComp.append(len(foundGeneNames))
        seqNumber += 1

    print('Sequence names: ')
    for name in inputSeqNames:
        print(name)

    print('Corresponding numbers of hits/genes: ')
    for hitCount in sequenceComp:
        print(hitCount)
    # for all input sequences (global), list of genes fitting criteria.

    print()
    print(f'Total unique hits/genes fitting criteria: {len(set(geneIndexListGlobal))}')
    # total number of unique genes fitting criteria among all the input sequences.
