# -*- coding: utf-8 -*-
import multiprocessing
import time
import re
import glob

def reverse_complement(inputSequence):
    reverseSeq = inputSequence[::-1].upper().replace('U','T')
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reverseComplement = ''.join(complement[base] for base in reverseSeq)
    return reverseComplement

def substring_match(string1, string2, segment_length, mismatches, gene_index, geneNames):
    substrings1 = [string1[i:i + segment_length] for i in range(len(string1) - (segment_length - 1))]
    substrings2 = [string2[j:j + segment_length] for j in range(len(string2) - (segment_length - 1))]

    for substring1 in substrings1:
        for substring2 in substrings2:
            sum1 = sum(1 for x in range(segment_length) if substring1[x] == substring2[x])
            if sum1 >= (segment_length - mismatches):
                a = string2.find(substring2) - string1.find(substring1)
                print(f'Input: {reverse_complement(substring1)}')
                print(f'RC Input: {substring1}')
                print(f'Gene Segment: {substring2}')
                print(f'Gene Segment Relative to Input: {string2[a:a+len(string1)]}')
                print(geneNames[gene_index], end='')
                print()
                return gene_index

if __name__ == "__main__":
    genome = open('Formatted_Genome', 'r').readlines()
    geneNames = open('Formatted_Gene_Names', 'r').readlines()

    allInputSeq = []
    inputSeqNames = []
    sourceFiles = []

    for file in glob.glob("*_sequence_name.txt"):
        with open(file, 'r') as f:
            for line in f:
                match = re.match(r'\[(.+?),\s*(.+?)\]', line.strip())
                if match:
                    name, seq = match.groups()
                    inputSeqNames.append(name.strip())
                    allInputSeq.append(seq.strip())
                    sourceFiles.append(file)

    print(f"Loaded {len(allInputSeq)} sequences from {len(set(sourceFiles))} files.\n")

    segmentLength = int(input('Segment Length? (nts of perfect complementarity): '))
    YNmismatch = input('Allow mismatches Y/N: ').upper()
    if YNmismatch == 'Y':
        segmentLengthWithMismatches = int(input('Segment Length? (with mismatches): '))
        allowedMismatches = int(input('Allowed Mismatches: '))

    geneIndexListGlobal = []
    sequenceComp = []

    for seqNumber, sequence in enumerate(allInputSeq):
        print(f'\nAnalyzing: {inputSeqNames[seqNumber]} (from file: {sourceFiles[seqNumber]})')
        seqAdjusted = reverse_complement(sequence)

        geneIndexList = []
        foundGeneNames = []

        segmentList = [seqAdjusted[i:i+segmentLength] for i in range(len(seqAdjusted) - segmentLength + 1)]

        for geneIndexPerfect, gene in enumerate(genome):
            exist = False
            for segment in segmentList:
                if segment in gene and not exist:
                    a = gene.find(segment) - seqAdjusted.find(segment)
                    print(f'Input: {reverse_complement(segment)} *perfect')
                    print(f'RC Input: {segment}')
                    print(f'Gene Segment Relative to Input: {gene[a:a + len(seqAdjusted)]}')
                    exist = True
            if exist:
                geneIndexList.append(geneIndexPerfect)
                print(geneNames[geneIndexPerfect])

        if YNmismatch == 'Y':
            print('*perfect complementarity search complete, transitioning to mismatched complementarity\n')
            pool = multiprocessing.Pool()
            processes = [pool.apply_async(substring_match, args=(seqAdjusted, gene, segmentLengthWithMismatches, allowedMismatches, idx, geneNames))
                         for idx, gene in enumerate(genome)]
            for z in processes:
                data = z.get()
                if data is not None:
                    geneIndexList.append(data)

        geneIndexList = set(geneIndexList)

        for gene in geneIndexList:
            foundGeneNames.append(geneNames[gene])
            geneIndexListGlobal.append(geneNames[gene])

        foundGeneNames = set(foundGeneNames)

        print('List of hits:')
        for entry in foundGeneNames:
            print(entry, end='')
        print(f'\nNumber of hits: {len(foundGeneNames)}\n')
        sequenceComp.append(len(foundGeneNames))

    print('\nSummary:')
    print('Sequence names and source files:')
    for name, file in zip(inputSeqNames, sourceFiles):
        print(f'{name} — from file: {file}')

    print('\nCorresponding numbers of hits/genes:')
    for hitCount in sequenceComp:
        print(hitCount)

    print(f'\nTotal unique hits/genes fitting criteria: {len(set(geneIndexListGlobal))}')

