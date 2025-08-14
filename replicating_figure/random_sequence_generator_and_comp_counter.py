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

# Randomly generates and prints sequences with a specified GC content and searches for perfect complementarity (no mismatches).
# Needed for use: properly formatted genome (see: genome formatter) [‘Formatted_Genome’]

from random import choice
from random import shuffle

def reverse_complement(inputSequence):
    reverseSeq = inputSequence[::-1].upper().replace('U','T')
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reverseComplement = ''.join(complement[base] for base in reverseSeq)
    return reverseComplement


def sequence_generator():
    # generates random sequences.
    sequence = [choice('GC') for b in range(round(GCcontent * seqLength))]
    # chooses 'G' or 'C' a number of times dependent on the specified GC content, and adds it to a list.
    sequence += [choice('TA') for c in range(seqLength - round(GCcontent * seqLength))]
    # chooses 'T' or 'A' a number of times dependent on the specified GC content, and adds it to a list.
    shuffle(sequence)
    shuffle(sequence)
    shuffle(sequence)
    generatedSequence = ''.join(sequence)
    # turns list of nucleotides into a string.
    return generatedSequence

seqLength = int(input('Desired Sequence Length: '))
GCcontent = float(input('GC Ratio (decimal): '))
loop = int(input('Number of Sequences to Generate: '))

generatedSequenceList = []

for _ in range(loop):
    generatedSequenceList.append(sequence_generator())

segmentLength = int(input('Desired Segment Length (nts perfect complementarity): '))

geneCountList = []

# Include the line below to read from a text file, instead of the newly generated sequences.
# This allows for complementarity and matching to be run on the same sequences.
# generatedSequenceList = open('Sequence_List','r').readlines()
for seq in generatedSequenceList:
    seqAdjusted = reverse_complement(seq)
    # get the reverse complement ssDNA of the random sequence(s).
    # for matching rather than complementarity, include the line below:
    # seqAdjusted = seq.upper().replace('U', 'T')
    n1 = 0
    n2 = segmentLength
    segmentList = []

    while n2 <= len(seqAdjusted):
        segmentList.append(seqAdjusted[n1:n2])
        n1 += 1
        n2 += 1
    # splits reverse complement ssDNA into all possible segments of a desired length

    geneCount = 0

    genome = open('Formatted_Genome', 'r').readlines()
    for gene in genome:
        exist = False
        for segment in segmentList:
            if segment in gene:
                exist = True
    # searches each line (gene/ORF/UTR/etc.) of the genome for each segment of the reverse complement ssDNA.
        if exist is True:
            geneCount += 1
    # tallies number of lines with positive hits.

    print(seq)
    geneCountList.append(geneCount)
    print('Gene #: ' + str(geneCount))
    # prints results and sequences.

print()
print('Output: ')
for sequence in generatedSequenceList:
    print(sequence)
for count in geneCountList:
    print(count)
    # prints results in a way that can be copied into a spreadsheet.
