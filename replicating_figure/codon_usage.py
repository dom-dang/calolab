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

# Determines ATA and ATG codon usage from a list of gene names.
# Needed for use: properly formatted genome (see: genome formatter) [‘Formatted_Genome’].
# Needed for use properly formatted gene names (see: gene name formatter) [‘Formatted_Gene_Names’].

def codon_usage(gene, allGeneNames, genome):
    geneIndex = allGeneNames.index(gene)
    # finds position of input gene in genome
    geneSequence = genome[geneIndex]
    # finds sequence of input gene
    totalCodons = 0
    ATAcount = 0
    ATGcount = 0

    for i in range(0, len(geneSequence), 3):
        codon = geneSequence[i:i+3]
        # seperates gene sequence into codons
        if codon == 'ATA':
            ATAcount+= 1
        # counts ATA codons in gene sequence
        if codon == 'ATG':
            ATGcount += 1
        # counts ATG codons in gene sequence
        totalCodons += 1

    ATAusage = ATAcount/totalCodons*100
    ATGusage = ATGcount/totalCodons*100
    return ATAusage, ATGusage

inputGeneName = input('Gene: ') + '\n'
allGeneNames = open('Formatted_Gene_Names','r').readlines()
genome = open('Formatted_Genome','r').readlines()
if inputGeneName in allGeneNames:
    ATA, ATG = codon_usage(inputGeneName, allGeneNames, genome)
    print(f'ATA codon usage: {round(ATA, 2)}%')
    print(f'ATG codon usage: {round(ATG, 2)}%')
else:
    print('Error: input gene name does not correspond to a name in the gene names file')
