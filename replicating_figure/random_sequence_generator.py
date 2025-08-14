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

# Generates random nucleotide sequences at a specific GC content for use in other programs.

from random import choice
from random import shuffle

seqLength = int(input('Desired Sequence Length: '))
GCcontent = float(input('GC Ratio (decimal): '))
Loop = int(input('Number of Sequences to Generate: '))

for a in range(Loop):
    sequence = [choice('GC') for b in range(round(GCcontent * seqLength))]
    # chooses 'G' or 'C' a number of times dependent on the specified GC content, and adds it to a list.
    sequence += [choice('TA') for c in range(seqLength - round(GCcontent * seqLength))]
    # chooses 'T' or 'A' a number of times dependent on the specified GC content, and adds it to a list.
    shuffle(sequence)
    shuffle(sequence)
    shuffle(sequence)
    # shuffles the list.
    generatedSequence = ''.join(sequence)
    # turns list into string.
    print(generatedSequence)

