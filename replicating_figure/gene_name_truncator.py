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

# Truncates formatted gene names from annotated genomes to only include the area of interest in each annotation.

with open('Formatted_Gene_Names', 'r') as d:
    e = d.readlines()
with open('Truncated_Gene_Names', 'w') as f:
    for line in e:
        g = line.split()
        f.write(g[1] + '\n')
        # writes the second word in the annotation, for S. cerevisiae, this corresponds to the gene name.
        # may be different for other organisms
        # depending on the annotations, it may be better to use a specific region of the annotation instead:
        # f.write(line[length1:length2] + '\n')

