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

# Creates a file with each gene name on a separate line. Each line on the formatted gene name file will correspond to the formatted genome file.
# Needed for use: annotated genome [‘Original_Genome’].
# Needed for use: blank file for formatted gene names to be written [‘Formatted_Gene_Names’].

with open('Original_Genome', 'r') as a:
    lines = a.readlines()
with open("Formatted_Gene_Names", 'w') as b:
    for line in lines:
        if line.startswith(">"):
            b.write(line.replace('\n', ''))
        else:
            b.write("\n")
        # deletes everything except annotations (any line not starting with '>').

with open("Formatted_Gene_Names", 'r+') as c:
    allLines = c.readlines()
    c.truncate(0)
    c.seek(0)
    blankLinesExcluded = [line for line in allLines if line.strip()]
    c.writelines(blankLinesExcluded)
    # removes empty/blank lines.
