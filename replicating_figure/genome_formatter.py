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

# Creates a file with each gene sequence (or other annotated sequence) on a separate line.
# Needed for use: annotated genome [‘Original_Genome’].
# Needed for use: blank file for formatted gene sequences to be written [‘Formatted_Genome’].

with open("Original_Genome", 'r') as a:
    lines = a.readlines()
with open("Formatted_Genome", 'w') as b:
    for line in lines:
        if line.startswith(">"):
            b.write("\n")
        else:
            b.write(line.replace('\n', ''))
        # deletes annotations (any line starting with '>').

with open("Formatted_Genome", 'r+') as c:
    allLines = c.readlines()
    c.truncate(0)
    c.seek(0)
    blankLinesExcluded = [line for line in allLines if line.strip()]
    c.writelines(blankLinesExcluded)
    # removes empty/blank lines.
