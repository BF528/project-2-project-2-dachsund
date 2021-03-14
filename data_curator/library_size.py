#!/usr/bin/env python3
"""
2021/03/14
Allison Nau

Python script to calculate average library size for the 8 samples looked at.

Will print the minimum, maximum, and average library size.

Expects info files to be in a subdirectory "prep_reads_info" in the current working directory.
Library size of P0_1 is loaded in by default WITHOUT an info file
"""

# Import packages:
import os
from statistics import mean

# Initialize list of library sizes with P0_1 already in:
library_sizes = [21577562]

# Get current working directory and append subfolder:
my_dir = os.getcwd()
# Raw string to deal with backslash:
my_dir += r"\prep_reads_info"
print("Expected directory containing info files:")
print(my_dir)

# Walk through info files:
for (root, dirs, files) in os.walk(my_dir):
    for my_file in files:
        open_me = my_dir + "\\" + my_file
        with open(open_me, "r") as the_file:
            lines = the_file.readlines()
            for line in lines:
                line = line.split("=")
                field, value = line[0], line[1]
                if field[0:13] == "left_reads_in":
                    value = int(value.strip())
                    library_sizes.append(value)

print("All library sizes:")
print(library_sizes)
print("Minimum: ", min(library_sizes))
print("Maximum: ", max(library_sizes))
print("Average: ", mean(library_sizes))

