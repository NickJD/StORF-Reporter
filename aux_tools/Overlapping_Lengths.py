import gzip
import glob
import sys
import collections
import textwrap
import argparse
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FuncFormatter




def overlap_length(range1, range2):

    overlap_start = max(range1[0], range2[0])
    overlap_end = min(range1[1], range2[1])
    overlap_length = max(0, overlap_end - overlap_start)
    return overlap_length



overlap_lengths = collections.defaultdict(int)

sequences = []
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Overlapping CDS Gene Lengths.')
    parser.add_argument('-d', action="store", dest='dir_in', required=True,
                        help='Input Directory with file extension separated by comma (./dir_of_cds_files,.fa.gz)')
    parser.add_argument('-o', action="store", dest='out_file', required=True,
                        help='Default - Without filetype - default appends \'_StORF-R\' to end of input gff filename (replaces \'.gff\')')

    options = parser.parse_args()

    directory_in = options.dir_in.split(',')[0]
    extension_in = '.49.gff3'#options.dir_in.split(',')[1]

    gff_list = list(glob.glob(directory_in+'/*.gff*'))
    counter = 0

    # Initialize variables to keep track of overlapping regions
    prev_gene = None
    prev_end = None
    overlap_start = None
    overlap_end = None

    for gff_file in gff_list:
        current_identifier = gff_file.split('/')[-1].split(extension_in)[0]
        try:
            gff_in = gzip.open(gff_file, 'rt')
            gff_in.read(1)
        except gzip.BadGzipFile:
            gff_in = open(gff_file, 'r')
        gff_in.seek(0)
        for line in gff_in:
            # Skip comment lines
            if line.startswith('#'):
                continue
            # Split the line into fields
            fields = line.strip().split('\t')

            feature = fields[2]
            start = int(fields[3])
            end = int(fields[4])
            # If this is a gene feature, check for overlap with the previous gene
            if feature == 'gene':
                # If there was a previous gene and it overlaps with this gene
                if prev_gene is not None:

                    range1 = (start, end)
                    range2 = (prev_gene['start'], prev_gene['end'])
                    overlap_len = overlap_length(range1, range2)
                    if overlap_len != 0:
                        overlap_lengths[overlap_len] +=1

                # Store the current gene as the previous gene for the next iteration
                prev_gene = {
                    'start': start,
                    'end': end,
                }
                prev_end = end
        # Output the lengths of the overlapping regions
        for key, value in overlap_lengths.items():
            print(key)
            print(value)






# Calculate the total count of all values
total_count = sum(overlap_lengths.values())

# Calculate the proportion of each value
proportions = {k: v / total_count for k, v in overlap_lengths.items()}

# Calculate the cumulative proportion of each value
cumulative_proportions = {}
cumulative_sum = 0
for k, v in proportions.items():
    cumulative_sum += v
    cumulative_proportions[k] = cumulative_sum

# Create a bar chart of the proportions
fig, ax = plt.subplots()
ax.bar(proportions.keys(), proportions.values(), label='Proportions')

# Create a line chart of the cumulative proportions
ax.plot(cumulative_proportions.keys(), cumulative_proportions.values(), label='Cumulative Proportions')

# Set the chart title and axis labels
ax.set_title('Proportions of Values and Cumulative Proportions')
ax.set_xlabel('Value')
ax.set_ylabel('Proportion')
plt.xlim([0,50])

# Add a legend to the chart
ax.legend()

# Display the chart
plt.show()















