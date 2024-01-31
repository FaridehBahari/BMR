import pybedtools
import pandas as pd

# Assume df1 and df2 are your pre-existing DataFrames and have a 'name' column corresponding to genomic intervals.
# Let's say bed_file1 and bed_file2 are the paths to the BED files for df1 and df2.

# Read the BED files
bed1 = pybedtools.BedTool(bed_file1)
bed2 = pybedtools.BedTool(bed_file2)

# Perform an intersection
intersection = bed1.intersect(bed2, wa=True, u=True)  # wa: Write the original entry in A for each overlap. u: Write original A entry once if any overlaps found.

# Create a set of names from df2 that intersect with df1
intersecting_names = set([feature.name for feature in intersection])

# Label df1: 1 if the row's name is in the intersection set; otherwise 0
df1['label'] = df1['name'].apply(lambda x: 1 if x in intersecting_names else 0)

# Now, df1 will have a new column labeled 'label' with 1s where names intersect with those in df2 based on the BED files, and 0s otherwise.
