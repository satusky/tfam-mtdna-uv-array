# Code to produce analysis similar to panels B and C of Figure 3 figure supplement 6.
#
# Author: Evan Corden
#
import numpy as np
import re
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Read in dataset
all_stats = pd.read_csv("data/tfam_array_all_stats.csv", index_col=[0,1], skipinitialspace=True, low_memory=False)

# Keep only the genomic probe data
genome = all_stats[all_stats["start_coord"] > 0].copy(deep=True)
genome = genome[[(not str(i[0]).lower().endswith("_a")) and (not str(i[0]).lower().endswith("_g")) and (not str(i[0]).lower().endswith("_t")) for i in genome.index]]
genes = [str(i[0]).split("_")[0] for i in genome.index]

# Separate forward (_f) and reverse (_r) probes
genome["gene"] = genes
f_genome = genome[["_o1" in str(i[0]) for i in genome.index]]
r_genome = genome[["_o2" in str(i[0]) for i in genome.index]]

# Separate by control vs UV and TFAM level
c30_f = f_genome.xs("ctrl_30", level=1).sort_values(by="mid_coord")
uv30_f = f_genome.xs("uv_30", level=1).sort_values(by="mid_coord")
c300_f = f_genome.xs("ctrl_300", level=1).sort_values(by="mid_coord")
uv300_f = f_genome.xs("uv_300", level=1).sort_values(by="mid_coord")

c30_r = r_genome.xs("ctrl_30", level=1).sort_values(by="mid_coord")
uv30_r = r_genome.xs("uv_30", level=1).sort_values(by="mid_coord")
c300_r = r_genome.xs("ctrl_300", level=1).sort_values(by="mid_coord")
uv300_r = r_genome.xs("uv_300", level=1).sort_values(by="mid_coord")

# We use control 30 TFAM level for this analysis
# Make a copy for safe editing
c30_r_copy = c30_r.copy()

# Reset index to remove orientations for later merging
c30_r_copy.index =  c30_r_copy.index.map(lambda x: re.sub("_o2", "", x))

c30_f_copy = c30_f.copy()
c30_f_copy.index =  c30_f_copy.index.map(lambda x: re.sub("_o1", "", x))

# Join on names
# _r suffix means reverse direction, and no suffix is _f
c30_joined = c30_f_copy.join(c30_r_copy, how='inner', rsuffix = "_r")

# Find maximum of forward and reverse probes for same name
c30_joined["max_norm_median"] = np.zeros((len(c30_joined),1))
c30_joined["max_sequence"] = c30_joined["sequence"] # Max sequence is the sequence of the f or r probe with maximum fluorescence
c30_joined["max_dir"] = c30_joined["sequence"] # Direction of probe with max fluorescence

for i, row in c30_joined.iterrows():
    if row["norm_median"] > row["norm_median_r"]:
        c30_joined.at[i, "max_norm_median"] = row["norm_median"]
        c30_joined.at[i,"max_sequence"] = row["sequence"][:33]
        c30_joined.at[i,"max_dir"] = "f"
    else:
        c30_joined.at[i,"max_norm_median"] = row["norm_median_r"]
        c30_joined.at[i,"max_sequence"] = row["sequence_r"][:33]
        c30_joined.at[i,"max_dir"] = "r"

# The bed file for the dgfs use coordinates starting at zero and the PBM coordinates
# start at 1. I am only using start and ends coords in this analysis, so I
# subtract 1 from both data series to match the bed file coordinate system.
c30_joined["start_coord"] += -1
c30_joined["end_coord"] += -1

# Read in coordinates of 29 common regions in the Blumberg paper
dgf_coords = [] # Coordinates of dgf binding sites from Blumberg

with open("data/footprinting_comparison/common_dgf_hg38.bedGraph", "r") as in_file:
    for line in in_file:
        line_array = line.strip().split("\t")
        co1 = int(line_array[1])
        co2 = int(line_array[2])
        
        dgf_coords.append([co1, co2])

# Iterate over c30 and designate values as in or out of TFAM intervals
in_vals_f = [] # Adjusted forward Alexa values for probes entirely overlapping bedGraph intervals
out_vals_f = [] # Adjusted forward Alexa values for probes not overlapping bedGraph intervals
in_vals_r = [] # Adjusted reverse Alexa values for probes entirely overlapping bedGraph intervals
out_vals_r = [] # Adjusted reverse Alexa values for probes not overlapping bedGraph intervals
in_vals_max = [] # Same data structure for max between f and r
out_vals_max = []
for index, row in c30_joined.iterrows():
    inb = False
    for interval in dgf_coords:
        # Normal case- can disregard circular case in both the probes and the bedGraph intervals because we know there are no circular bedGraph intervals that span the end and beginning of the mitochondrial genome coordinates, so it is not possible for any circular probes to be entirely contained in a bedGraph interval.
        if interval[0] < interval[1]:
            if (row['start_coord'] >= interval[0] and row['end_coord'] <= interval[1]) or (interval[0] >= row['start_coord'] and interval[1] <= row['end_coord']) and row['start_coord'] < row['end_coord']:
                inb = True
                break
            
    if inb:
        in_vals_f.append(row['norm_median'])
        in_vals_r.append(row['norm_median_r'])
        in_vals_max.append(row['max_norm_median'])
    else:
        out_vals_f.append(row['norm_median'])
        out_vals_r.append(row['norm_median_r'])
        out_vals_max.append(row['max_norm_median'])

# Some checks
print(len(in_vals_f))
print(len(out_vals_f))
print(in_vals_f[:5])
        
print(len(in_vals_r))
print(len(out_vals_r))
print(in_vals_r[:5])

print(len(in_vals_max))
print(len(out_vals_max))
print(in_vals_max[:5])

# Create null distribution with blocks similar to those in the DGF regions
import random

# Get length distribution of dgf coordinates
length_dist = []
for i in dgf_coords:
    length_dist.append(i[1] + 1 - i[0])

# Find minimum separation between dgf regions
min_sep = 2000
for i in range(len(dgf_coords)):
    if i == 0:
        prev_co = dgf_coords[i][1]
        continue
    else:
        diff = dgf_coords[i][0] - prev_co
        prev_co = dgf_coords[i][1]
        
        # Update min if necessary
        if diff < min_sep:
            min_sep = diff

# handle circular case
diff = (dgf_coords[0][0]) + (16569 - prev_co)
if diff < min_sep:
            min_sep = diff

# Run simulations
num_trials = 1000 # Number of null trials to run
null_means_max = [] # Contains the null distribution of mean values for max fluorescence
null_means_f = [] # Same for forward 
null_means_r = [] 
for j in range(num_trials):
    # Construct simulated data set
    sample = [] 
    for i in range(len(length_dist)):
        
        # Need region that is not overlapping other regions and not within min separation
        good_pick = False
        while good_pick == False:
            good_pick = True
            region_len = length_dist[i] # Length of this simulated region

            index1 = random.randint(0,16569 - 1)
            index2 = index1 + region_len - 1
            
            # Handle circular case
            if index2 > 16568:
                index2 = index2 - 16568 - 1
            
            # Check for overlap by looking over regions already in sample
            for t in sample:
                
                 # Index and region are circular- not possible to not overlap
                if t[0] > t[1] and index1 > index2:
                    good_pick = False
                    break
                    
                # Both non-circular
                elif t[0] < t[1] and index1 < index2:
                    # Overlap cases
                    if index1 >= t[0] - min_sep and index1 <= t[1] + min_sep:
                        good_pick = False
                        break
                    elif index2 >= t[0] - min_sep and index2 <= t[1] + min_sep:
                        good_pick = False
                        break

                    # Existing region is fully inside of new region
                    elif index1 <= t[0] and index2 >= t[1]:
                        good_pick = False
                        break
                
                # Region is non-circular and index is circular
                elif t[0] < t[1] and index1 > index2:
                    # Overlap cases
                    if index1 >= t[0] and index1 <= t[1] + min_sep:
                        good_pick = False
                        break
                        
                    elif index2 >= t[0] - min_sep and index2 <= t[1]:
                        good_pick = False
                        break

                    # Existing region is fully inside of new region
                    elif (index1 <= t[0] and t[1] <= 16568) or (t[0] >= 0 and t[1] <= index2):
                        good_pick = False
                        break
                
                # Region is circular and index is non-circular
                elif t[0] > t[1] and index1 < index2:
                    # Overlap cases
                    if (index1 >= t[0] - min_sep and index1 <= 16568) or (index1 >= 0 and index1 <= t[1] + min_sep):
                        good_pick = False
                        break
                        
                    elif (index2 >= t[0] - min_sep and index2 <= 16568):
                        good_pick = False
                        break

                    # Existing region is fully inside of new region- not possible

        sample.append([index1, index2])

    # Now find probes that fully overlap these regions
    in_vals_max_null = [] # Adjusted Alexa values for max probes entirely overlapping bedGraph intervals
    out_vals_max_null = [] # Adjusted Alexa values for max probes not overlapping bedGraph intervals
    in_vals_f_null = [] # Adjusted Alexa values for f probes entirely overlapping bedGraph intervals
    out_vals_f_null = [] # Adjusted Alexa values for f probes not overlapping bedGraph intervals
    in_vals_r_null = [] # Adjusted Alexa values for r probes entirely overlapping bedGraph intervals
    out_vals_r_null = [] # Adjusted Alexa values for r probes not overlapping bedGraph intervals

    for index, row in c30_joined.iterrows():
        inb = False
        for interval in sample:
            # Normal interval, normal probe
            if interval[0] < interval[1] and row['start_coord'] < row['end_coord']:
                if (row['start_coord'] >= interval[0] and row['end_coord'] <= interval[1]) or (interval[0] >= row['start_coord'] and interval[1] <= row['end_coord']):
                    inb = True
                    break
            
            # Circular interval, circular probe
            elif interval[1] < interval[0] and row['start_coord'] > row['end_coord']:
                # Probe inside interval
                if row['start_coord'] >= interval[0] and row['end_coord'] <= interval[1]:
                    inb = True
                    break
                    
                # Interval inside probe- not possible in this dataset because no interval is less than 33 nt
                if interval[0] >= row['start_coord'] and interval[1] <= row['end_coord']:
                    inb = True
                    break
            
            # Circular interval, normal probe
            elif interval[1] < interval[0] and row['start_coord'] < row['end_coord']:
                # Probe inside of interval on high end
                if row['start_coord'] >= interval[0] and row['start_coord'] <= 16568:
                    inb = True
                    break
                    
                if row['end_coord'] >= 0 and row['end_coord'] <= interval[1]:
                    inb = True
                    break
                    
            # Normal interval, circular probe is not possible if the probe must be fully inside of the interval

        if inb:
            in_vals_max_null.append(row['max_norm_median'])
            in_vals_f_null.append(row['norm_median'])
            in_vals_r_null.append(row['norm_median_r'])
        else:
            out_vals_max_null.append(row['max_norm_median'])
            out_vals_f_null.append(row['norm_median'])
            out_vals_r_null.append(row['norm_median_r'])
            
    null_means_max.append(np.mean(in_vals_max_null))
    null_means_f.append(np.mean(in_vals_f_null))
    null_means_r.append(np.mean(in_vals_r_null))

# In this analysis we are prioritizing the null distribution of the max of the two orientations
print(np.percentile(null_means_max, [90,95,99]))
print(np.mean(in_vals_max))

# Save data to files
#with open("null_means_max.txt", "w") as out_file:
#    for x in null_means_max1:
#        out_file.write(str(x) + "\n")
#
#with open("in_probes_vals.txt", "w") as out_file:
#    for x in in_vals_max:
#        out_file.write(str(x) + "\n")

# Make a boxplot
plt.boxplot([in_vals_max, out_vals_max], labels = ["In", "Out"])
plt.title("Footprints Mean Max Fluorescence")
plt.show()
#plt.savefig("figure1.png")

# Print some additional statistics
print(np.mean(in_vals_max))
print(np.mean(out_vals_max))
print(np.mean(in_vals_max) - np.mean(out_vals_max))

# Make main histogram plot
sns.histplot(null_means_max)
plt.axvline(np.mean(in_vals_max),color='r')
plt.xlabel("Mean Sample Fluorescence")
plt.title("Maximum Median Adjusted Alexa488 Value")
plt.legend(["Mean of In Probes", "Null Distribution"])
plt.show()
#plt.savefig("figure2.png")

