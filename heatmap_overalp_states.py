#!/home/ferrari/anaconda3/bin/python3

# This script computes the overlap between one chromatin segmentation and a second chromatin segmentation  

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("segmentation_file1", help='Segmentation file in bed format \
                    output from ChromHMM LearnModel command')
parser.add_argument("segmentation_file2", help='Segmentation file in bed format \
                    you want to compute the overlap with')
parser.add_argument("-b", "--bin_size", dest="bin_size", type = int, default = 200,
help="binning size of the genome. This parameter should equal the binning size \
        chosen for ChromHMM. default = 200")
parser.add_argument("-o", "--output", dest="output", help="file name", default="heatmap_chromatin_states_overlap.png")
parser.add_argument("-e", "--extension", dest= "extension", help="file extension", default="png", choices=["png","pdf","jpg"])

args = parser.parse_args()

import os
import pandas as pd
import itertools
import matplotlib.pyplot as plt
import seaborn as sns
import math
plt.switch_backend('agg')

def bin_region(lista, bin_size):

    '''
    This function takes as input a list of the format [chr, start, end] and a bin size (int)
    and returns a list of lists with the coordinates of the binned regions
    '''

    out_lista=[]

    if (int(lista[2])-int(lista[1])) % bin_size == 0:
        n_bins = (int(lista[2])-int(lista[1]))/bin_size
        start = int(lista[1])
        end = int(lista[1]) + bin_size
        out_lista.append([lista[0], start, end])
        for i in range(int(n_bins - 1)):
            out_lista.append([lista[0],out_lista[-1][2],out_lista[-1][2]+bin_size])
    else:
        print('Warning! The chosen bin size does not match the segmentation bin size')
        exit(1)

    return out_lista

# get absolute paths of input files
seg_file1 = os.path.abspath(args.segmentation_file1)
seg_file2 = os.path.abspath(args.segmentation_file2)
bin_size = args.bin_size

# usefull data
states_1 = []
states_2 = []

print("Let's get started dude!")

#import segmentation data from file 1
with open(seg_file1) as seg1:
    dixio_1 = dict()
    for line in seg1:
        lista = line.strip().split()
        if lista[3].startswith("E"):
            lista[3] = lista[3][1:]
        if not lista[3] in states_1:
            states_1.append(lista[3])
        for k in bin_region(lista, bin_size):
            dixio_1['\t'.join([str(j) for j in k])] = lista[3]

print('Segmentation file 1 successfully imported. Identified {} bins'.format(len(dixio_1)))

# import segmentation data from file 2
with open(seg_file2) as seg2:
    dixio_2 = dict()
    for line in seg2:
        lista = line.strip().split()
        if lista[3].startswith("E"):
            lista[3] = lista[3][1:]
        if not lista[3] in states_2:
            states_2.append(lista[3])
        for k in bin_region(lista, bin_size):
            dixio_2['\t'.join([str(j) for j in k])] = lista[3]

print('Segmentation file 2 successfully imported. Identified {} bins'.format(len(dixio_2)))

# check if the number of states are equal in both files
if len(states_1) != len(states_2):
    print('Error! Differring number of states in the two segmentation files provided!')
    exit(1)

final_states = ['{}'.format(j) for j in sorted([int(l) for l in states_1])]

print('Chromatin states identified: {}'.format(final_states))

# initialize final dictionary
dixio_def = {}
for i in final_states:
    dixio_def[i]={'count_{}'.format(i):0}
    for k in final_states:
        dixio_def[i][k]=0

# count overlaps
for i in dixio_1:
    A_state=dixio_1[i]
    B_state=dixio_2[i]
    dixio_def[A_state]['count_{}'.format(A_state)]+=1
    dixio_def[A_state][B_state]+=1

# compute percentage of overlap
for k in dixio_def:
    for l in dixio_def[k]:
        if not l.startswith('count_'):
            dixio_def[k][l]=round(((dixio_def[k][l]/dixio_def[k]['count_{}'.format(k)])*100),2)
    del dixio_def[k]['count_{}'.format(k)]


# create dataframe from dictionary
plot_df=pd.DataFrame.from_dict(dixio_def)
# order columns
plot_df = plot_df[final_states]


# order rows
ordered_rows=[]
dup=[]
for i in final_states:
    if not plot_df[i].idxmax() in ordered_rows:
        ordered_rows.append(plot_df[i].idxmax())
    else:
        dup.append(i)

left=[i for i in final_states if not i in ordered_rows]
ordered_rows=ordered_rows+left
plot_df = plot_df.reindex(ordered_rows)


# plot
fig, ax = plt.subplots(figsize=(7,6))

sns.heatmap(plot_df.T, cmap='bone_r', square=True, cbar_kws={'label': "\n% of overlap"}, )
ax.set(ylabel = "{}\n".format(os.path.basename(args.segmentation_file1)) , xlabel="\n{}".format(os.path.basename(args.segmentation_file2)))
plt.savefig(".".join([args.output, args.extension]), dpi=100)
