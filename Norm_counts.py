#!/usr/bin/env python

### normalize counts.tsv from RNA-Seq snakepipe ###

import argparse
import pandas as pd
import sys
import os

def parse_args(defaults={
                         "normalization":["RPKM","TPM"],
                         "out_dir":"./"
                         }):

    """parse arguments from the command line"""

    parser = argparse.ArgumentParser(
        prog=sys.argv[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        #description=textwrap.dedent(__description__),
        add_help= False
    )


    # positional - required
    parser.add_argument("counts",
                       metavar="COUNTS",
                       help="counts.tsv output of RNA-Seq snakepipe")

    parser.add_argument("gene_lengths",
                       metavar="LENGTHS",
                       help="feature lengths file")

    # general arguments
    general = parser.add_argument_group('general arguments')
    general.add_argument("-o", "--output_directory",
                         dest="out_dir",
                         help="output directory",
                         default=defaults["out_dir"])

    general.add_argument("-norm","--normalization",
                         dest="NORM",
                         nargs="*",
                         choices=["CPM","RPKM","TPM","TMM","RLE","quantile"],
                         default=defaults["normalization"])


    return parser


def readin_data(counts_tsv, lengths):

    df_counts = pd.read_csv(counts_tsv,sep="\t",index_col=0)
    df_lengths = pd.read_csv(lengths,sep="\t",index_col=0,comment="#")

    return df_counts, df_lengths

def CPM(df_counts):

    return (df_counts / df_counts.sum()) * 1000000

def RPKM(df_counts, df_lengths):

    norm_fact = df_counts.sum() / 1000000

    df_counts_n = df_counts.divide(norm_fact, axis="columns")

    df_counts_rpkm = df_counts_n.divide(df_lengths[list(df_lengths)[0]]/1000, axis="index")

    return df_counts_rpkm

def TPM(df_counts, df_lengths):

    df_counts_tpm = df_counts.divide(df_lengths[list(df_lengths)[0]]/1000, axis="index")

    perMil_sf = df_counts_tpm.sum() / 1000000

    df_counts_tpm = df_counts_tpm.divide(perMil_sf, axis="columns")

    return df_counts_tpm




def main():

    parser = parse_args()
    arg = parser.parse_args()

    df_counts, df_lengths = readin_data(arg.counts, arg.gene_lengths)

    base = ".".join(os.path.basename(arg.counts).split('.')[:-1])
    outpath = os.path.join(arg.out_dir,base)
    #print(outpath)

    for n in arg.NORM:

        if n == "RPKM":
            df_rpkm = RPKM(df_counts, df_lengths)
            df_rpkm.to_csv("{}_RPKM.tsv".format(outpath),sep="\t")

        elif n == "TPM":
            df_tpm = TPM(df_counts, df_lengths)
            df_tpm.to_csv("{}_TPM.tsv".format(outpath),sep="\t")



if __name__ == "__main__":
    main()
