import pandas as pd
import time
from functools import reduce
import argparse


def get_args():
    parser = argparse.ArgumentParser("cfDNA Feature Extraction\n")
    parser.add_argument("-d", "--directory",
                        help="Directory with BED files", required=True)
    parser.add_argument("-i", "--input", nargs='+',
                        help="Sample names", required=True)
    return parser


# args = get_args().parse_args()
# workdir = args.directory
# samples = args.input

# workdir = "/cluster/dataset/medinfmk/cfDNA-radiooncology/output/newly_demultiplexed-20210407/BED/"
# samples = ["HN003-d8", "OMD001-3m"]


# @profile
def stratify_reads(workdir, sample, chrs):
    print(">>>>>>> Sample: " + sample)
    start_time = time.time()

    input_f = workdir + sample + ".bed"
    outlier_f = workdir + sample + "_long_short.csv"
    right_f = workdir + sample + "_point.csv"
    right_f_len = workdir + sample + "_len.csv"

    unmatch_c = 0
    sv_c = 0
    long_c = 0
    short_c = 0
    chrX_c = 0
    chrY_c = 0

    with open(input_f) as f1, open(right_f_len, 'w') as f2, open(right_f, 'w'
                                                                 ) as f3, open(outlier_f, 'w') as f4:
        for line in f1:
            fields = line.split('\t')
            # filtering out invalid items
            if '-1' in fields[1:6]:
                unmatch_c += 1
                continue
            # filtering out structure variance
            if fields[0] != fields[3]:
                sv_c += 1
                continue
            # filtering out chrX, chrY:
            if fields[0] not in chrs:
                if fields[0] == 'X':
                    chrX_c += 1
                elif fields[0] == 'Y':
                    chrY_c += 1
                continue
            # compute the fragment length
            l = int(fields[5]) - int(fields[1]) + 1
            fields.append(str(l))

            # filtering out length outliers:
            if l >= 1000:
                long_c += 1
            elif l < 76:
                short_c += 1
            else:
                # write to a new file
                f2.write(fields[10] + '\n')
                f3.write('\t'.join([fields[0][3:], fields[1], fields[10]]) + '\n')
                continue
            f4.write(fields[10] + '\n')
    time1 = time.time()
    print("--- Write to new files: %2f seconds ---" % (time1 - start_time))

    d_nonSV_right = pd.read_csv(right_f_len,
                                sep="\t",
                                names=["length"],
                                index_col=False)

    set_l = d_nonSV_right["length"].unique()
    mean_l = d_nonSV_right["length"].mean()
    median_l = d_nonSV_right["length"].median()

    sample_stats = {'sample': sample,
                    'unmatched': unmatch_c,
                    'SV': sv_c,
                    'l>=1000': long_c,
                    'l<76': short_c,
                    'Y': chrY_c,
                    'X': chrX_c,
                    'mean': mean_l,
                    'median': median_l,
                    'unique': set_l}
    time2 = time.time()
    print("--- Calculate statistics: %2f seconds ---" % (time2 - time1))
    return sample_stats


# @profile
def binning(workdir, sample, binsize):
    print(">>>>>>> Sample: " + sample)
    start_time = time.time()

    df = pd.read_csv(workdir + sample + "_point.csv",
                     sep="\t",
                     names=(("c1", "e1", "length")),
                     index_col=False,
                     dtype={"c1": 'int8', "e1": 'int32', "length": 'int16'})  ## 85% reduction in memory usage!
    df = df.sort_values(by=['c1', 'e1'])
    df["bin"] = (df["e1"] // binsize).astype('int16')

    df_binned = df.groupby(['c1', 'bin'])
    binned_stats_0 = df_binned['length'].agg(f=(lambda x: (x < 140).sum()),
                                             g=(lambda x: ((x < 220) & (x >= 140)).sum()),
                                             h=(lambda x: ((x < 400) & (x >= 220)).sum()),
                                             i=(lambda x: (x >= 400).sum())).reset_index()
    print(binned_stats_0)
    binned_stats_1 = df_binned['length'].agg(['mean', 'median', 'min',
                                              'max', 'count']).reset_index()
    print(binned_stats_1)
    binned_stats = binned_stats_0.merge(binned_stats_1)
    df.rename(columns={"f": "[76,140)", "[140,220)": "b", "h": "[220,400)",
                       "i": "[400,1000)"})

    binned_stats.to_csv(workdir + sample + "_binned.csv",
                        index=False,
                        header=True,
                        sep="\t")

    # for ML matrix, fragment lengths
    binned_lengths = df.groupby(['c1', 'bin']).length.apply(list).reset_index()
    binned_lengths.to_csv(workdir + sample + "_binned_lengths.csv",
                          index=False,
                          header=True,
                          sep="\t")
    time3 = time.time()
    print("--- Binning: %2f seconds ---" % (time3 - start_time))
    # binned_stats_0 = df_binned['length'].agg(f=(lambda x: (x < 140).sum()),
    #                           g=(lambda x: ((x < 220)&(x >=140)).sum()),
    #               h=(lambda x: ((x < 400)&(x >=220)).sum()),
    #               i=( lambda x: (x >= 400).sum())).reset_index()
    # binned_stats_1 = df_binned['length'].agg(['mean','median','min',
    #                                           'max','count']).reset_index()
    # binned_stats = binned_stats_0.merge(binned_stats_1)
    # df.rename(columns={"f": "[76,140)", "[140,220)": "b", "h": "[220,400)",
    # "i": "[400,1000)"})

    # binned_stats.to_csv(workdir + sample + "_binned.csv",
    #                  index=False,
    #                  header=True,
    #                  sep="\t")
    # time3 = time.time()
    # print("--- Binning: %2f seconds ---" % (time3 - start_time))


def main():
    args = get_args().parse_args()
    workdir = args.directory
    samples = args.input

    chrs = []
    for number in range(1, 23):
        chrs += ("".join(("", str(number))),)

    stats = []
    for sample in samples:
        sample_stats = stratify_reads(workdir, sample, chrs)
        stats.append(sample_stats)

    stats = pd.DataFrame(stats)
    stats.to_csv(workdir + "Stats.csv",
                 index=False,
                 sep="\t")

    binsize = 1000000
    for sample in samples:
        binning(workdir, sample, binsize)

    dfs = []
    for sample in samples:
        print(">>>>>>> Sample: " + sample)
        start_time = time.time()
        df = pd.read_csv(workdir + sample + "_binned.csv",
                         sep="\t",
                         names=(("c1", "bin", "f", "g", "h", "i", "mean", "median",
                                 "min", "max", "cou")),
                         index_col=False,
                         skiprows=1)
        nc = "normcount" + sample
        fg = "fperg" + sample
        gh = "gperh" + sample
        df[nc] = df.cou / df.cou.mean()
        df[fg] = df.f / df.g
        df[gh] = df.g / df.h
        dfraw = df[['c1', 'bin', "normcount" + sample, "fperg" + sample,
                    "gperh" + sample]]
        dfs.append(dfraw)
    df_merged = reduce(lambda left, right: pd.merge(left, right, on=["c1", "bin"],
                                                    how='inner'), dfs)
    df_merged.to_csv(workdir + "mergeddf.csv",
                     index=False,
                     sep="\t")


if __name__ == "__main__":
    main()
