#! /usr/bin/env python

import argparse
import multiprocessing as mp
from datetime import datetime

import numpy as np
import pandas as pd
import sys

global sz_min       # :param min: minimum fragment size
global sz_max       # :param max: maximum fragment size

mp.set_start_method('fork')


def cli_parser():
    p = argparse.ArgumentParser()
    p = argparse.ArgumentParser(description=f'Go over all fragments in a BEDpe file, representing paired reads '
                                            f'and filter fragments by size selection.',
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-i", "--input", help="Input BEDpe file to be filtered.")
    p.add_argument("-o", "--output", help="Output BEDpe file filtered by fragments size selection.")
    p.add_argument("-min", "--min", type=int, help="Minimum fragments size.")
    p.add_argument("-max", "--max", type=int, help="Maximum fragments size.")
    p.add_argument("-n", "--n_cores", default=None, type=int, help="number of cores")
    p.add_argument("-p", "--n_partitions", default=None, type=int, help="number of partitions to split the data")
    return p


def size_selection_df(bedpe_df):
    global sz_min
    global sz_max
    filtered_bed_df = pd.DataFrame(columns=bedpe_df.columns)
    df_len = len(bedpe_df.index)
    print(f'size_selection_df : len={df_len}')
    for i, fragment in bedpe_df.iterrows():
        fragment_chrom1 = fragment['chrom1']
        fragment_chrom2 = fragment['chrom2']
        fragment_start = fragment['start1']
        fragment_end = fragment['end2']
        fragment_len = fragment_end - fragment_start
        # print(f'{fragment["chrom1"]}, {fragment_start}, {fragment_end}: fragment_len={fragment_len}')
        if (fragment_chrom1 == fragment_chrom2) and (sz_min <= fragment_len <= sz_max):
            filtered_bed_df = filtered_bed_df.append(fragment, ignore_index=True)
    return filtered_bed_df


def parallelize_dataframe(df, func, n_cores: int = None, n_partitions: int = None):
    df_len = len(df.index)
    print(f'dataframe len={df_len}')
    n_processes = n_cores if n_cores else min(mp.cpu_count(), df_len)

    if n_partitions is None:
        n_partitions = n_processes

    print(f'Using {n_processes} number of processes...')
    print(f'Splitting data in {n_partitions} partitions...')

    df_split = np.array_split(df, n_partitions)

    # try:
    with mp.Pool(processes=n_processes) as pool:
        results = pool.map(func, df_split)
        df = pd.concat(results)
    # except err:
    #    print('oops')
    #    df = None
    return df


def bed_size_selection(bed_file: str, filtered_bed_file: str, n_cores: int, n_partitions):
    """
    Go over all fragments in a BEDpe file, representing paired reads and filter fragments by size selection.

    :param bed_file: the input BED file to be filtered by size selection fragments
    :param filtered_bed_file: the output BED file filtered by size selection.
    :param n_cores: number of cpu core to be used
    :param n_partitions: number of partitions to split the input BEE file
    """
    bedpe_df = pd.read_table(bed_file, names=['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'name', 'score',
                                            'strand1', 'strand2'],
                             dtype={'start1': int, 'end1': int, 'start2': int, 'end2': int, 'score': int})

    filtered_bed_df = parallelize_dataframe(bedpe_df, size_selection_df, n_cores=n_cores, n_partitions=n_partitions)
    # filtered_bed_df = size_selection_df(bedpe_df)  # without parallelization
    filtered_bed_df.to_csv(filtered_bed_file, sep='\t', na_rep='.', header=False, index=False)
    print(f'Results saved in file {filtered_bed_file}')


def main():
    global sz_min
    global sz_max

    args = cli_parser().parse_args(sys.argv[1:])

    # global params
    sz_min = args.min
    sz_max = args.max

    start = datetime.now()
    current_time = start.strftime("%y-%d-%m %H:%M:%S")
    print("start time: ", current_time)

    bed_size_selection(args.input, args.output, n_cores=args.n_cores, n_partitions=args.n_partitions)

    end = datetime.now()
    time_elapsed = end - start
    print('time elapsed (hh:mm:ss.ms) {}'.format(time_elapsed))
    print("end time: ", end)


if __name__ == '__main__':
    main()

