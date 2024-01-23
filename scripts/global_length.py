import pandas as pd
import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
import glob
import argparse
import ast

def get_args():
    parser = argparse.ArgumentParser("cfDNA Feature Extraction\n")
    parser.add_argument("-i", "--input", nargs='+',
                        help="List with paths to sample.csv files with fragment lengths", required=True)
    parser.add_argument("-o", "--output", nargs='+',
                        help="Outputs .csv file with global length and Mouliere's features", required=True)
    return parser


def main():
    args = get_args().parse_args()
    files = args.input
    output_csv = args.output

    tmp = []
    for file in files:

        df = pd.read_csv(file, names=["length"])
        data = df[(df["length"] >= 100) & (df["length"] <= 220)]
        # Mean, median and standard deviation
        mu, std = norm.fit(data)
        med = np.median(data)

        # Mouliere's length features
        p100_150 = len(df[(df["length"] >= 100) & (df["length"] <= 150)]) / len(df)
        p160_180 = len(df[(df["length"] >= 160) & (df["length"] <= 180)]) / len(df)
        p180_220 = len(df[(df["length"] >= 180) & (df["length"] <= 220)]) / len(df)
        p250_320 = len(df[(df["length"] >= 250) & (df["length"] <= 320)]) / len(df)
        p380_500 = len(df[(df["length"] >= 380) & (df["length"] <= 500)]) / len(df)
        p500_1000 = len(df[(df["length"] >= 500) & (df["length"] <= 1000)]) / len(df)
        div_tmp = len(df[(df["length"] >= 163) & (df["length"] <= 169)])
        if div_tmp > 0:
            r_100150_163169 = len(df[(df["length"] >= 100) & (df["length"] <= 150)]) / div_tmp
        else:
            r_100150_163169 = 0.0
        r_160180_180220 = len(df[(df["length"] >= 160) & (df["length"] <= 180)]) / len(
            df[(df["length"] > 180) & (df["length"] <= 220)])

        # binning
        df.loc[:, 'binned'] = pd.cut(df['length'], np.arange(70, 1000, 10).tolist())
        #df.loc[:, 'binned'] = pd.cut(df['length'], [70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200,
                                                    #210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 310, 320, 330,
                                                    #340, 350, 360, 370, 380, 400])
        bin_df = df.groupby('binned').count()
        bin_df['norm_length'] = bin_df['length'] / len(df)
        bin_dict = {k.left: v for k, v in bin_df.to_dict()['norm_length'].items()}

        # saving sample features
        sample_dict = {"sample": file.split("/")[-1].split("_len")[0], "mean": mu, "std": std, "med": med,
                       "p100_150": p100_150, "p160_180": p160_180, "p180_220": p180_220, "p250_320": p250_320,
                       "p380_500": p380_500, "p500_1000": p500_1000,
                       "r_100150_163169": r_100150_163169, "r_160180_180220": r_160180_180220}
        sample_dict.update(bin_dict)
        tmp.append(sample_dict)
        print(f"{file} done.")

    # Save output to tsv file
    output = pd.DataFrame(tmp)
    output = output.add_prefix('globlen_')
    output = output.rename(columns={"globlen_sample": "sample"})
    output.to_csv(output_csv[0], sep="\t", index=False)

if __name__ == "__main__":
    main()