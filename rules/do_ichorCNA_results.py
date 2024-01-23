import argparse
import pandas as pd


def cli_parser():
    p = argparse.ArgumentParser(description=f'Perform Histogram-based Gaussian Mixture Model with fixed means.',
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-i", "--input", nargs='+', help="Input ichorCNA/{sample}/{sample}.params.txt files.")
    p.add_argument("-o", "--output", help=f"Output results file.")
    return p


def main(input_file, output):
    results_df = pd.DataFrame(columns=["sample", "tumor_fraction", "ploidy", "gender"])
    for sample_input in input_file:
        sample_df = pd.read_table(sample_input, sep="\t", header=0, names=["sample", "tumor_fraction", "ploidy"],
                                  skipfooter=18, engine='python')
        # sample_ichor = sample_df[0]
        sample_df["gender"] = ""

        with open(sample_input, 'r') as fp:
            lines = fp.readlines()
            for line in lines:
                if line.find("Gender:") != -1:
                    sample_df["gender"] = line.split(":")[1].strip()
                    break
            results_df = pd.concat((results_df, sample_df), axis=0)
    results_df.to_csv(f'{output}', sep="\t", index=False)


if __name__ == '__main__':
    args = cli_parser().parse_args()
    main(args.input, args.output)
