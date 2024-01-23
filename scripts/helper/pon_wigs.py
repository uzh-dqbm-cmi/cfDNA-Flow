import argparse
import os

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-wd', '--wig_dir_path', default='',
                        help='Directory with wig files, to use as a prefix.')
    parser.add_argument('-cs', '--control_samples', default='',
                        help='Path to the file with a list of control samples.')
    parser.add_argument('-od', '--output_destination', default='',
                        help='Path to the output destination, usually main ichorCNA folder.')
    return parser.parse_args()

def main(args):
    samples = []
    with open(args.control_samples) as fin:
        for line in fin:
            samples.append(line.rstrip())
    fout = open(os.path.join(args.output_destination,"pon_wigs.txt"), 'w')
    for sample in samples:
        full_path = os.path.join(args.wig_dir_path, sample+'.wig')
        fout.write(full_path + '\n')
    fout.close()

if __name__=="__main__":
    args = get_args()
    main(args)