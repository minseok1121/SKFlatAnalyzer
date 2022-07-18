import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--target", "-d", required=True, type=str, help="target directory")
parser.add_argument("--output_filename", "-o", required=True, type=str, help="output filename")
parser.add_argument("--era", "-e", required=True, type=str, help="era")
args = parser.parse_args()
target = args.target
output_filename = args.output_filename

main_dirs = os.listdir(mother)
with open(f"{era}/ForSNU/{output_filename}.txt", "w" ) as f:
		for d in main_dirs:
				sub_dirs = os.listdir(f"{target}/{d}")
				for s in sub_dirs:
						f.write(f"{target}/{d}/{s}\n")

