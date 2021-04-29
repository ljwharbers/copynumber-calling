#!/usr/bin/env python2

import sys
import time
from operator import itemgetter


def main():

	infilename = "/mnt/AchTeraD/Documents/references/hg19_cnv/chrY.fa"
	outfilename = "/mnt/AchTeraD/Documents/references/hg19_cnv/chrY.psr.fa"

	INFILE = open(infilename, "r")
	OUTFILE = open(outfilename, "w")
	
	chrY = ""
	x = INFILE.readline()
	header_line = x.rstrip()
	for x in INFILE:
		aline = x.rstrip()
		chrY += aline

	newChrY = chrY[0:10000] + "N" * (2649520 - 10000) + chrY[2649520:59034049] + "N" * (59363566 - 59034049) + chrY[59363566:]
	print len(chrY)
	print len(newChrY)

	chrY = newChrY	

	OUTFILE.write(header_line)
	OUTFILE.write("\n")
	lines = len(chrY) / 50
	remainder = len(chrY) - (lines * 50)
	for i in range(lines):
		OUTFILE.write(chrY[(i * 50):((i * 50) + 50)])
		OUTFILE.write("\n")
	OUTFILE.write(chrY[(lines * 50):len(chrY)])
	OUTFILE.write("\n")


if __name__ == "__main__":
	main()
