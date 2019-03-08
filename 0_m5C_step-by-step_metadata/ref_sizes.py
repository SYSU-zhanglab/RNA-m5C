#!bin/usr/env python

#Jianheng Liu @ Zhanglab, SYSU
#Feb, 2018
#Email: liujh26@mail2.sysu.edu.cn
#Usage: Extract reference lengths from fasta file(s)
#Input: [.fa]

import argparse
from Bio import SeqIO

if __name__ == "__main__":
	description = """
	"""
	parser = argparse.ArgumentParser(prog="ref_dict",fromfile_prefix_chars='@',description=description,formatter_class=argparse.RawTextHelpFormatter)
	#Require
	group_required = parser.add_argument_group("Required")
	group_required.add_argument("--input","-i",nargs="*",dest="input",required=True,help="input")
	group_required.add_argument("--output","-o",dest="output",required=False,help="output")
	options = parser.parse_args()
	if options.output:
		with file(options.output,'w') as output:
			for fn in options.input:
				for seq in SeqIO.parse(fn,'fasta'):
					length = len(seq.seq)
					output.write(seq.id+"\t"+str(length)+"\n")
	else:
		for fn in options.input:
			for seq in SeqIO.parse(fn,'fasta'):
				length = len(seq.seq)
				print "\t".join([seq.id,str(length)])