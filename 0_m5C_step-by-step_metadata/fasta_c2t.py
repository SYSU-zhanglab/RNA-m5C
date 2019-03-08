#!/usr/bin/env python

#Jianheng Liu @ Zhanglab, SYSU
#Feb, 2018
#Email: liujh26@mail2.sysu.edu.cn
#Usage: C2T convert fasta files
#Input: [.fa]

import sys
import os
import getopt
from Bio import SeqIO

opts,args=getopt.getopt(sys.argv[1:],"i:h")
for op,value in opts:
	if op == "-i":
		file = value
	elif op == "-h":
		print "Usage: python fasta_c2t.py -i input.fa > output.fa"
		sys.exit()
	else:
		print "Usage: python fasta_c2t.py -i input.fa > output.fa"
		sys.exit()

for seq_record in SeqIO.parse(file,"fasta"):
	print ">"+seq_record.description
	print str(seq_record.seq).upper().replace("C","T")
