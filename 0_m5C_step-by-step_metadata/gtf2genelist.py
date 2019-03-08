#!bin/usr/env python
import sys
import os
from optparse import OptionParser
from Bio import SeqIO
	
if __name__ == "__main__":
	#Parser
	usage = "USage: python %prog -i <gtf> -f <fasta>"
	parser = OptionParser(usage=usage)
	parser.add_option("-i",dest="input",help="GTF file")
	parser.add_option("-f",dest="fasta",help="Fasta file")
	(options,args) = parser.parse_args()
	dictLength = {}
	
	for seq in SeqIO.parse(options.fasta,"fasta"):
		dictLength[seq.id] = len(seq.seq)
	
	with open(options.input,'r') as gtf:
		line = gtf.readline()
		list = []
		while(line):
			if line.startswith("#"):
				line = gtf.readline()
				continue
			row = line.split(	)
			type = row[2]
			if type == "transcript":
				chr = row[0]
				dir = row[6]
				start = row[3]
				end = row[4]
				gene_id = row[row.index("gene_id")+1].replace('\"','').strip(";")
				try:
					gene_name = row[row.index("gene_name")+1].replace('\"','').strip(";")
				except ValueError:
					gene_name = gene_id
				trans_id = row[row.index("transcript_id")+1].replace('\"','').strip(";")
				try:
					trans_name = row[row.index("transcript_name")+1].replace('\"','').strip(";")
				except ValueError:
					trans_name = trans_id
				try:
					gene_biotype = row[row.index("gene_biotype")+1].replace('\"','').strip(";")
				except ValueError:
					gene_biotype = "unknown"
				length = str(dictLength.get(trans_id))
				list.append((trans_name,gene_name,trans_id,gene_id,gene_biotype,chr,start,end,dir,length))
			line = gtf.readline()
	list_sorted = sorted(list,key=lambda x:(x[0],x[1]))
	for i in list_sorted:
		print "\t".join(i)
