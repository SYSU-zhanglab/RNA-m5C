#!bin/usr/env python

#Jianheng Liu @ Zhanglab, SYSU
#Feb, 2018
#Email: liujh26@mail2.sysu.edu.cn
#Usage: Extract GTF annotation to annotation table (UCSC table like)
#Input: [.gtf]

import sys
import os
from optparse import OptionParser

if __name__ == "__main__":
	#Parser
	usage = "Usage: python %prog -i <gtf> > output.anno"
	parser = OptionParser(usage=usage)
	parser.add_option("-i",dest="input",help="Input file")
	(options,args) = parser.parse_args()
	#format:
	#transcript_id gene_id gene_name chr dir 5_UTR_len CDS_len 3_UTR_len 5_UTR start_codon CDS stop_codon 3_UTR 5_len CDS_len 3_len Exon
	dictInfo = {}
	with open(options.input,'r') as gtf:
		line = gtf.readline()
		while(line):
			if line.startswith("#"):	
				line = gtf.readline()
				continue
			row = line.split(	)
			chr = row[0]
			dir = row[6]
			start,end = int(row[3])-1,int(row[4])
			try:
				trans_id = row[row.index("transcript_id")+1].replace('\"','').strip(";")
				gene_id = row[row.index("gene_id")+1].replace('\"','').strip(";")
				try:
					gene_name = row[row.index("gene_name")+1].replace('\"','').strip(";")
				except ValueError:
					gene_name = gene_id
				info = (trans_id,gene_id,gene_name,chr,dir)
				type = row[2]
				if not dictInfo.has_key(trans_id):
					dictInfo[trans_id]={'start':[],
										'end':[],
										'dir':"",
										'chr':"",
										'gene_id':""}
					dictInfo[trans_id]['dir'] = dir
					dictInfo[trans_id]['chr'] = chr
					dictInfo[trans_id]['gene_id'] = gene_id
				if type == "exon" or type == "tRNAscan":
					dictInfo[trans_id]['start'] = dictInfo[trans_id]['start'] + [start]
					dictInfo[trans_id]['end'] = dictInfo[trans_id]['end'] + [end]
				line = gtf.readline()
			except ValueError:
				line = gtf.readline()
	
	for trans_id in dictInfo.keys():
		exonCount = str(len(dictInfo[trans_id]['start']))
		exonStarts = ','.join([str(i) for i in sorted(dictInfo[trans_id]['start'])]) + ","
		exonEnds = ','.join([str(i) for i in sorted(dictInfo[trans_id]['end'])]) + ","
		print "\t".join(['.',trans_id,dictInfo[trans_id]['chr'],dictInfo[trans_id]['dir'],'txStart','txEnd','cdsStart','cdsEnd',exonCount,exonStarts,exonEnds,".",dictInfo[trans_id]['gene_id'],'.','.','.'])
