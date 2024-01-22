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
	dictCDS_start = {}
	dictCDS_end = {}
	dictBiotype = {}
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
					gene_biotype = row[row.index("gene_type")+1].replace('\"','').strip(";")
				except ValueError:
					gene_biotype = "NA"
				
				try:
					transcript_biotype = row[row.index("transcript_type")+1].replace('\"','').strip(";")
				except ValueError:
					transcript_biotype = gene_biotype

				try:
					gene_name = row[row.index("gene_name")+1].replace('\"','').strip(";")
				except ValueError:
					gene_name = gene_id
				dictBiotype[trans_id] = gene_biotype+","+transcript_biotype
			
				info = (trans_id,gene_id,gene_name,chr,dir)
				type = row[2]
				if not dictInfo.has_key(trans_id):
					dictInfo[trans_id]={'start':[],
										'end':[],
										'dir':"",
										'chr':"",
										'gene_id':"",
										'gene_name':""}
					dictInfo[trans_id]['dir'] = dir
					dictInfo[trans_id]['chr'] = chr
					dictInfo[trans_id]['gene_id'] = gene_id
					dictInfo[trans_id]['gene_name'] = gene_name
				if type == "exon" or type == "Exon" or type == "EXON" or type == "tRNAscan":
					dictInfo[trans_id]['start'] = dictInfo[trans_id]['start'] + [start]
					dictInfo[trans_id]['end'] = dictInfo[trans_id]['end'] + [end]
				elif type == "CDS":
					if trans_id not in dictCDS_start:
						dictCDS_start[trans_id] = start
					elif start < dictCDS_start[trans_id]:
							dictCDS_start[trans_id] = start
					if trans_id not in dictCDS_end:
						dictCDS_end[trans_id] = end
					elif dictCDS_end[trans_id] < end:
						dictCDS_end[trans_id] = end
					
				line = gtf.readline()
			except ValueError:
				line = gtf.readline()
	
	for trans_id in dictInfo.keys():
		exonCount = str(len(dictInfo[trans_id]['start']))
		exonStarts = ','.join([str(i) for i in sorted(dictInfo[trans_id]['start'])]) + ","
		exonEnds = ','.join([str(i) for i in sorted(dictInfo[trans_id]['end'])]) + ","
		txStart = sorted(dictInfo[trans_id]['start'])[0]
		txEnd  = sorted(dictInfo[trans_id]['end'])[-1]
		cdsStart= dictCDS_start.get(trans_id)
		cdsEnd = dictCDS_end.get(trans_id)
		if cdsStart is None:
			cdsStart = "."
		if cdsEnd is None:
			cdsEnd  = "."
		ident = str(dictBiotype.get(trans_id))
		print "\t".join([ident,trans_id,dictInfo[trans_id]['chr'],dictInfo[trans_id]['dir'],str(txStart),str(txEnd),str(cdsStart),str(cdsEnd),str(exonCount),exonStarts,exonEnds,".",dictInfo[trans_id]['gene_id'],dictInfo[trans_id]['gene_name'],'.','.'])
