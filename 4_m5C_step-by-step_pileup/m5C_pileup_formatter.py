#!/usr/bin/env python
import os,sys
import argparse
from collections import defaultdict,Counter
import copy
import sqlite3
import time
from time import gmtime, strftime
import numpy as np
def write_CR():
	#ALL	900365017	2875105	0.00319326600403
	#ELSE	85530991	352474	0.00412100919069
	#ENSG00000000003	25890	148	0.00571649285438
	overall = []
	lines = []
	with file(options.conversion,'w') as CR_output:
		CR_output.write("\t".join(["#ALL",str(all_cov),str(all_C),str(all_C/(all_cov+0.0))]))
		CR_output.write("\n")
		for key in sorted(result_cov.keys()):
			C = result_C.get(key)
			cov = result_cov.get(key)
			if cov is not None and cov > 0:
				if C is None:
					C = 0
				nonCR = C/(cov+0.0)
				line = "\t".join([key,str(cov),str(C),str(nonCR)]) + "\n"
				lines.append(line)
				if cov > 1000:
					overall.append(nonCR)
		CR_output.write("#Median\t%.8f\n" % np.median(overall))
		CR_output.write("#Mean\t%.8f\n" % np.mean(overall))
		CR_output.write("#90%%\t%.8f\n" % np.percentile(overall,90))
		CR_output.write("#75%%\t%.8f\n" % np.percentile(overall,75))
		CR_output.write("#50%%\t%.8f\n" % np.percentile(overall,50))
		CR_output.write("#25%%\t%.8f\n" % np.percentile(overall,25))
		CR_output.write("#10%%\t%.8f\n" % np.percentile(overall,10))
		for line in lines:
			CR_output.write("#")
			CR_output.write(line)

def create_connection(db_file):
	if os.path.isfile(db_file):
		conn = sqlite3.connect(db_file)
	else:
		raise IOError("Database file not exist.\n")
	return conn

def search(cursor,chr,dir,pos):
	cursor.execute("SELECT gene,name,trans,isoform,biotype FROM locations WHERE chr=? AND dir=? AND pos=?",(chr,dir,pos))
	rows = cursor.fetchall()
	if rows:
		return [str(i) for i in rows[0]]
	else:
		return None

def Tideup(chr,pos_1,dir,PIPLEUPS,SURROUNDINGS):
	global cursor,output,all_C,all_cov
	#print 2,len(PIPLEUPS),len(SURROUNDINGS)
	query = search(cursor,chr,dir,pos_1)
	if query is not None:
		gene,name,trans,isoform,biotype = query
		GENE = gene
	else:
		gene,name,trans,isoform,biotype = ["NA"] * 5
		GENE = "ELSE"
	bases = sorted(Counter(zip(PIPLEUPS,SURROUNDINGS)).items(),key=lambda x:x[0][1])
	counts = {'A':0,'T':0,'C':0,'G':0,'total':0,'CT':0}
	bins = {i: {'A':0,'T':0,'C':0,'G':0,'total':0,'CT':0} for i in [1,2,3,4,5,6,7,8,9,10,15,20]}
	limits = [99999, 20, 15, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1]
	for (base,surr),count in bases:
		if surr > limits[-1]:
			while surr > limits[-1]:
				index = limits.pop()
				bins[index] = copy.copy(counts)
		
		if base != "N":
			counts["total"] += count
			counts[base] += count
		if base == "C" or base == "T":
			counts["CT"] += count
	for item in limits:
		bins[item] = copy.copy(counts)
	#output.write("\t".join([str(i) for i in limits])+"\n")
	#write to disk
	output.write("\t".join([chr,str(pos_1),dir,gene,name,trans,isoform,biotype])) #information
	output.write("\t")
	output.write("\t".join([str(counts['total']),str(counts['CT']),str(counts['A']),str(counts['T']),str(counts['C']),str(counts['G'])])) #Overall count
	for i in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 99999][:-1]:
		output.write("\t")
		output.write(str(i))
		output.write(";")
		output.write(",".join([str(bins[i]['total']),str(bins[i]['CT']),str(bins[i]['C'])])) #SP=bins[i]['total']/counts['total']
	output.write("\n")
	
	#Add values to result_cov and result_C
	result_C[GENE] += counts['C']
	result_cov[GENE] += counts['CT']
	all_C += counts['C']
	all_cov += counts['CT']
	
def check(chr,pos_1,dir,pileup,surrounding):
	global PIPLEUPS,SURROUNDINGS,last
	if not last: #first in
		last = (chr,pos_1,dir)
		PIPLEUPS.extend(pileup)
		SURROUNDINGS.extend(surrounding)
	else:
		if last == (chr,pos_1,dir): #keep putting items
			PIPLEUPS.extend(pileup)
			SURROUNDINGS.extend(surrounding)
		else: #To output
			#print 1,len(PIPLEUPS),len(SURROUNDINGS)
			Tideup(last[0],last[1],last[2],PIPLEUPS,SURROUNDINGS)
			PIPLEUPS = []
			SURROUNDINGS = []
			last = (chr,pos_1,dir)
			PIPLEUPS.extend(pileup)
			SURROUNDINGS.extend(surrounding)
	
if __name__ == "__main__":
	description = """
	"""
	parser = argparse.ArgumentParser(prog="",fromfile_prefix_chars='@',description=description,formatter_class=argparse.RawTextHelpFormatter)
	#Require
	group_require = parser.add_argument_group("Required")
	group_require.add_argument("-i","--input",dest="input",required=True,help="input pileup, sorted")
	group_require.add_argument("-o","--output",dest="output",required=True,help="output") #chr start_0 start
	group_require.add_argument("--CR",dest="conversion",required=True,help="gene conversion rate")
	group_require.add_argument("--db",dest="database",required=True,help="database, base-gene annotation")
	group_optional = parser.add_argument_group("Optional")
	group_optional.add_argument("-c","--coverage",dest="cov",default=10,type=int,help="coverage, default=10")
	group_optional.add_argument("-r","--ratio",dest="ratio",default=0.1,type=float,help="m5C level/ratio, default=0.1")
	group_optional.add_argument("-p","--pvalue",dest="pvalue",default=0.05,type=float,help="binomial pvalue, default=0.05")
	group_other = parser.add_argument_group("Other")
	options = parser.parse_args()
	
	sys.stderr.write("[%s] Reading input...\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
	
	database = create_connection(options.database)
	cursor = database.cursor()
	# cursor.execute("PRAGMA cache_size =8000")
	# cursor.execute("PRAGMA page_size = 8192")
	# cursor.execute("PRAGMA temp_store =2")
	# cursor.execute("PRAGMA synchronous=OFF")
	cursor.execute("PRAGMA journal_mode = WAL") #do not need a write lock while reading, so enable multiple reading
	
	#init
	result_cov = defaultdict(int)
	result_C = defaultdict(int)
	all_cov = 0
	all_C = 0
	PIPLEUPS= []
	SURROUNDINGS = []
	last = None
	with open(options.input,'r') as input,file(options.output+"_tmp",'w') as output:
		#10      100007581(1-based)       -       C       ENST00000260702 3519    T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
		#10      100007581       -       C       Genome  100007581       T,T,T,T,T,T,T,T,T,T,T   0,0,0,0,0,0,0,0,0,0,0
		line = input.readline()
		while (line):
			row = line.strip().split("\t")
			chr = row[0]
			pos_1 = row[1]
			dir = row[2]
			base = row[3] # as transcript
			pileup = row[6].split(",")
			surrounding = [int(i) for i in row[7].split(",")]
			if base == "C":
				check(chr,pos_1,dir,pileup,surrounding)
			line = input.readline()
		Tideup(last[0],last[1],last[2],PIPLEUPS,SURROUNDINGS)
	sys.stderr.write("[%s] Finished reading input.\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
	#get conversion rates
	write_CR()
	sys.stderr.write("[%s] Calculating conversion rates...\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
	os.system("cat %s %s > %s" % (options.conversion,options.output+"_tmp",options.output))
	os.remove(options.output+"_tmp")
	sys.stderr.write("[%s] All finished!\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))