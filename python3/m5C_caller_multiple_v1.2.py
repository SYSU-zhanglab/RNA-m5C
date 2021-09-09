#!/usr/bin/env python
import os,sys
import argparse
from collections import defaultdict,Counter
import time
from time import gmtime, strftime
import numpy as np
import pandas as pd
import math
import scipy.stats
import signal
import multiprocessing
from multiprocessing import Process,Pool,Manager
from itertools import chain
import gzip
import csv

class database_gzip_file(object):
	def __init__(self, data, mode="txt"):
		self.mode = mode
		if self.mode == "gzip":
			self.file = gzip.open(data, "rt")
			self.reader = csv.reader(self.file, delimiter="\t")
		elif self.mode == "txt":
			self.file = open(data, "r")
	def write(self):
		pass
	def flush(self):
		pass
	def close(self):
		self.file.close()
	def __exit__(self, exctype, value, traceback):
		if self.file != None:
			self.file.close()
	def __enter__(self):
		return self
	def readline(self):
		try:
			if self.mode == "gzip":
				line = next(self.reader)
				if not line:
					return None
				else:
					return line
			elif self.mode == "txt":
				line = self.file.readline()
				if not line:
					return None
				else:
					return line.strip().split("\t")
		except StopIteration:
			return None

def open_database(fn):
	if fn.endswith(".gz"):
		database_obj = database_gzip_file(fn, mode="gzip")
	else:
		database_obj = database_gzip_file(fn, mode="txt")
	return database_obj

def read_CR(input,nonCRs,control):
	overall = []
	n = 0
	line = input.readline()
	while (line):
		if line[0].startswith("#"):
			# line = line.strip().split("\t")
			if len(line) == 4:
				n += 1
				gene,total,Cs,nonCR = line
				gene = gene.replace("#","")
				if int(total) < options.CT_number or float(nonCR) >= options.gene_CR:
					nonCRs[gene] = 1.0
				else:
					nonCRs[gene] = float(nonCR)
					if gene != "ALL":
						overall.append(float(nonCR))
				if gene in control:
					if int(total) > options.CT_number:
						control[gene] = float(nonCR)
			else:
				stat,value = line
		else:
			nonCRs["Median"] = np.median(overall)
			nonCRs["Mean"] = np.mean(overall)
			nonCRs["discard"] = 1.0
			break
		line = input.readline() #Note that the first line of data was read
	control_nonCR = []
	if control:
		for key,value in control.items():
			if value != "None":
				control_nonCR.append(value)
		nonCRs["control"] = np.median(control_nonCR)
	return line,input,nonCRs,control

def call_m5C(line,nonCRs,CR,C_cutoff,control):
	col = C_cutoff_cols.get(C_cutoff)
	chr,pos_1,dir,gene,name,trans,isoform,biotype,total_cov,CT,A_count,T_count,C_count,G_count = line[0:14]
	total_cov = int(total_cov)
	CT = int(CT)
	if col == -1:
		total_cov_col = total_cov
		CT_col = CT
		C_count_col = int(C_count)
	else:
		total_cov_col,CT_col,C_count_col = line[col].split(";")[1].split(",")
		total_cov_col = int(total_cov_col)
		CT_col = int(CT_col)
		C_count_col = int(C_count_col)
	
	if CT_col >= options.coverage and C_count_col >= options.count:
		ratio = C_count_col / (CT_col + 0.0)
		if ratio < options.ratio or CT/(total_cov+0.0) < options.var_ratio or CT_col/(total_cov_col+0.0) < options.var_ratio:
			return None
			
		nonCR_gene = nonCRs.get(gene)
		
		if nonCR_gene is None:
			nonCR_gene = nonCRs[options.non_anno]
			
		if CR == "gene":
			nonCR = nonCR_gene
		elif CR == "overall":
			if nonCR_gene == 1.0:
				nonCR = 1.0
			else:
				nonCR = nonCRs.get("ALL")
		elif CR == "control":
			if nonCR_gene == 1.0:
				nonCR = 1.0
			else:
				nonCR = nonCRs.get("control")

		if options.method == "binomial":
			pvalue = scipy.stats.binom_test(C_count_col, n=CT_col, alternative='greater', p=nonCR)
		elif options.method == "fisher":
			pass
		elif options.method == "poisson":
			pvalue = scipy.stats.poisson.sf(C_count_col, int(ceil(nonCR*CT_col)))

		if pvalue < options.pvalue:
			Signal = total_cov_col/(total_cov+0.0)
			if Signal >= options.signal:
				if nonCR_gene < 0.00001:
					nonCR_format = format(nonCR_gene,'.1e')
				else:
					nonCR_format = str(round(nonCR_gene,5))
				#LINE = (chr,pos_1,dir,gene,name,trans,isoform,biotype,nonCR_format,total_cov_col,CT_col,C_count_col,ratio,pvalue,Signal)

				return (chr,pos_1,dir)
	else:
		return None

		
def fetch_m5C(line,nonCRs,CR,C_cutoff,control):
	col = C_cutoff_cols.get(C_cutoff)
	chr,pos_1,dir,gene,name,trans,isoform,biotype,total_cov,CT,A_count,T_count,C_count,G_count = line[0:14]
	if (chr,pos_1,dir) in candidates:
		total_cov = int(total_cov)
		CT = int(CT)
		if col == -1:
			total_cov_col = total_cov
			CT_col = CT
			C_count_col = int(C_count)
		else:
			total_cov_col,CT_col,C_count_col = line[col].split(";")[1].split(",")
			total_cov_col = int(total_cov_col)
			CT_col = int(CT_col)
			C_count_col = int(C_count_col)
		nonCR_gene = nonCRs.get(gene)
		if nonCR_gene is None:
				nonCR_gene = nonCRs[options.non_anno]
		if CR == "gene":
			nonCR = nonCR_gene
		elif CR == "overall":
			if nonCR_gene == 1.0:
				nonCR = 1.0
			else:
				nonCR = nonCRs.get("ALL")
		elif CR == "control":
			if nonCR_gene == 1.0:
				nonCR = 1.0
			else:
				nonCR = nonCRs.get("control")
		if CT_col < 10:
			ratio = np.nan
		else:
			ratio = C_count_col / (CT_col + 0.0)
		
		try:
			if CT/(total_cov+0.0) < options.var_ratio or CT_col/(total_cov_col+0.0) < options.var_ratio:
				nonCR = 1.0
		except ZeroDivisionError:
			nonCR = np.nan
		if ratio != np.nan:
			if options.method == "binomial":
				pvalue = scipy.stats.binom_test(C_count_col, n=CT_col, alternative='greater', p=nonCR)
			elif options.method == "fisher":
				pass
			elif options.method == "poisson":
				pvalue = scipy.stats.poisson.sf(C_count_col, int(ceil(nonCR*CT_col)))
		else:
			pvalue = np.nan
		try:
			Signal = total_cov_col/(total_cov+0.0)
		except ZeroDivisionError:
			Signal = np.nan
			
		if nonCR_gene < 0.00001:
			nonCR_format = format(nonCR_gene,'.1e')
		else:
			nonCR_format = str(round(nonCR_gene,5))
		LINE = (chr,pos_1,dir,gene,name,trans,isoform,biotype,nonCR_format,total_cov_col,CT_col,C_count_col,ratio,pvalue,Signal)
		return LINE
	else:
		return None
		
def search_for_candidates(args):
	fn,name,CR,cutoff = args
	
	# with open(fn,'r') as input:
	with open_database(fn) as input:
		nonCRs = {}
		results = []
		control = {}
		if CR == "control":
			if not options.control:
				raise IOError("control list not given, exit.\n")
			else:
				with open(options.control,'r') as input:
					for line in input.readlines():
						# control[line.strip()] = "None"
						control[line[0]] = "None"
		line,input,nonCRs,control = read_CR(input,nonCRs,control)
		while (line):
			# line = line.strip().split("\t")
			result = call_m5C(line,nonCRs,CR,cutoff,control)
			if result is not None:
				results.append(result)
			line = input.readline()
	sys.stderr.write("[%s] %s called: %d\n" % (strftime("%Y-%m-%d %H:%M:%S", time.localtime()),name,len(results)))

	return results

def fetch_candidates_with_list(args):
	fn,name,CR,cutoff = args
	with open_database(fn) as input:
		nonCRs = {}
		results = []
		control = {}
		if CR == "control":
			if not options.control:
				raise IOError("control list not given, exit.\n")
			else:
				with open(options.control,'r') as input:
					for line in input.readlines():
						# control[line.strip()] = "None"
						control[line[0]] = "None"
		line,input,nonCRs,control = read_CR(input,nonCRs,control)
		while (line):
			# line = line.strip().split("\t")
			result = fetch_m5C(line,nonCRs,CR,cutoff,control)
			if result is not None:
				results.append(result)
			line = input.readline()
	#init
	data = {}
	data["info"] = {"gene_id":{}, 'gene':{}}
	data[name]={"coverage":{},
				"C count":{},
				"m5C level":{},
				"P-value":{},
				"signal":{},
				"nonCR":{}
				}

	for records in results:
		chr,pos_1,dir,gene,gene_name,trans,isoform,biotype,nonCR_format,total_cov_col,CT_col,C_count_col,ratio,pvalue,Signal = records
		key = "@".join([chr,pos_1,dir])
		
		data["info"]["gene_id"][key] = gene
		data["info"]["gene"][key] = gene_name
		data[name]["coverage"][key] = CT_col
		data[name]["C count"][key] = C_count_col
		data[name]["m5C level"][key] = ratio
		data[name]["P-value"][key] = pvalue
		data[name]["signal"][key] = Signal
		data[name]["nonCR"][key] = nonCR_format

	df = pd.DataFrame.from_dict({(i,j): data[i][j] 
								  for i in data.keys() 
								  for j in data[i].keys()},
								  orient='columns')
	csv_name = name+"_temp_"+str(os.getpid())+".csv"
	df.to_csv(csv_name)
	return csv_name
	
def reorder_columns(cols):
	order = {'coverage':1,
			 'C count' :2,
			 'm5C level':3,
			 'P-value':4,
			 'signal':5,
			 'nonCR' :6,
	         }
	new_cols = sorted(cols,key=lambda x:[x[0],order.get(x[1])])
	return new_cols
	
def merge_sperated_lists(csvs):
	dict_csv = {}
	for csv_name in csvs:
		dict_csv[csv_name] = pd.read_csv(csv_name,index_col=0,header=[0,1])
		dict_csv[csv_name].set_index([('info','gene_id'),('info','gene')],append=True,inplace=True)
	df = pd.concat(dict_csv.values(),axis=1)
	df.index.rename(["","",""],inplace=True)
	df = df[reorder_columns(df.columns)]
	if options.withTemp == False:
		for csv_name in csvs:
			os.remove(csv_name)
	return df
	
def read_candidate_list():
	temp = []
	with open(options.list,'r') as input:
		try:
			for line in input.readlines():
				line = line.strip().split("\t")
				chr,pos,dir = line
				temp.append((chr,pos,dir))
		except IndexError:
			raise Warning("Site list should in format (TAB seperated): chr pos(1-based) dir\n")
	return set(temp)

def read_fn_list():
	search_list = []
	fetch_list = []
	with open(options.input,'r') as input:
		for line in input.readlines():
			line = line.strip().split("\t") #name	directory	CR mode	C-cutoff
			name,dir,CR,cutoff = line
			search_list.append([dir,name,CR,cutoff])
	return search_list
	
def signal_handler(sig,frame):
	pool.terminate()
	sys.exit()
	
if __name__ == "__main__":
	description = """
	Step1: read through all of the files, find m5C sites, and merge the list (optional)

	Step2: Use the site list (or a given one) to fetch sites
	
	Step3: Merge tables
	"""
	parser = argparse.ArgumentParser(prog="m5C_caller_multiple",fromfile_prefix_chars='@',description=description,formatter_class=argparse.RawTextHelpFormatter)
	#Require
	group_required = parser.add_argument_group("Required")
	group_required.add_argument("-i","--input",dest="input",required=True,help="a list. Format as: name\tdir\tCR\tC-cutoff")
	group_required.add_argument("-o","--output",dest="output",required=False,default="multiple.csv",help="output name, default=multiple.csv")
	#Running
	group_optional = parser.add_argument_group("Optional")
	group_optional.add_argument("-P","--processors",dest="processors",default=1,type=int,help="Processor number, default=1")
	group_optional.add_argument("-l","--list",dest="list",required=False,help="site list")
	group_optional.add_argument("--with-temp",dest="withTemp",default=False,help="do not delete temp files")
	#Filter
	group_site = parser.add_argument_group("m5C filter")
	group_site.add_argument("-c","--coverage",dest="coverage",default=20,type=int,help="C+T coverage, default=20")
	group_site.add_argument("-C","--count",dest="count",default=3,type=int,help="C count, below which the site will not count, default=3")
	group_site.add_argument("-r","--ratio",dest="ratio",default=0.1,type=float,help="m5C level/ratio, default=0.1")
	group_site.add_argument("-p","--pvalue",dest="pvalue",default=0.05,type=float,help="pvalue, default=0.05")
	group_site.add_argument("-s","--signal",dest="signal",default=0.9,type=float,help="signal ratio, equals coverage(under C-cutoff)/coverage, default=0.9")
	group_site.add_argument("-R","--var_ratio",dest="var_ratio",default=0.8,type=float,help="the ratio cutoff of CT/Total to filter sequencing/mapping errors, default=0.8")
	group_site.add_argument("-g","--gene-CR",dest="gene_CR",default=0.05,type=float,help="conversion rate, below which a gene will be discarded, default=0.05")
	group_site.add_argument("-N","--CT",dest="CT_number",default=1000,type=int,help="CT count, below which a gene will be discarded, default=1000")
	#Statistics
	group_stat = parser.add_argument_group("Statistic method")
	group_stat.add_argument("--method",dest="method",default="binomial",choices=['binomial', 'fisher', 'poisson'],help="statistical method: binomial, fisher exact test, or poisson, default=binomial")
	#group_stat.add_argument("--CR",dest="conversion_rate",default="gene",choices=['gene','overall','control'],help="conversion rate used: gene or overall, default=gene")
	group_stat.add_argument("--NA",dest="non_anno",default="ELSE",choices=['ELSE','Median','Mean','ALL',"discard"],help="which CR to use if no aene annotation, default=ELSE")
	group_stat.add_argument("--control",dest="control",help="control list, median non-conversion rate will be used")
	#Version
	group_other = parser.add_argument_group("Other")
	group_other.add_argument("--version",action="version",version="%(prog)s 1.0")
	options = parser.parse_args()
	
	sys.stderr.write("[%s] CMD: %s\n" % (strftime("%Y-%m-%d %H:%M:%S", time.localtime())," ".join(sys.argv)))

	sys.stderr.write("[%s] Running with %d processors, pid [%s]\n" % (strftime("%Y-%m-%d %H:%M:%S", time.localtime()),options.processors,str(os.getpid())))


	# C_cutoff = options.C_cutoff

	# if C_cutoff not in ['1','2','3','4','5','6','7','8','9','10','15','20','None']:
		# raise ValueError("%s not in C-cutoffs, exit.\n" % C_cutoff) 
			
	# sys.stderr.write("[%s] C-cutoff = [%s]\n" % (strftime("%Y-%m-%d %H:%M:%S", time.localtime()),C_cutoff))
	
	C_cutoff_cols = {}
	col_start = 14
	n = 0
	for i in ['1','2','3','4','5','6','7','8','9','10','15','20','None'][:-1]:
		C_cutoff_cols[i] = 14 + n
		n += 1
	C_cutoff_cols["None"] = -1
	candidates = {}
	#init
	search_list = read_fn_list()
	
	#init pool
	signal.signal(signal.SIGINT,signal_handler)
	
	
	sys.stderr.write("[%s] Searching for m5C candidates...\n" % (strftime("%Y-%m-%d %H:%M:%S", time.localtime())))
	#Step 1
	if not options.list: #search for candidates first
		pool = multiprocessing.Pool(options.processors)
		try:
			candidates_list = pool.map(search_for_candidates,search_list)
			candidates_set = set(chain(*candidates_list))
		finally:
			pool.terminate()
	else:
		candidates_set = read_candidate_list()
	
	if not candidates_set:
		sys.stderr.write("No candidates found, exit.\n")
		pool.terminate()
		sys.exit()
	else:
		candidates = {i:1 for i in candidates_set}

	sys.stderr.write("[%s] Getting common sites...\n" % (strftime("%Y-%m-%d %H:%M:%S", time.localtime())))
	#Step2, re-build a pool
	pool = multiprocessing.Pool(options.processors)
	csvs = []
	try:
		csvs = pool.map(fetch_candidates_with_list,search_list) #directly output to table
	finally:
		pool.terminate()
	
	sys.stderr.write("[%s] Merge tables...\n" % (strftime("%Y-%m-%d %H:%M:%S", time.localtime())))
	#Step3, merge dataframe
	df = merge_sperated_lists(csvs)
	df.to_csv(options.output)
	sys.stderr.write("[%s] Finished successfully!\n" % (strftime("%Y-%m-%d %H:%M:%S", time.localtime())))
