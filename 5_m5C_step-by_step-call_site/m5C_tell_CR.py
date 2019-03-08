#!/usr/bin/env python
import os,sys
import argparse
from collections import defaultdict,Counter
import time
from time import gmtime, strftime
import numpy as np
import pandas as pd

def read_CR(name,input,nonCRs,control):
	overall = []
	n = 0
	control_all = 0
	control_C = 0
	line = input.readline()
	while (line):
		if line.startswith("#"):
			line = line.strip().split("\t")
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
					control_all += int(total)
					control_C  += int(Cs)
					if int(total) > options.CT_number:
						control[gene] = float(nonCR)
				if gene == "ELSE":
					#sys.stdout.write("#No annotation\t%.6f\n" % (float(nonCR)))
					statistics["Non-annotated"][name] = float(nonCR)
				elif gene == "ALL":
					#sys.stdout.write("#Overall\t%s\n" % (float(nonCR)))
					statistics["Overall"][name] = float(nonCR)
			else:
				stat,value = line
		else:
			#sys.stdout.write("#%d of %d genes have CT coverage >= %d\n" % (len(overall),n,options.CT_number))
			statistics["Genes passed"][name] = len(overall)
			statistics["Genes"][name] = n
			nonCRs["Median"] = np.median(overall)
			statistics["Median"][name] = np.median(overall)
			nonCRs["Mean"] = np.mean(overall)
			statistics["Mean"][name] = np.mean(overall)
			nonCRs["discard"] = 1.0
			#sys.stdout.write("#Median\t%.8f\n" % nonCRs["Median"])
			#sys.stdout.write("#Mean\t%.8f\n" % nonCRs["Mean"])
			#sys.stdout.write("#90%%\t%.8f\n" % np.percentile(overall,90))
			statistics["90%"][name] = np.percentile(overall,90)
			#sys.stdout.write("#75%%\t%.8f\n" % np.percentile(overall,75))
			statistics["75%"][name] = np.percentile(overall,75)
			#sys.stdout.write("#50%%\t%.8f\n" % np.percentile(overall,50))
			statistics["50%"][name] = np.percentile(overall,50)
			#sys.stdout.write("#25%%\t%.8f\n" % np.percentile(overall,25))
			statistics["25%"][name] = np.percentile(overall,25)
			#sys.stdout.write("#10%%\t%.8f\n" % np.percentile(overall,10))
			statistics["10%"][name] = np.percentile(overall,10)
			break
		line = input.readline() #Note that the first line of data was read
	control_nonCR = []
	if control:
		for key,value in control.items():
			if value != "None":
				control_nonCR.append(value)
		nonCRs["control"] = np.median(control_nonCR)
		statistics["control median"][name] = np.mean(control_nonCR)
		statistics["control mean"][name] = np.median(control_nonCR)
		try:
			statistics["control overall"][name] = control_C/(control_all + 0.0)
		except ZeroDivisionError:
			statistics["control overall"][name] = np.nan
		statistics["control used"][name] = len(control_nonCR)
		statistics["control total"][name] = len(control)
		#sys.stdout.write("Control number: %d of %d; median: %.6f; mean: %.6f\n" \		% (len(control_nonCR),len(control),np.median(control_nonCR),np.mean(control_nonCR)))
	return line,input


def read_fn_list():
	search_list = []
	fetch_list = []
	
	with open(options.input,'r') as input:
		line = input.readline()
		if line.startswith("#"):
			if not options.name:
				search_list.append([options.input,options.input])
			else:
				search_list.append([options.name,options.input])
		else:
			while (line):
				line = line.strip().split("\t") #name	directory	CR mode	C-cutoff
				name,dir = line[0],line[1]
				search_list.append([name,dir])
				line = input.readline()
			
	return search_list
	
if __name__ == "__main__":
	description = """
	"""
	parser = argparse.ArgumentParser(prog="m5C_caller_multiple",fromfile_prefix_chars='@',description=description,formatter_class=argparse.RawTextHelpFormatter)
	#Require
	group_required = parser.add_argument_group("Required")
	group_required.add_argument("-i","--input",dest="input",required=True,help="a pileup file or a list. If a list is given format as: name\tdir\tCR\tC-cutoff")
	#Running
	group_optional = parser.add_argument_group("Optional")
	group_optional.add_argument("-n","--name",dest="name",help="name to use, if single file as input")
	group_optional.add_argument("--control",dest="control",help="control list, median non-conversion rate will be used")
	group_optional.add_argument("-o","--output",dest="output",required=False,default="multiple.csv",help="output name, default=multiple.csv")
	#Filter
	group_site = parser.add_argument_group("m5C filter")
	group_site.add_argument("-g","--gene-CR",dest="gene_CR",default=0.05,type=float,help="conversion rate, below which a gene will be discarded, default=0.05")
	group_site.add_argument("-N","--CT",dest="CT_number",default=1000,type=int,help="CT count, below which a gene will be discarded, default=1000")
	
	#Statistics
	group_stat = parser.add_argument_group("Statistic method")
	#Version
	group_other = parser.add_argument_group("Other")
	group_other.add_argument("--version",action="version",version="%(prog)s 1.0")
	options = parser.parse_args()
	
	# sys.stderr.write("[%s] CMD: %s\n" % (strftime("%Y-%m-%d %H:%M:%S", time.localtime())," ".join(sys.argv)))

	C_cutoff_cols = {}
	col_start = 14
	n = 0
	for i in ['1','2','3','4','5','6','7','8','9','10','15','20','None'][:-1]:
		C_cutoff_cols[i] = 14 + n
		n += 1
	C_cutoff_cols["None"] = -1
	
	control = {}
	if options.control:
		with open(options.control,'r') as input:
			for line in input.readlines():
				control[line.strip()] = "None"
	
	search_list = read_fn_list()
	
	statistics = defaultdict(dict)
	#name	directory	CR mode	C-cutoff
	for name,directory in search_list:
		nonCRs = {}
		#sys.stdout.write("Name: %s \n" % name)
		with open(directory,'r') as input:
			read_CR(name,input,nonCRs,control)
	df = pd.DataFrame(statistics)
	pd.options.display.max_rows = 999
	df = df[df.columns[::-1]]
	print df
	if options.output:
		df.to_csv(options.output,sep="\t")