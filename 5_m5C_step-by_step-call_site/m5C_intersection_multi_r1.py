#!bin/usr/env python

#Jianheng Liu @ Zhanglab, SYSU
#Feb, 2018
#Email: liujh26@mail2.sysu.edu.cn
#Usage: This program is used to intersect m5C candidates in a candidate marked csv
#Input: [.csv]

import os,sys
import argparse
from collections import defaultdict
from collections import OrderedDict
import pandas as pd
import numpy as np
import scipy.stats
import time
from time import gmtime, strftime
import scipy.stats

def P_values_combine(pvalues,combine = True):
	pvalue_used = []
	if combine == True:
		for p in pvalues:
			if np.isnan(p) == False:
				if p != 0.0:
					pvalue_used.append(p)
				else:
					pvalue_used.append(1E-300)
		if len(pvalue_used) > 1:
			return scipy.stats.combine_pvalues(pvalue_used,method=options.p_combine_method)[1]
		elif len(pvalue_used) == 1:
			return pvalue_used[0]
		else:
			return np.nan
	else:
		return pvalues[0]
		
		
def count_replicates(pvalues,combine = True):
	pvalue_used = []
	if combine == True:
		for p in pvalues:
			if np.isnan(p) == False:
				if p != 0.0:
					pvalue_used.append(p)
				else:
					pvalue_used.append(1E-300)
		if len(pvalue_used) >= 1:
			return len(pvalue_used)
		else:
			return 0
	else:
		return len(pvalues)
		
def count_present(pvalues):
	a = 0
	for p in pvalues:
		if np.isnan(p) == False:
			a += 1
	return a

def subdf_formation(subs,format_df):
	if len(subs) == 1:
		min_C = options.single_count
	else:
		min_C = options.count
	for item in subs:
		not_m5C_indexes = format_df[(format_df[(item,"coverage")] >= options.coverage) & (format_df[(item,"C count")] < min_C)].index
		format_df.loc[not_m5C_indexes,(item,"P-value")] = 1.0
		#bad sites, all info = 0
		low_cover_indexes = format_df[format_df[(item,"coverage")] < options.coverage_nan].index
		format_df.loc[low_cover_indexes,(item,"m5C level")] = np.nan
		format_df.loc[low_cover_indexes,(item,"coverage")] = 0
		format_df.loc[low_cover_indexes,(item,"C count")] = 0
		format_df.loc[low_cover_indexes,(item,"coverage")] = 0
		format_df.loc[low_cover_indexes,(item,"P-value")] = np.nan
		
		bad_indexes = format_df[(format_df[(item,"signal")] < options.signal) | (format_df[(item,"nonCR")]>=options.gene_CR)].index
		format_df.loc[bad_indexes,(item,"m5C level")] = np.nan
		format_df.loc[bad_indexes,(item,"coverage")] = 0
		format_df.loc[bad_indexes,(item,"C count")] = 0
		format_df.loc[bad_indexes,(item,"coverage")] = 0
		format_df.loc[bad_indexes,(item,"P-value")] = 1.0

	return format_df

def handling_groups(df,sample_groups):
	dfs = []
	
	# Reorder indexes
	# coverage	C count	m5C level	P-combined	passed	replicates
	reorder_columns = []

	for group, subs in sample_groups.iteritems():
		subdf = df.iloc[:, df.columns.get_level_values(0).isin(subs)].copy()
		
		reorder_columns.append((group,"coverage"))
		reorder_columns.append((group,"C count"))
		reorder_columns.append((group,"m5C level"))
		reorder_columns.append((group,"P-combined"))
		reorder_columns.append((group,"passed"))
		reorder_columns.append((group,"replicates"))
		
		# for sub in subs:
			# indexes = subdf[(subdf[(sub,"coverage")]>=options.coverage_nan)&(subdf[(sub,"C count")]==0)].index
			# subdf.loc[indexes,(sub,"m5C level")] = np.nan
			# subdf.loc[indexes,(sub,"C count")] = np.nan
			# subdf.loc[indexes,(sub,"coverage")] = np.nan
			# subdf.loc[indexes,(sub,"")] = 1
			
			# indexes = subdf[(subdf[(sub,"coverage")]<options.coverage)&(subdf[(sub,"P-value")]!=1)].index
			# subdf.loc[indexes,(sub,"m5C level")] = np.nan
			# subdf.loc[indexes,(sub,"C count")] = np.nan
			# subdf.loc[indexes,(sub,"coverage")] = np.nan
			# subdf.loc[indexes,(sub,"P-value")] = np.nan
			
			# indexes = subdf[(subdf[(sub,"coverage")]<options.coverage)&(subdf[(sub,"P-value")]!=1)].index
		
		# format_df = subdf_formation(subs,subdf)
		
		indexes = subdf[subdf.xs("passed",axis=1,level=1).any(axis=1)==True].index
		
		df_new = pd.DataFrame(columns=pd.MultiIndex.from_product([[group],["coverage","C count","m5C level","P-combined","passed","replicates"]]))
		
		#Coverages, C count and m5C level
		df_new[(group,"coverage")] = subdf.xs("coverage",axis=1,level=1).sum(axis=1)
		df_new[(group,"C count")] = subdf.xs("C count",axis=1,level=1).sum(axis=1)
		df_new[(group,"m5C level")] = df_new[(group,"C count")]/df_new[(group,"coverage")]
		
		#P value
		df_new[(group,"P-combined")] = 1.0
		df_new.loc[indexes,(group,"P-combined")] =  df_new.apply(lambda row: P_values_combine(subdf.xs("P-value",level=1,axis=1).loc[tuple(row.name)].tolist(),combine = True), axis=1)
		
		#Replicates
		df_new[(group,"replicates")] = 0
		# df_new[(group,"present")] = 0
		
		df_new[(group,"replicates")] = subdf.xs("coverage",axis=1,level=1).notnull().sum(axis=1)
		# print subdf.xs("coverage",axis=1,level=1).isnull().sum(axis=1)
		#df_new.apply(lambda row: count_replicates(format_df.xs("P-value",level=1,axis=1).loc[tuple(row.name)].tolist(),combine = True), axis=1)
		# df_new[(group,"present")] = df_new.apply(lambda row: count_present(format_df.xs("P-value",level=1,axis=1).loc[tuple(row.name)].tolist()), axis=1)
		
		#passed
		df_new[(group,"passed")] = False
		# df_new.loc[df_new[(group,"coverage")]<options.coverage_nan,(group,"m5C level")] = np.nan
		df_new.loc[(df_new[(group,"P-combined")]<options.p_combine)&(df_new[(group,"coverage")]>=options.coverage),(group,"passed")] = True
		df_new.loc[(df_new[(group,"P-combined")]>=options.p_combine),(group,"replicates")] = 0
		df_new.loc[(df_new[(group,"m5C level")]<options.merged_ratio),(group,"passed")] = False
		df_new.loc[(df_new[(group,"m5C level")]<options.merged_ratio),(group,"P-combined")] = 1
		#drop not enough replicates
		sub_len = len(subs)
		
		# df_new.loc[df_new[(group,"C count")]<options.count*df_new[(group,"present")],(group,"passed")] = False
		df_new.loc[(df_new[(group,"passed")]==True) & ((df_new[(group,"replicates")] == 1))& (df_new[(group,"C count")] < options.single_count),(group,"passed")] = False
		
		dfs.append(df_new)
		
		
	df = pd.concat(dfs,axis=1)
	if options.no_discard == False:
		df = df[df.xs("passed",axis=1,level=1).any(axis=1)==True] #discard all False sites
	
	# reorder columns
	df = df[reorder_columns]
	df.to_csv(options.output)
		
def handling_inputs():
	files = options.input.split(",")
	if len(files) == 1:
		df = pd.read_csv(files[0],index_col=[0,1,2],header=[0,1])
		return df
	else:
		dfs = {}
		for df_name in files:
			dfs[df_name] = pd.read_csv(df_name,index_col=[0,1,2],header=[0,1])
		df = pd.concat(dfs.values(),axis=1)
		del dfs
		return df
		
def handling_list():
	sample_groups = OrderedDict()
	with open(options.list,'r') as input:
		for line in input.readlines():
			line = line.strip().split("\t")
			group = line[0]
			name = line[1]
			# if len(line) > 2:
				# order = int(line[2])
			# else:
				# order = 999
			if group not in sample_groups:
				sample_groups[group] = []
			sample_groups[group].append(name)#(name,order))
	return sample_groups

if __name__ == "__main__":
	description = """
	Intersect replicates by a given list.
	Note: 
	(1) when computing combined P-value, zero P-values will be treated as 1E-300
	(2) singleton samples were not combined when mixed with paired samples
	"""
	parser = argparse.ArgumentParser(prog="m5C_caller_multiple",fromfile_prefix_chars='@',description=description,formatter_class=argparse.RawTextHelpFormatter)
	#Require
	group_required = parser.add_argument_group("Required")
	group_required.add_argument("-i","--input",dest="input",required=True,help="CSV tables, seperated by comma")
	group_required.add_argument("-l","--list",dest="list",required=True,help="a list. Format as (TAB seperated): group name order")
	group_required.add_argument("-o","--output",dest="output",required=False,default="intersection.csv",help="output name, default=intersection.csv")
	#Single sample
	group_sample = parser.add_argument_group("Replice inner filter")
	group_sample.add_argument("-c","--coverage",dest="coverage",default=20,type=int,help="minimal coverage for each site, default=20")
	group_sample.add_argument("-C","--count",dest="count",default=3,type=int,help="minimal C count for each site, default=3")
	group_sample.add_argument("--single-C",dest="single_count",default=5,type=int,help="minimal C count for sites without replicates, default=5")
	# group_sample.add_argument("-r","--ratio",dest="ratio",default=0.0,type=float,help="minimal ratio for each site, default=0.0")
	group_sample.add_argument("-p","--pvalue",dest="pvalue",default=0.001,type=float,help="pvalue, default=0.001")
	group_sample.add_argument("--na-coverage",dest="coverage_nan",default=10,type=int,help="below this value, set coverage as nan, default=10")
	group_sample.add_argument("-s","--signal",dest="signal",default=0.9,type=float,help="signal ratio, equals coverage(under C-cutoff)/coverage, default=0.9")
	group_sample.add_argument("-g","--gene-CR",dest="gene_CR",default=0.05,type=float,help="conversion rate, below which a gene will be discarded, default=0.05")
	#Merged sample
	group_merged = parser.add_argument_group("Merged sample filter")
	group_merged.add_argument("-P","--pcomb",dest="p_combine",default=0.001,type=float,help="pvalue, default=0.001")
	group_merged.add_argument("--P-method",dest="p_combine_method",default="stouffer",choices=["stouffer","fisher"],help="method used in combining pvalues, see scipy.stats.combine_pvalues")
	# group_merged.add_argument("-a","--all",dest="all_present",default=False,action="store_true",help="A site should present in all detected samples")
	# group_merged.add_argument("-n","--replicate",dest="replicate",default=1,type=int,help="Minmimal replicates passed filter")
	# group_merged.add_argument("--merged-coverage",dest="merged_coverage",type=int,help="minimal coverage for merged site, if not given, equals minimal replicate * minimal coverage")
	# group_merged.add_argument("--merged-count",dest="merged_count",type=float,help="minimal ratio for merged site, if not given, equals minimal replicate * minimal C count")
	group_merged.add_argument("-R","--merged-ratio",dest="merged_ratio",default=0.1,type=float,help="minimal ratio for merged site, default=0.1")
	group_merged.add_argument("--no-discard",dest="no_discard",default=False,action="store_true",help="Do not discard all False sites")
	#Version
	group_other = parser.add_argument_group("Other")
	group_other.add_argument("--version",action="version",version="%(prog)s 1.0")
	options = parser.parse_args()
	
	sys.stderr.write("[%s] CMD: %s\n" % (strftime("%Y-%m-%d %H:%M:%S", time.localtime())," ".join(sys.argv)))
	
	sample_groups = handling_list()
	df = handling_inputs()
	handling_groups(df,sample_groups)
	
	sys.stderr.write("[%s] Finished.\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))