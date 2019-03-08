#!/usr/bin/env python
import os,sys
import argparse
from collections import defaultdict
import pandas as pd
import numpy as np
import scipy.stats
import time
from time import gmtime, strftime
import scipy.stats

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

def filter(row):
	cov,C_count,pvalue,signal,nonCR = row["coverage"],row["C count"],row["P-value"],row["signal"],row["nonCR"]
	# print cov,C_count,pvalue,signal,nonCR
	if cov >= options.coverage and C_count >= options.count and pvalue < options.pvalue and signal >= options.signal and nonCR < options.gene_CR:
		return True
	else:
		return False
	

def handling_groups(df):
	samples = set(df.columns.get_level_values(0))
	for sample in samples:
		df[(sample,"passed")] = False
		subdf = df.xs(sample,axis=1,level=0)
		subdf = subdf[(subdf["coverage"] >= options.coverage ) & (subdf["C count"]>= options.count) & (subdf["P-value"]< options.pvalue) & (subdf["signal"]>= options.signal) & (subdf["nonCR"] < options.gene_CR) & (subdf["m5C level"] >= options.ratio)]
		# df[(sample,"passed")] = df.apply(lambda row: filter(df.xs(sample,axis=1,level=0).loc[tuple(row.name)]), axis=1)
		df.loc[subdf.index,(sample,"passed")] = True
		if options.no_mask_low == False:
			#low coverage
			indexes = df[df[(sample,"coverage")]<options.coverage_nan].index #|(df[(sample,"nonCR")]>=options.gene_CR)|(df[(sample,"signal")]<options.signal)
			df.loc[indexes,(sample,"coverage")] = np.nan
			df.loc[indexes,(sample,"m5C level")] = np.nan
			df.loc[indexes,(sample,"C count")] = np.nan
			df.loc[indexes,(sample,"P-value")] = np.nan
			df.loc[indexes,(sample,"signal")] = np.nan
			df.loc[indexes,(sample,"nonCR")] = np.nan
		if options.no_mask_bad == False:
			# bad conversion
			indexes = df[(df[(sample,"nonCR")]>=options.gene_CR)|(df[(sample,"signal")]<options.signal)].index
			df.loc[indexes,(sample,"coverage")] = np.nan
			df.loc[indexes,(sample,"m5C level")] = np.nan
			df.loc[indexes,(sample,"C count")] = np.nan
			df.loc[indexes,(sample,"P-value")] = 1
			df.loc[indexes,(sample,"signal")] = np.nan
			df.loc[indexes,(sample,"nonCR")] = np.nan
			
	order = {'coverage':0,'C count':1,'m5C level':2,'passed':3,'P-value':4,'signal':5,'nonCR':6}
	columns = sorted(df.columns.tolist(), key = lambda x:(x[0],order[x[1]]))
	df = df[columns]
	if options.no_discard == False:
		df = df[df.xs("passed",axis=1,level=1).any(axis=1)==True] #discard all False sites
		
	# if options.discard_na == True:
		# indexes = df.index.tolist()
		# index_used = []
		# for idx in indexes:
			# x,y,z = idx
			# print x,y,z 
			# print type(y)
			# if y is np.nan == False and  z is np.nan == False:
				# index_used.append(idx)
		# df = df.loc[index_used]
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

if __name__ == "__main__":
	description = """
	"""
	parser = argparse.ArgumentParser(prog="m5C_caller_multiple",fromfile_prefix_chars='@',description=description,formatter_class=argparse.RawTextHelpFormatter)
	#Require
	group_required = parser.add_argument_group("Required")
	group_required.add_argument("-i","--input",dest="input",required=True,help="CSV tables, seperated by comma")
	# group_required.add_argument("-l","--list",dest="list",required=True,help="a list. Format as (TAB seperated): group name order")
	group_required.add_argument("-o","--output",dest="output",required=False,default="intersection.csv",help="output name, default=intersection.csv")
	#Single sample
	group_sample = parser.add_argument_group("Replice inner filter")
	group_sample.add_argument("-c","--coverage",dest="coverage",default=20,type=int,help="minimal coverage for each site, default=20")
	group_sample.add_argument("-C","--count",dest="count",default=3,type=int,help="minimal C count for each site, default=3")
	group_sample.add_argument("-r","--ratio",dest="ratio",default=0.1,type=float,help="minimal ratio for each site, default=0.0")
	group_sample.add_argument("-p","--pvalue",dest="pvalue",default=0.001,type=float,help="pvalue, default=0.001")
	group_sample.add_argument("--na-coverage",dest="coverage_nan",default=10,type=int,help="below this value, set coverage as nan, default=10")
	group_sample.add_argument("-s","--signal",dest="signal",default=0.9,type=float,help="signal ratio, equals coverage(under C-cutoff)/coverage, default=0.9")
	group_sample.add_argument("-g","--gene-CR",dest="gene_CR",default=0.05,type=float,help="conversion rate, below which a gene will be discarded, default=0.05")
	group_sample.add_argument("--no-discard",dest="no_discard",default=False,action="store_true",help="Do not discard all False sites")
	group_sample.add_argument("--no-mask-low",dest="no_mask_low",default=False,action="store_true",help="Do not mask low coverage sites")
	group_sample.add_argument("--no-mask-bad",dest="no_mask_bad",default=False,action="store_true",help="Do not mask hard conversion sites")
	# group_sample.add_argument("--discard-na",dest="discard_na",default=False,action="store_true",help="Discard sites without annotation")
	#Version
	group_other = parser.add_argument_group("Other")
	group_other.add_argument("--version",action="version",version="%(prog)s 1.0")
	options = parser.parse_args()
	
	sys.stderr.write("[%s] CMD: %s\n" % (strftime("%Y-%m-%d %H:%M:%S", time.localtime())," ".join(sys.argv)))
	
	# sample_groups = handling_list()
	df = handling_inputs()
	handling_groups(df)
	
	sys.stderr.write("[%s] Finished.\n" % strftime("%Y-%m-%d %H:%M:%S"))