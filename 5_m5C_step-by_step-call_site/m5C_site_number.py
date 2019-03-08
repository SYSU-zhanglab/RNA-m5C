#!bin/usr/env python
import sys
import os
import argparse
import pandas as pd
import numpy as np
from collections import defaultdict

if __name__ == "__main__":
	description = """
	"""
	parser = argparse.ArgumentParser(prog="m5C_site_number",fromfile_prefix_chars='@',description=description,formatter_class=argparse.RawTextHelpFormatter)
	#Require
	group_required = parser.add_argument_group("Required")
	group_required.add_argument("-i","--input",dest="input",required=True,help="input table")	#Version
	group_other = parser.add_argument_group("Other")
	group_other.add_argument("--version",action="version",version="%(prog)s 1.0")
	options = parser.parse_args()
	
	df = pd.read_csv(options.input,header=[0,1],index_col=[0,1,2])
	
	samples = set(df.columns.get_level_values(0))
	
	for sample in samples:
		subdf = df.xs(sample,axis=1,level=0)
		subdf = subdf[(subdf["m5C level"] >= 0.1) & (subdf["C count"] >= 3) & (subdf["P-value"] < 0.001) & (subdf["coverage"] >= 20) & (subdf["signal"] >= 0.9) & (subdf["nonCR"] < 0.05)]
		indexes = subdf.index.tolist()
		indexes = [i for i in indexes if "ERCC" not in i[0] and "SIRV" not in i[0] and "vitro" not in i [0]]
		print sample,len(indexes)