#!/usr/bin/env python
import os,sys
import argparse
from collections import defaultdict,Counter
import copy
import sqlite3
import time
from time import gmtime, strftime
import numpy as np
import gzip
import csv

def create_table(cursor): #id text PRIMARY KEY
	sql = """ CREATE TABLE IF NOT EXISTS 
			locations (
				chr text NOT NULL,
				pos integer NOT NULL,
				dir text NOT NULL,
				gene text NOT NULL,
				trans text NOT NULL,
				name text NOT NULL,
				isoform text NOT NULL,
				biotype text NOT NULL
			); 
		   """
	cursor.execute(sql)

def insert_var(cursor,chr,pos,dir,gene,trans,name,isoform,biotype):
	cursor.execute("insert into locations (chr, pos, dir, gene, trans, name, isoform, biotype) values ( ?, ?, ?, ?, ?, ?, ?, ?)", (chr,pos,dir,gene,trans,name,isoform,biotype)) #id

def insert_many_var(cursor,rows):
	cursor.executemany("""
	insert into locations (chr, pos, dir, gene, trans, name, isoform, biotype) values ( ?, ?, ?, ?, ?, ?, ?, ?)
	""",rows)

def create_index(cursor):
	cursor.execute("CREATE INDEX locations_index ON locations (chr,dir,pos)")

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

def build_database():
	global database,cursor,conn
	# with open(options.database,'r') as input:
	with open_database(options.database) as input:
		if options.method == "sql":
			# conn = sqlite3.connect(":memory:")
			conn = sqlite3.connect("file::memory:?cache=shared")  # for python 2.7.16, add uri=True for python >=3.6
			cursor = conn.cursor()
			create_table(cursor)
			cursor.execute("PRAGMA cache_size=65536")
			cursor.execute("PRAGMA page_size=65536")
			cursor.execute('PRAGMA temp_store=MEMORY')
			cursor.execute("PRAGMA synchronous=OFF")
			cursor.execute('PRAGMA journal_mode=MEMORY')
		elif options.method == "dict":
			# cursor = None
			# database = {}
			pass
		n = 0
		rows = []
		line = input.readline()
		while(line):
			# line = line.strip().split("\t")
			chr = line[0]
			pos_1 = line[2]
			dir = line[3]
			gene = line[4]
			trans = line[5]
			name = line[6]
			isoform = line[7]
			biotype = line[8]
			#id = "@".join([chr,pos_1,dir])
			if options.method == "sql":
				#slow
				#insert_var(cursor,chr,int(pos_1),dir,gene,trans,name,isoform,biotype)
				rows.append([chr,int(pos_1),dir,gene,trans,name,isoform,biotype])
			elif options.method == "dict":
				database[(chr,int(pos_1),dir)] = [gene,name,trans,isoform,biotype] #SELECT gene,name,trans,isoform,biotype
			n += 1
			if n%1000000 == 0:
				if options.method == "sql":
					insert_many_var(cursor,rows)
					rows = []
				sys.stderr.write("[%s] %d items processed...\n" % (strftime("%Y-%m-%d %H:%M:%S", time.localtime()),n))
			line = input.readline()
		if options.method == "sql":
			if rows:
				insert_many_var(cursor,rows)
				rows = []
	if options.method == "sql":
		sys.stderr.write("[%s] All loaded. Creating index...\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
		create_index(cursor)
		sys.stderr.write("[%s] All finished!\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
		# return database,cursor
	elif options.method == "dict":
		pass
		# return database,cursor

def write_CR():
	#ALL	900365017	2875105	0.00319326600403
	#ELSE	85530991	352474	0.00412100919069
	#ENSG00000000003	25890	148	0.00571649285438
	overall = []
	lines = []
	with open(options.conversion,'w') as CR_output:
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
'''
def create_connection(db_file):
	if os.path.isfile(db_file):
		conn_old = sqlite3.connect(db_file)
		conn = sqlite3.connect(':memory:')
		query = ""
		for line in conn_old.iterdump():
			if line.endswith(";")
				conn.executescript(query)
				query = ""
			else:
				query = "".join([query,line])
		# query = "".join(line for line in conn_old.iterdump())
		# conn.executescript(query)
		conn_old.close()
	else:
		raise IOError("Database file not exist.\n")
	return conn
'''
def search(chr,dir,pos,database=None,cursor=None):
	if options.method == "sql":
		cursor.execute("SELECT gene,name,trans,isoform,biotype FROM locations WHERE chr=? AND dir=? AND pos=?",(chr,dir,pos))
		rows = cursor.fetchall()
		if rows:
			return [str(i) for i in rows[0]]
		else:
			return None
	else:
		return database.get((chr,dir,pos))
	
def Tideup(chr,pos_1,dir,PIPLEUPS,SURROUNDINGS,database=None,cursor=None):
	global output,all_C,all_cov
	#print 2,len(PIPLEUPS),len(SURROUNDINGS)
	if options.method == "sql":
		query = search(chr,dir,pos_1,cursor=cursor)
	elif options.method=="dict":
		query = search(chr,dir,pos_1,database=database)
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
			counts[base] += count
			counts["total"] += count
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
	
def check(chr,pos_1,dir,pileup,surrounding,database,cursor):
	PIPLEUPS = pileup
	SURROUNDINGS = surrounding
	Tideup(chr,pos_1,dir,PIPLEUPS,SURROUNDINGS,database=database,cursor=cursor)

if __name__ == "__main__":
	description = """
	Input should be no redundance
	"""
	parser = argparse.ArgumentParser(prog="",fromfile_prefix_chars='@',description=description,formatter_class=argparse.RawTextHelpFormatter)
	#Require
	group_require = parser.add_argument_group("Required")
	group_require.add_argument("-i","--input",dest="input",required=True,help="input pileup, sorted")
	group_require.add_argument("-o","--output",dest="output",required=True,help="output") #chr start_0 start
	group_require.add_argument("--CR",dest="conversion",required=True,help="gene conversion rate")
	group_require.add_argument("--db",dest="database",required=True,help="database, base-gene annotation")
	group_require.add_argument("--method",dest="method",default="sql",choices=["sql","dict"],help="database type, sql or dictionary")
	# group_optional = parser.add_argument_group("Optional")
	group_other = parser.add_argument_group("Other")
	options = parser.parse_args()
	
	sys.stderr.write("[%s] Reading input...\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
	
	if options.method == "sql":
		database = None
		cursor = None
		conn = None
	elif options.method == "dict":
		database = {}
		cursor = None
		conn = None
	build_database()
	sys.stderr.write("[%s] In-memory database connection setup\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
	#init
	result_cov = defaultdict(int)
	result_C = defaultdict(int)
	all_cov = 0
	all_C = 0

	with open_database(options.input) as input,open(options.output+"_tmp",'w') as output:  # open(options.input,'r') as input
		#10      100007581(1-based)       -       C       ENST00000260702 3519    T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
		#10      100007581       -       C       Genome  100007581       T,T,T,T,T,T,T,T,T,T,T   0,0,0,0,0,0,0,0,0,0,0
		line = input.readline()
		while (line):
			row = line # .strip().split("\t")
			chr = row[0]
			pos_1 = row[1]
			dir = row[2]
			base = row[3] # as transcript
			pileup = row[6].split(",")
			surrounding = [int(i) for i in row[7].split(",")]
			if base == "C":
				check(chr,pos_1,dir,pileup,surrounding,database,cursor)
				# if options.method == "sql":
					# check(chr,pos_1,dir,pileup,surrounding,database=database,cursor=cursor)
				# elif options.method == "dict":
					# check(chr,pos_1,dir,pileup,surrounding,database=database,cursor=cursor)
			line = input.readline()

	sys.stderr.write("[%s] Finished reading input.\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
	#get conversion rates
	write_CR()
	sys.stderr.write("[%s] Calculating conversion rates...\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
	
	# for better performance
	with open(options.output, "w") as output:
		with open(options.conversion, "r") as input:
			line = input.readline()
			while (line):
				output.write(line)
				line = input.readline()
		with open(options.output+"_tmp", "r") as input:
			line = input.readline()
			while (line):
				output.write(line)
				line = input.readline()
				
	# huge memory usage, deprecated
	# os.system("cat %s %s > %s" % (options.conversion,options.output+"_tmp",options.output))
	
	os.remove(options.output+"_tmp")
	if options.database == "sql":
		database.close()
	else:
		del database
	sys.stderr.write("[%s] All finished!\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
