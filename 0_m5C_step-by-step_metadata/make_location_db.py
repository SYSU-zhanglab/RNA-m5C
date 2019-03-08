#!bin/usr/env python

#Jianheng Liu @ Zhanglab, SYSU
#Feb, 2018
#Email: liujh26@mail2.sysu.edu.cn
#Usage: Turn single base annotation to Sqlite3 database. Redundance should be removed.
#Input: [.base]

import sqlite3
import sys
import time
from time import gmtime, strftime
import argparse

def create_connection(db_file):
	conn = sqlite3.connect(db_file)
	return conn

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

def create_index(cursor):
	cursor.execute("CREATE INDEX locations_index ON locations (chr,dir,pos)")

description = """
Create sqlite database
columns:
chr,pos_1,dir,gene,trans,name,isoform,biotype
index: chr,pos_1,dir
"""
parser = argparse.ArgumentParser(prog="make_location_db",fromfile_prefix_chars='@',description=description,formatter_class=argparse.RawTextHelpFormatter)
#Require
group_require = parser.add_argument_group("Required")
group_require.add_argument("-i","--input",dest="input",required=True,help="input file")
group_require.add_argument("-o","--output",dest="output",required=True,help="database output") #chr start_0 start
group_other = parser.add_argument_group("Other")
options = parser.parse_args()

database = create_connection(options.output)
sys.stderr.write("[%s] Connection created\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
cursor = database.cursor()
cursor.execute("PRAGMA cache_size =8000")
cursor.execute("PRAGMA page_size = 8192")
cursor.execute("PRAGMA journal_mode = WAL") #do not need a write lock while reading, so enable multiple reading

create_table(cursor)
sys.stderr.write("[%s] Begin inserting data...\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
with open(options.input,'r') as input:
	line = input.readline()
	n = 0
	while(line):
		line = line.strip().split("\t")
		chr = line[0]
		pos_1 = line[2]
		dir = line[3]
		gene = line[4]
		trans = line[5]
		name = line[6]
		isoform = line[7]
		biotype = line[8]
		line = input.readline()
		#id = "@".join([chr,pos_1,dir])
		insert_var(cursor,chr,int(pos_1),dir,gene,trans,name,isoform,biotype)
		n += 1
		if n%1000000 == 0:
			sys.stderr.write("[%s] %d items processed...\n" % (strftime("%Y-%m-%d %H:%M:%S", time.localtime()),n))
sys.stderr.write("[%s] All loaded. Creating index...\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
create_index(cursor)
sys.stderr.write("[%s] All finished!\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
database.close()

