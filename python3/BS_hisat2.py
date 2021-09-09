#!bin/usr/env python

#Jianheng Liu @ Zhanglab, SYSU
#Feb, 2018
#Email: liujh26@mail2.sysu.edu.cn
#Usage: Map reads to HISAT2 indexed genome
#Input: [.fastq] [HISAT2 index]

__doc__ = """
Usage: Map reads to HISAT2 indexed genome.
Note that the input should be stranded. If the conversion
rates are very low, this program might take too much RAM.
In such a situation, please split the fastq input and merge
BAMs later.
"""


import os,sys,time,signal,argparse,pysam,multiprocessing
from Bio import SeqIO
from Bio.Seq import reverse_complement
from time import gmtime, strftime
from collections import defaultdict
from pysam import qualities_to_qualitystring
import subprocess


def read_parameters(para_file):
	parameters = {}
	if para_file:
		with open(para_file, 'r') as input:
			for line in input.readlines():
				if line.startswith("#") == False:
					line = line.split(" ")
					parameters[line[0]] = line[1]
	return parameters

def fastq_converter_1(fin,fout,method):
	sys.stderr.write("[%s]Converting %s, %s...\n" % (strftime("%Y-%m-%d %H:%M:%S", time.localtime()),fin,method))
	total_reads = 0
	output_dict = {}
	if method == "C2T":
		original = "C"
		convert = "T"
	elif method == "G2A":
		original = "G"
		convert = "A"
	with open(fin,'r') as input, open(fout,'w') as output:
		line = input.readline()
		i = 0
		while (line):
			if i == 0:
				output.write(line)
				name = line.strip("@").strip("\n").split( )[0]
				i += 1
			elif i == 1:
				line = line.strip()
				if original in line:
					output_dict[name] = line
					output.write(line.replace(original,convert))
					output.write("\n")
				else:
					output.write(line)
					output.write("\n")
				i += 1
			elif i == 2:
				output.write(line)
				i +=1
			elif i == 3:
				output.write(line)
				i = 0
				total_reads += 1
			line = input.readline()
	return output_dict,total_reads
	
def read_C_position(fwd,rev=None):
	total_reads = 0
	fwd_read_dict = {}
	rev_read_dict = {}
	if rev is not None:
		with open(fwd,'r') as fastq_fwd, open(rev,'r') as fastq_rev:
			line_fwd = fastq_fwd.readline()
			line_rev = fastq_rev.readline()
			n = 0
			while (line_fwd):
				if n == 0:
					seq_name = line_fwd.strip("@").strip("\n").split( )[0]
					n=n+1
					total_reads += 1
				elif n == 1:
					fwd_C = line_fwd.strip("\n").count("C")
					rev_G = line_rev.strip("\n").count("G")
					if fwd_C != 0:
						fwd_read_dict[seq_name] = line_fwd.strip("\n")
					if rev_G != 0:
						rev_read_dict[seq_name] = line_rev.strip("\n")
					n=n+1
				elif n == 2:
					n=n+1
				elif n == 3:
					n=0
				line_fwd = fastq_fwd.readline()
				line_rev = fastq_rev.readline()
	else:
		with open(fwd,'r') as fastq_fwd:
			line_fwd = fastq_fwd.readline()
			n = 0
			while (line_fwd):
				if n == 0:
					seq_name = line_fwd.strip("@").strip("\n").split( )[0]
					n=n+1
					total_reads += 1
				elif n == 1:
					fwd_C = line_fwd.strip("\n").count("C")
					if fwd_C != 0:
						fwd_read_dict[seq_name] = line_fwd.strip("\n")
					n=n+1
				elif n == 2:
					n=n+1
				elif n == 3:
					n=0
				line_fwd = fastq_fwd.readline()
	return fwd_read_dict,rev_read_dict,total_reads
	
def mapping_string(parameters,fwd,rev):
	if "-p" not in parameters:
		parameters["-p"] = '5'
	if "-k" not in parameters:
		parameters["-k"] = '10'
	if "--reorder" not in parameters:
		parameters["--reorder"] = None
	if fwd and rev:
		parameters["-1"] = fwd
		parameters["-2"] = rev
		if "--rna-strandness" not in parameters:
			parameters["--rna-strandness"] = "FR"
		if "--no-mixed" not in parameters:
			parameters["--no-mixed"] = None
		if "--fr" not in parameters:
			parameters["--fr"] = None
	elif fwd and not rev:
		parameters["-U"] = fwd
		if "--rna-strandness" not in parameters:
			parameters["--rna-strandness"] = "F"
		
	# if options.avoid_psedu and "--avoid-pseudogene" not in parameters:
		# parameters["--avoid-pseudogene"] = None
		
	para = []
	for key,values in parameters.items():
		if values is not None:
			para.extend([key,values])
		else:
			para.extend([key])
	
	para = "|".join(para)
	return para
	
def mapping(index,method,hisat2_path,parameters,preifx):
	if parameters:
		parameters = parameters.split("|")
	else:
		parameters = []
	summary_name = method + "hisat2.summary"

	if method == "C2T":
		sam = preifx + ".C2T.sam"
	elif method == "G2A":
		sam = preifx + ".G2A.sam"
	cmd = [hisat2_path+"/hisat2", \
		    "-x",index,\
			"-S",sam] + \
		   parameters

	hisat2_out = subprocess.Popen(cmd,\
								  stdin=subprocess.PIPE,\
								  stdout=subprocess.PIPE,\
								  stderr=subprocess.PIPE,\
								  shell=False)
	stdout,stderr = hisat2_out.communicate()
	# print stdout
	sys.stderr.write("[%s]%s report:\n" % (strftime("%Y-%m-%d %H:%M:%S", time.localtime()),method))
	sys.stderr.write(stderr)
	sys.stderr.write("\n")
	#stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=False
	# hisat2.wait()
def set_tag(segment_tag,genome):
	segment_tag.set_tag("YG",genome)
	if segment_tag.is_read1:
		segment_tag.set_tag("YR","C2T")
	elif segment_tag.is_read2:
		segment_tag.set_tag("YR","G2A")
	else:
		if segment_tag.is_reverse:
			segment_tag.set_tag("YR","G2A")
		else:
			segment_tag.set_tag("YR","C2T")
	return segment_tag
	
def check_mapping(segment_check):
	if segment_check.get_tag("YG") == "C2T" and segment_check.get_tag("XS") == "+" and segment_check.get_tag("AS") >= options.min_AS:
		if segment_check.is_read1 and segment_check.is_reverse == False:
			return True,segment_check
		elif segment_check.is_read2 and segment_check.is_reverse == True:
			return True,segment_check
		elif segment_check.is_paired == False and segment_check.is_reverse == False:
			return True,segment_check
		else:
			return False,segment_check
	elif segment_check.get_tag("YG") == "G2A" and segment_check.get_tag("XS") == "-" and segment_check.get_tag("AS") >= options.min_AS:
		if segment_check.is_read1 and segment_check.is_reverse == True:
			return True,segment_check
		elif segment_check.is_read2 and segment_check.is_reverse == False:
			return True,segment_check
		elif segment_check.is_paired == False and segment_check.is_reverse == True:
			return True,segment_check
		else:
			return False,segment_check
	else:
		return False,segment_check
		
def check_list(TMP,genome,discordant_C2T,discordant_G2A):
	passed = []
	for segment in TMP:
		segment = set_tag(segment,genome)
		status, segment = check_mapping(segment)
		if status == True:
			passed.append(segment)
		else:
			if segment.get_tag("YG") == "C2T":
				discordant_C2T.append(segment)
			elif segment.get_tag("YG") == "G2A":
				discordant_G2A.append(segment)
	return passed,discordant_C2T,discordant_G2A
	
def to_unmapped(segment,unal_read1,unal_read2):
	if segment.is_read1:
		if segment.query_name in fwd_read_dict:
			seq = fwd_read_dict[segment.query_name]
			qual = qualities_to_qualitystring(segment.query_qualities)
			if segment.is_reverse:
				qual = qual[::-1]
		else:
			seq = segment.query_sequence
			qual = qualities_to_qualitystring(segment.query_qualities)
			if segment.is_reverse:
				seq = reverse_complement(seq)
				qual = qual[::-1]
		unal_read1.write("".join(["@",segment.query_name,"\n",seq,'\n+\n',qual,'\n']))
	elif segment.is_read2:
		if segment.query_name in rev_read_dict:
			seq = rev_read_dict[segment.query_name]
			qual = qualities_to_qualitystring(segment.query_qualities)
			if segment.is_reverse:
				qual = qual[::-1]
		else:
			seq = segment.query_sequence
			qual = qualities_to_qualitystring(segment.query_qualities)
			if segment.is_reverse:
				seq = reverse_complement(seq)
				qual = qual[::-1]
		unal_read2.write("".join(["@",segment.query_name,"\n",seq,'\n+\n',qual,'\n']))
	else: #single end and fully unmapped
		if segment.query_name in fwd_read_dict:
			seq = fwd_read_dict[segment.query_name]
			qual = qualities_to_qualitystring(segment.query_qualities)
			if segment.is_reverse:
				qual = qual[::-1]
		else:
			seq = segment.query_sequence
			qual = qualities_to_qualitystring(segment.query_qualities)
			if segment.is_reverse:
				seq = reverse_complement(seq)
				qual = qual[::-1]
		unal_read1.write("".join(["@",segment.query_name,"\n",seq,'\n+\n',qual,'\n']))

def to_unmapped_discordant(discordant,unal_read1,unal_read2):
	read1 = None
	read2 = None
	discordant_iter = iter(discordant)
	while read1 is None or read2 is None:
		segment = next(discordant_iter)
		if not read1 and not read2:
			if segment.is_read1:
				to_unmapped(segment,unal_read1,unal_read2)
				read1 = 1
			elif segment.is_read2:
				to_unmapped(segment,unal_read1,unal_read2)
				read2 = 1
			else:
				to_unmapped(segment,unal_read1,unal_read2)
				read1 = 1
				read2 = 1
		elif not read1 and read2:
			if segment.is_read1:
				to_unmapped(segment,unal_read1,unal_read2)
				read1 = 1
		elif not read2 and read1:
			if segment.is_read2:
				to_unmapped(segment,unal_read1,unal_read2)
				read2 = 1
		else:
			break
	
def to_bam_genome(segments,output):
	for segment in segments:
		if segment.is_read1:
			seq = fwd_read_dict.get(segment.query_name)
		elif segment.is_read2:
			seq = rev_read_dict.get(segment.query_name)
		else:
			seq = fwd_read_dict.get(segment.query_name)
		if seq is not None:
			if segment.is_reverse:
				seq = reverse_complement(str(seq))
			qual = segment.query_qualities
			mpq = segment.mapping_quality
			segment_output = segment
			segment_output.query_sequence = seq
			segment_output.mapping_quality = mpq
			segment_output.query_qualities = qual
			output.write(segment_output)
		else:
			output.write(segment)
	
def check_pairs(has_read1,has_read2,segment):
	if segment.is_read1:
		has_read1 = True
	elif segment.is_read2:
		has_read2 = True
	return has_read1,has_read2
	
def validate_read_pairs(a,b):
	if (a.is_read1 and b.is_read2) or (a.is_read2 and b.is_read1):
		if a.reference_name == b.reference_name:
			return True
		else:
			return False
	else:
		return False

def remove_high_poly(segment):
	A = 0.0
	T = 0.0
	C = 0.0
	G = 0.0
	for i in segment.query_sequence:
		if i == "A":
			A += 1
		elif i == "T":
			T += 1
		elif i == "C":
			C += 1
		elif i == "G":
			G += 1
	Total = A + T + C +G
	if A/Total >= 0.8 or T/Total >= 0.8 or C/Total >= 0.8 or G/Total >= 0.8:
		return True
	else:
		return False


def mapping_result(unique_C2T_tmp,unique_G2A_tmp,unmapped_C2T,unmapped_G2A,multi_C2T_tmp,multi_G2A_tmp,\
                   output,unal_read1,unal_read2,multi_output,mode):
	global multimapper_num,unmapped_num,C2T_num,G2A_num
	discordant_C2T,discordant_G2A = [],[]
	unique_C2T,discordant_C2T,discordant_G2A = check_list(unique_C2T_tmp,"C2T",discordant_C2T,discordant_G2A)
	unique_G2A,discordant_C2T,discordant_G2A = check_list(unique_G2A_tmp,"G2A",discordant_C2T,discordant_G2A)
	multi_C2T,discordant_C2T,discordant_G2A = check_list(multi_C2T_tmp,"C2T",discordant_C2T,discordant_G2A)
	multi_G2A,discordant_C2T,discordant_G2A = check_list(multi_G2A_tmp,"G2A",discordant_C2T,discordant_G2A)
	
	if mode == "PE": #no more multimapping
		if len(unique_C2T) == 1:
			discordant_C2T.extend(unique_C2T)
			unique_C2T = []
		if len(unique_G2A) == 1:
			discordant_G2A.extend(unique_G2A)
			unique_G2A = []
		if len(multi_C2T) == 2:
			has_read1,has_read2 = False,False
			# has_read1,has_read2 = check_pairs(has_read1,has_read2,multi_C2T[0])
			#has_read1,has_read2 = check_pairs(has_read1,has_read2,multi_C2T[1])
			if validate_read_pairs(multi_C2T[0],multi_C2T[1]) == True:#has_read1 and has_read2:
				while multi_C2T:
					segment = multi_C2T.pop()
					segment.is_secondary = False
					segment.set_tag("NH",1)
					unique_C2T.append(segment)
			else:
				discordant_C2T.extend(multi_C2T)
				multi_C2T = []
		if len(multi_G2A) == 2:
			#has_read1,has_read2 = False,False
			#has_read1,has_read2 = check_pairs(has_read1,has_read2,multi_G2A[0])
			#has_read1,has_read2 = check_pairs(has_read1,has_read2,multi_G2A[1])
			if validate_read_pairs(multi_G2A[0],multi_G2A[1]) == True:#has_read1 and has_read2:
				while multi_G2A:
					segment = multi_G2A.pop()
					segment.is_secondary = False
					segment.set_tag("NH",1)
					unique_G2A.append(segment)
			else:
				discordant_G2A.extend(multi_G2A)
				multi_G2A = []
	elif mode == "SE":
		if len(multi_C2T) == 1:
			while multi_C2T:
				segment = multi_C2T.pop()
				segment.is_secondary = False
				segment.set_tag("NH",1)
				unique_C2T.append(segment)
		if len(multi_G2A) == 1:
			while multi_G2A:
				segment = multi_G2A.pop()
				segment.is_secondary = False
				segment.set_tag("NH",1)
				unique_G2A.append(segment)

	if unmapped_C2T and unmapped_G2A:
		for segment in unmapped_C2T:
			to_unmapped(segment,unal_read1,unal_read2)
		unmapped_num += 1
	elif multi_C2T or multi_G2A:
		to_bam_genome(multi_C2T,multi_output)
		to_bam_genome(multi_G2A,multi_output)
		to_bam_genome(unique_C2T,multi_output)
		to_bam_genome(unique_G2A,multi_output)
		if options.no_multi_fastq == False:
			if multi_C2T:
				to_unmapped_discordant(multi_C2T+discordant_C2T,unal_read1,unal_read2)
			elif multi_G2A:
				to_unmapped_discordant(multi_G2A+discordant_G2A,unal_read1,unal_read2)
		multimapper_num += 1
	elif unique_C2T and unique_G2A:
		''' experimental function '''
		if options.fully_unique == False:
			if mode == "PE":
				if unique_C2T[0].is_proper_pair == True and unique_G2A[0].is_proper_pair == False:
					for i in unique_C2T:
						is_poly = remove_high_poly(i)
						if is_poly == True:
							to_unmapped_discordant(unique_C2T,unal_read1,unal_read2)
							unmapped_num += 1
							return None
					to_bam_genome(unique_C2T,output)
					C2T_num += 1
				elif unique_C2T[0].is_proper_pair == False and unique_G2A[0].is_proper_pair == True:
					for i in unique_G2A:
						is_poly = remove_high_poly(i)
						if is_poly == True:
							to_unmapped_discordant(unique_G2A,unal_read1,unal_read2)
							unmapped_num += 1
							return None
					to_bam_genome(unique_G2A,output)
					G2A_num += 1
				else:
					C2T_AS = unique_C2T[0].get_tag("AS") + unique_C2T[1].get_tag("AS")
					G2A_AS = unique_G2A[0].get_tag("AS") + unique_G2A[1].get_tag("AS")
					if C2T_AS - G2A_AS > options.diff_AS:
						for i in unique_C2T:
							is_poly = remove_high_poly(i)
							if is_poly == True:
								to_unmapped_discordant(unique_C2T,unal_read1,unal_read2)
								unmapped_num += 1
								return None
						to_bam_genome(unique_C2T,output)
						C2T_num +=1 
					elif G2A_AS - C2T_AS > options.diff_AS:
						for i in unique_G2A:
							is_poly = remove_high_poly(i)
							if is_poly == True:
								to_unmapped_discordant(unique_G2A,unal_read1,unal_read2)
								unmapped_num += 1
								return None
						to_bam_genome(unique_G2A,output)
						G2A_num +=1 
					else:
						to_bam_genome(unique_C2T,multi_output)
						to_bam_genome(unique_G2A,multi_output)
						if options.no_multi_fastq == False:
							to_unmapped_discordant(unique_C2T,unal_read1,unal_read2)
						multimapper_num += 1
			elif mode == "SE":
				C2T_AS = unique_C2T[0].get_tag("AS")
				G2A_AS = unique_G2A[0].get_tag("AS")
				if C2T_AS - G2A_AS > options.diff_AS:
					for i in unique_C2T:
						is_poly = remove_high_poly(i)
						if is_poly == True:
							to_unmapped_discordant(unique_C2T,unal_read1,unal_read2)
							unmapped_num += 1
							return None
					to_bam_genome(unique_C2T,output)
					C2T_num +=1 
				elif C2T_AS - G2A_AS > options.diff_AS:
					for i in unique_G2A:
						is_poly = remove_high_poly(i)
						if is_poly == True:
							to_unmapped_discordant(unique_G2A,unal_read1,unal_read2)
							unmapped_num += 1
							return None
					to_bam_genome(unique_G2A,output)
					G2A_num +=1 
				else:
					to_bam_genome(unique_C2T,multi_output)
					to_bam_genome(unique_G2A,multi_output)
					if options.no_multi_fastq == False:
						to_unmapped_discordant(unique_C2T,unal_read1,unal_read2)
					multimapper_num += 1
		else:
			to_bam_genome(unique_C2T,multi_output)
			to_bam_genome(unique_G2A,multi_output)
			if options.no_multi_fastq == False:
				to_unmapped_discordant(unique_C2T,unal_read1,unal_read2)
			multimapper_num += 1
	elif unique_C2T and not unique_G2A:
		for i in unique_C2T:
			is_poly = remove_high_poly(i)
			if is_poly == True:
				to_unmapped_discordant(unique_C2T,unal_read1,unal_read2)
				unmapped_num += 1
				return None
		to_bam_genome(unique_C2T,output)
		C2T_num += 1
	elif unique_G2A and not unique_C2T:
		for i in unique_G2A:
			is_poly = remove_high_poly(i)
			if is_poly == True:
				to_unmapped_discordant(unique_G2A,unal_read1,unal_read2)
				unmapped_num += 1
				return None
		to_bam_genome(unique_G2A,output)
		G2A_num += 1
	elif any([unique_C2T,unique_G2A,multi_C2T,multi_G2A]) == False:
		if discordant_C2T:
			to_unmapped_discordant(discordant_C2T,unal_read1,unal_read2)
		elif discordant_G2A:
			to_unmapped_discordant(discordant_G2A,unal_read1,unal_read2)
		unmapped_num += 1

def itering_groups(C2T_iter,G2A_iter,C2T_next,G2A_next,\
				   output,unal_read1,unal_read2,multi_output,mode):
		unique_C2T = []
		unique_G2A = []
		unmapped_C2T = []
		unmapped_G2A = []
		multi_C2T = []
		multi_G2A = []
		#first read
		try:
			if not C2T_next and not G2A_next:
				C2T_read = C2T_iter.__next__()
				G2A_read = G2A_iter.__next__()
			else:
				C2T_read = C2T_next
				G2A_read = G2A_next
		except StopIteration:
			return None,None,None,None
			
		last_read = C2T_read.query_name
		C2T_group = C2T_read.get_tag("NH") if C2T_read.has_tag("NH") else None

		try:
			while C2T_read.query_name == last_read:
				if C2T_group == 1: #C2T unique
					unique_C2T.append(C2T_read)
				elif C2T_group is None: #Unmapped
					unmapped_C2T.append(C2T_read)
				else: #multimapper:
					multi_C2T.append(C2T_read)
				C2T_read = C2T_iter.__next__()
			C2T_next = C2T_read
		except StopIteration:
			pass
			
		G2A_group = G2A_read.get_tag("NH") if G2A_read.has_tag("NH") else None
		try:
			while G2A_read.query_name == last_read:
				if G2A_group == 1: #G2A unique
					unique_G2A.append(G2A_read)
				elif G2A_group is None: #Unmapped
					unmapped_G2A.append(G2A_read)
				else: #multimapper:
					multi_G2A.append(G2A_read)
				G2A_read = G2A_iter.__next__()
			G2A_next = G2A_read
			mapping_result(unique_C2T,unique_G2A,unmapped_C2T,unmapped_G2A,multi_C2T,multi_G2A,\
						   output,unal_read1,unal_read2,multi_output,mode)
			return C2T_iter,G2A_iter,C2T_next,G2A_next
		except StopIteration:
			mapping_result(unique_C2T,unique_G2A,unmapped_C2T,unmapped_G2A,multi_C2T,multi_G2A,\
						   output,unal_read1,unal_read2,multi_output,mode)
			C2T_next = None
			G2A_next = None
			return C2T_iter,G2A_iter,C2T_next,G2A_next

def sam_handler_PE():
	unal_1 = options.fwd.replace(".fastq",".unmapped.fastq")
	unal_2 = options.rev.replace(".fastq",".unmapped.fastq")
	if not options.continue_prefix:
		hisat2_prefix =  "hisat2_" + str(os.getpid())
		hisat2_C2T_SAM = hisat2_prefix + ".C2T.sam"
		hisat2_G2A_SAM = hisat2_prefix + ".G2A.sam"
	else:
		hisat2_C2T_SAM = options.continue_prefix+".C2T.sam"
		hisat2_G2A_SAM = options.continue_prefix+".G2A.sam"
		
	# with pysam.AlignmentFile("hisat2.C2T.sam","r") as C2T_SAM, \
		 # pysam.AlignmentFile("hisat2.G2A.sam","r") as G2A_SAM, \
		 # pysam.AlignmentFile(options.output+".bam","wb",template=C2T_SAM) as output, \
		 # open(unal_1,'w') as unal_read1, \
		 # open(unal_2,'w') as unal_read2, \
		 # pysam.AlignmentFile(options.output+".multimappers.bam","wb",template=C2T_SAM) as multi_output:
	with pysam.AlignmentFile(hisat2_C2T_SAM,"r") as C2T_SAM, \
		 pysam.AlignmentFile(hisat2_G2A_SAM,"r") as G2A_SAM, \
		 pysam.AlignmentFile(options.output+".bam","wb",template=C2T_SAM) as output, \
		 open(unal_1,'w') as unal_read1, \
		 open(unal_2,'w') as unal_read2, \
		 pysam.AlignmentFile(options.output+".multimappers.bam","wb",template=C2T_SAM) as multi_output:
		C2T_iter = C2T_SAM.fetch(until_eof=True)
		G2A_iter = G2A_SAM.fetch(until_eof=True)
		C2T_next = None
		G2A_next = None
		while (C2T_iter and G2A_iter):
			C2T_iter,G2A_iter,C2T_next,G2A_next = itering_groups(C2T_iter,G2A_iter,C2T_next,G2A_next,output,unal_read1,unal_read2,multi_output,"PE")

def itering_groups_SE(C2T_iter,C2T_next,\
				      output,unal_read1,multi_output,mode):
		unique_C2T = []
		unique_G2A = []
		unmapped_C2T = []
		unmapped_G2A = []
		multi_C2T = []
		multi_G2A = []
		#first read
		try:
			if not C2T_next and not G2A_next:
				C2T_read = C2T_iter.__next__()
			else:
				C2T_read = C2T_next
		except StopIteration:
			return None,None

		last_read = C2T_read.query_name
		C2T_group = C2T_read.get_tag("NH") if C2T_read.has_tag("NH") else None
		
		try:
			while C2T_read.query_name == last_read:
				if C2T_group == 1: #C2T unique
					unique_C2T.append(C2T_read)
				elif C2T_group is None: #Unmapped
					unmapped_C2T.append(C2T_read)
				else: #multimapper:
					multi_C2T.append(C2T_read)
				C2T_read = C2T_iter.__next__()
			C2T_next = C2T_read
		except StopIteration:
			mapping_result(unique_C2T,unique_G2A,unmapped_C2T,unmapped_G2A,multi_C2T,multi_G2A,\
						   output,unal_read1,unal_read2,multi_output,mode)
			C2T_next = None
			return C2T_iter,C2T_next
			
def sam_handler_SE():
	unal_1 = options.fwd.replace(".fastq",".unmapped.fastq")
	unal_2 = "single_end.ignore"
	if not options.continue_prefix:
		hisat2_prefix = "hisat2_" + str(os.getpid())
		hisat2_C2T_SAM = hisat2_prefix + ".C2T.sam"
		hisat2_G2A_SAM = hisat2_prefix + ".G2A.sam"
	else:
		hisat2_C2T_SAM = options.continue_prefix+".C2T.sam"
		hisat2_G2A_SAM = options.continue_prefix+".G2A.sam"
	with pysam.AlignmentFile(hisat2_C2T_SAM,"r") as C2T_SAM, \
		 pysam.AlignmentFile(hisat2_G2A_SAM,"r") as G2A_SAM, \
		 pysam.AlignmentFile(options.output+".bam","wb",template=C2T_SAM) as output, \
		 open(unal_1,'w') as unal_read1, \
		 open(unal_2,'w') as unal_read2, \
		 pysam.AlignmentFile(options.output+".multimappers.bam","wb",template=C2T_SAM) as multi_output:
		C2T_iter = C2T_SAM.fetch(until_eof=True)
		G2A_iter = G2A_SAM.fetch(until_eof=True)
		C2T_next = None
		G2A_next = None
		while (C2T_iter and G2A_iter):
			C2T_iter,G2A_iter,C2T_next,G2A_next = itering_groups(C2T_iter,G2A_iter,C2T_next,G2A_next,output,unal_read1,unal_read2,multi_output,"SE")
	os.remove(unal_2)
	
def signal_handler(sig,frame):
	try:
		pool.terminate()
	except NameError:
		pass
	sys.exit()

if __name__ == "__main__":
	parser = argparse.ArgumentParser(prog="BS_hisat2",fromfile_prefix_chars='@',description=__doc__,formatter_class=argparse.RawTextHelpFormatter)
	#Required
	group_required = parser.add_argument_group("Required")
	group_required.add_argument("--fwd","-F",dest="fwd",required=True,help="Forward read (C less read)")
	group_required.add_argument("--rev","-R",dest="rev",required=False,help="Reverse read (G less read)")
	group_required.add_argument("--index","-I",dest="index",required=True,help="Path to index folder")
	group_required.add_argument("--output","-o",dest="output",required=True,help="Output prefix")
	#Hisat2
	group_hisat2 = parser.add_argument_group("Hisat2")
	group_hisat2.add_argument("--index-prefix",dest="index_prefix",required=False,default="HISAT2",help="Index will looks like /path/C2T/[prefix]_C2T and /path/G2A/[prefix]_G2A")
	group_hisat2.add_argument("--hisat2-path",dest="hisat2_path",help="Path to hisat2 bin")
	group_hisat2.add_argument("--hisat2-param",dest="hisat2_param",help="Index building parameters")
	#Filter
	group_filter = parser.add_argument_group("Filter")
	group_filter.add_argument("--min-AS",dest="min_AS",type=int,default=-10,help="minimal AS")
	group_filter.add_argument("--diff-AS",dest="diff_AS",type=int,default=9,help="AS different between best and secondary mapping")
	group_filter.add_argument("--fully-unique",dest="fully_unique",default=False,action="store_true",help="Reads should be aligned to C2T or G2A by only ONE time")
	#Steps
	group_steps = parser.add_argument_group("Steps")
	group_steps.add_argument("--continue-prefix",dest="continue_prefix",help="If given, process the sam files: [preifx].C2T.sam and [prefix].G2A.sam. If not given, prefix=hisat2_[pid]")
	group_steps.add_argument("--skip-convert",dest="skip_convert",default=False,action="store_true",help="Skip converting fastq")
	group_steps.add_argument("--skip-mapping",dest="skip_mapping",default=False,action="store_true",help="Skip mapping")
	group_steps.add_argument("--sorted",dest="sorted_bam",default=False,action="store_true",help="Do not sort bam")
	group_steps.add_argument("--no-index",dest="no_sorted_bam_index",default=False,action="store_true",help="Do not index bam")
	#Files
	group_files = parser.add_argument_group("Files")
	group_files.add_argument("--del-convert",dest="del_convert",default=False,action="store_true",help="Do ot delete the converted fastq files")
	group_files.add_argument("--del-sam",dest="del_sam",default=False,action="store_true",help="Delete the sam files")
	group_files.add_argument("--del-bam",dest="del_bam",default=False,action="store_true",help="Delete the unsorted bam")
	group_files.add_argument("--no-multi-fastq",dest="no_multi_fastq",default=False,action="store_true",help="No multimapping in unmapped fastq")
	#Other
	group_other = parser.add_argument_group("Other")
	group_other.add_argument("--version",action="version",version="%(prog)s 1.0")
	options = parser.parse_args()
	
	signal.signal(signal.SIGINT,signal_handler)
	multimapper_num = 0
	unmapped_num = 0
	unique_num = 0
	C2T_num = 0
	G2A_num = 0
	
	if options.continue_prefix:
		hisat2_prefix = options.continue_prefix
	else:
		hisat2_prefix = "hisat2_" + str(os.getpid())
	fwd_c2t = None
	rev_g2a = None
	
	if options.fwd:
		if options.fwd.endswith(".fastq") == False:
			raise ValueError("Forward file do not end with .fastq")
		fwd_c2t = options.fwd.replace(".fastq",".c2t.fastq")
	if options.rev:
		if options.rev.endswith(".fastq") == False:
			raise ValueError("Reverse file do not end with .fastq")
		rev_g2a = options.rev.replace(".fastq",".g2a.fastq")
	
	#Step1 convert transcripts
	if options.skip_convert == False:
		if options.fwd:
			fwd_read_dict,total_reads_fwd = fastq_converter_1(options.fwd,fwd_c2t,"C2T")
		if options.rev:
			rev_read_dict,total_reads_rev = fastq_converter_1(options.rev,rev_g2a,"G2A")
			if total_reads_fwd != total_reads_rev:
				raise Warning("Read number not equal!")
		else:
			total_reads = total_reads_fwd
	
	#Step2 mapping
	if options.skip_mapping == False:
		sys.stderr.write("[%s]Mapping with hisat2, TEMP prefix: %s\n" % (strftime("%Y-%m-%d %H:%M:%S", time.localtime()),hisat2_prefix))
		hisat2 = options.hisat2_path
		hisat2_parameters = read_parameters(options.hisat2_param)
		hisat2_parameters = mapping_string(hisat2_parameters,fwd_c2t,rev_g2a)
		fwd_index = options.index + "/C2T/" + options.index_prefix + "_C2T"
		rev_index = options.index + "/G2A/" + options.index_prefix + "_G2A"
		
		if options.fwd:
			pool = multiprocessing.Pool(2)
			try:
				pool.apply_async(mapping, args=(fwd_index,"C2T",hisat2,hisat2_parameters,hisat2_prefix,))
				pool.apply_async(mapping, args=(rev_index,"G2A",hisat2,hisat2_parameters,hisat2_prefix,))
				pool.close()
				pool.join()
			finally:
				pool.terminate()
	else:
		if not options.continue_prefix:
			raise Warning("Please provide a SAM prefix.")
	
	if options.skip_convert:
		#Step get C-contained reads
		sys.stderr.write("[%s]Finding non-converted C containing reads...\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
		fwd_read_dict,rev_read_dict,total_reads = read_C_position(options.fwd,options.rev)
	
	if options.del_convert == True:
		os.remove(fwd_c2t)
		if rev_g2a:
			os.remove(rev_g2a)
	
	#Step Get unique mapping and report unmapped and multimappers
	sys.stderr.write("[%s]Handling SAM outputs...\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
	if options.fwd and options.rev:
		sam_handler_PE()
	elif options.fwd and not options.rev:
		sam_handler_SE()

	if options.del_sam == True:
		if options.fwd:
			os.remove(hisat2_prefix + ".C2T.sam")
		if options.rev:
			os.remove(hisat2_prefix + ".G2A.sam")
		# os.remove("hisat2.C2T.sam")
		# os.remove("hisat2.G2A.sam")

	unique_num = C2T_num + G2A_num
	total_reads = unique_num + multimapper_num + unmapped_num
	sys.stderr.write("[%s]Completed successfully:\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
	sys.stderr.write(" Total reads: %d\n" % total_reads)
	sys.stderr.write(" Unique mapping: %d (%.3f%%)\n" % (unique_num,100*unique_num/(total_reads+0.0)))
	sys.stderr.write("   C2T: %d (%.2f%%)\n" % (C2T_num,100*C2T_num/(total_reads+0.0)))
	sys.stderr.write("   G2A: %d (%.2f%%)\n" % (G2A_num,100*G2A_num/(total_reads+0.0)))
	sys.stderr.write(" Multiple mapping: %d (%.3f%%)\n" % (multimapper_num,100*multimapper_num/(total_reads+0.0)))
	sys.stderr.write(" Unmapped: %d (%.3f%%)\n" % (unmapped_num,100*unmapped_num/(total_reads+0.0)))
	
	if options.sorted_bam == True:
		sys.stderr.write("[%s]Sorting bam...\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
		pysam.sort("-o", options.output+".sorted.bam", options.output+".bam")
		if options.no_sorted_bam_index == False:
			sys.stderr.write("[%s]Indexing bam...\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
			pysam.index(options.output+".sorted.bam")

	if options.del_bam == True:
		os.remove(options.output+".bam")