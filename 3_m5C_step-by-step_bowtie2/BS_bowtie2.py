#!bin/usr/env python

#Jianheng Liu @ Zhanglab, SYSU
#Feb, 2018
#Email: liujh26@mail2.sysu.edu.cn
#Usage: Map reads to C2T references with bowtie2
#Input: [.fa] [.genelist]

__doc__ = """
Map reads to C2T references with bowte2. If the conversion
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
import copy

def read_transcripts_info(fn):
	output = {}
	gene_biotype = {}
	with open(fn,'r') as input:
		for line in input.readlines():
			line = line.strip().split("\t")
			enst = line[2]
			ensg = line[3]
			biotype = line[4]
			strand = line[8]
			length = line[9]
			try:
				length = int(length)
			except ValueError:
				length = 0
			except TypeError:
				length = 0
			output[enst] = {'ENSG':ensg,
							'Biotype':biotype,
							'length':length,
							'strand':strand,
			                }
			gene_biotype[ensg] = biotype
	return output,gene_biotype
	
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
	with open(fin,'r') as input, file(fout,'w') as output:
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
	# if "--reorder" not in parameters:
		# parameters["--reorder"] = None
	if "-P" not in parameters:
		parameters["-p"] = '8'
	if "-k" not in parameters and "-a" not in parameters:
		parameters["-k"] = '10'
		# parameters["-a"] = None
	if "--end-to-end" or "--local" not in parameters:
		parameters["--end-to-end"] = None
		# parameters["--local"] = None
	if "--gbar" not in parameters:
		parameters["--gbar"] = '2'
	if "--mp" not in parameters:
		parameters["--mp"] = '5'
	if "-R" not in parameters:
		parameters["-R"] = '2'
	if "-D" not in parameters:
		parameters["-D"] = '5'
	if fwd and rev:
		parameters["-1"] = fwd
		parameters["-2"] = rev
		if "--fr" not in parameters:
			parameters["--fr"] = None
	if fwd and not rev:
		parameters["-U"] = fwd
		parameters["--norc"] = None
		

	para = []
	for key,values in parameters.items():
		if values is not None:
			para.extend([key,values])
		else:
			para.extend([key])
	
	para = "|".join(para)
	return para
	
def mapping(index,method,bowtie2_path,parameters):
	if parameters:
		parameters = parameters.split("|")
	else:
		parameters = []
	
	sam = "bowtie2_" + str(os.getpid()) + ".sam"
	
	cmd = [bowtie2_path+"/bowtie2", \
		    "-x",index,\
			# "--no-mixed",\
			"-S",sam] + \
		   parameters

	# print cmd
	bowtie2_out = subprocess.Popen(cmd,\
								   stdin=subprocess.PIPE,\
								   stdout=subprocess.PIPE,\
								   stderr=subprocess.PIPE,\
								   shell=False)
	stdout,stderr = bowtie2_out.communicate()
	# print stdout
	sys.stderr.write("[%s]%s report:\n" % (strftime("%Y-%m-%d %H:%M:%S", time.localtime()),method))
	sys.stderr.write(stderr)
	sys.stderr.write("\n")

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
	else: #single end
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

def to_unmapped_discordant(discordant,unal_read1,unal_read2,mode):
	read1 = None
	read2 = None
	discordant_iter = iter(discordant)
	
	while read1 is None or read2 is None:
		try:
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
		except StopIteration:
			print discordant
			print discordant[0]
			print read1,read2,mode
			sys.exit()
			break

def to_bam_transcriptome(segments,output,unique=False):
	for segment in segments:
		if unique == True:
			segment.is_secondary = False
		# segment.set_tag("YG",genome)
		if segment is not None:
			if segment.is_read1:
				seq = fwd_read_dict.get(segment.query_name)
			elif segment.is_read2:
				seq = rev_read_dict.get(segment.query_name)
			else:
				seq = fwd_read_dict.get(segment.query_name)
			if seq is not None:
				if segment.is_reverse:
					seq = reverse_complement(str(seq))
				segment.set_tag("CT",seq.count("C"))
				qual = segment.query_qualities
				mpq = segment.mapping_quality
				segment_output = segment
				segment_output.query_sequence = seq
				segment_output.mapping_quality = mpq
				segment_output.query_qualities = qual
				output.write(segment_output)
			else:
				segment.set_tag("CT",0)
				output.write(segment)
				
def check_mapping(segment_check):
	if segment_check.is_read1 and segment_check.is_reverse == False and segment_check.get_tag("AS") >= options.AS_min_SE:
		return True,segment_check
	elif segment_check.is_read2 and segment_check.is_reverse == True and segment_check.get_tag("AS") >= options.AS_min_SE:
		return True,segment_check
	elif segment_check.is_paired == False and segment_check.is_reverse == False and segment_check.get_tag("AS") >= options.AS_min_SE:
		return True,segment_check
	else:
		return False,segment_check
		
def set_tag(segment_tag):
	try:
		strand = transcripts_info[segment_tag.reference_name]["strand"]
		segment_tag.set_tag("TS",strand)
		if strand == "+":
			segment_tag.set_tag("YG","C2T")
		elif strand == "-":
			segment_tag.set_tag("YG","G2A")
		# else:
			# raise ValueError("Strand not correct: %s" % segment_tag.reference_name)
	except KeyError:
		segment_tag.set_tag("YG","C2T")

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
	
def classifer(multi_C2T,mode):
	read1 = {}
	read2 = {}
	comb = {}
	comb_best = {}
	multi_C2T_tmp = copy.copy(multi_C2T)
	while multi_C2T_tmp:
		segment = multi_C2T_tmp.pop()
		if segment.is_read1:
			try:
				read1[segment.reference_name].append(copy.copy(segment))
			except KeyError:
				read1[segment.reference_name] = []
				read1[segment.reference_name].append(copy.copy(segment))
		elif segment.is_read2:
			try:
				read2[segment.reference_name].append(copy.copy(segment))
			except KeyError:
				read2[segment.reference_name] = []
				read2[segment.reference_name].append(copy.copy(segment))
		else:
			try:
				read1[segment.reference_name].append(copy.copy(segment))
			except KeyError:
				read1[segment.reference_name] = []
				read1[segment.reference_name].append(copy.copy(segment))
	for read1_key in read1.keys():
		info = transcripts_info.get(read1_key)
		if mode == "PE":
			read2_values = read2.get(read1_key)
			if read2_values is not None:
				read1_values = read1.get(read1_key)
				if info:
					gene = info["ENSG"]
					biotype = info["Biotype"]
					length = info["length"]
				else:
					gene = "Unkonwn"
					biotype = "Unkonwn"
					length = 0
				for i in read1_values:
					for j in read2_values:
						if j.reference_start >= i.reference_start and j.reference_start - i.reference_end <= 500:
							AS = i.get_tag("AS") + j.get_tag("AS")
							if AS >= options.AS_min_PE: #end-to-end -20
								try:
									comb[read1_key].append((copy.copy(i),copy.copy(j),AS,gene,biotype,length,read1_key))
								except KeyError:
									comb[read1_key] =[]
									comb[read1_key].append((copy.copy(i),copy.copy(j),AS,gene,biotype,length,read1_key))
								if read1_key not in comb_best:
									comb_best[read1_key] = AS
								elif AS > comb_best[read1_key]:
									comb_best[read1_key] = AS
						elif i.reference_start < j.reference_end-1:
							AS = i.get_tag("AS") + j.get_tag("AS")
							if AS >= options.AS_min_PE: #end-to-end
								overlap = (j.reference_end - i.reference_start)/(i.reference_end-j.reference_start+0.0)
								if overlap >= 0.1:
									try:
										comb[read1_key].append((copy.copy(i),copy.copy(j),AS,gene,biotype,length,read1_key))
									except KeyError:
										comb[read1_key] = []
										comb[read1_key].append((copy.copy(i),copy.copy(j),AS,gene,biotype,length,read1_key))
									if read1_key not in comb_best:
										comb_best[read1_key] = AS
									elif AS > comb_best[read1_key]:
										comb_best[read1_key] = AS
		elif mode == "SE":
			read1_values = read1[read1_key]
			for i in read1_values:
				if info:
					gene = info["ENSG"]
					biotype = info["Biotype"]
					length = info["length"]
				else:
					gene = "Unkonwn"
					biotype = "Unkonwn"
					length = 0
				AS = i.get_tag("AS")
				if AS >= options.AS_min_SE:
					try:
						comb[read1_key].append((copy.copy(i),None,AS,gene,biotype,length,read1_key))
					except KeyError:
						comb[read1_key] = []
						comb[read1_key].append((copy.copy(i),None,AS,gene,biotype,length,read1_key))
					if read1_key not in comb_best:
						comb_best[read1_key] = AS
					elif AS > comb_best[read1_key]:
						comb_best[read1_key] = AS

	comb_best_AS = -99999
	comb_second_AS = -99999
	for key,value in comb_best.items():
		if value > comb_best_AS:
			comb_best_AS = value
		# elif value > comb_second_AS:
			# comb_second_AS = value
		# if comb_best_AS - comb_second_AS < options.AS_diff:
			# cannot_distinguish += 1
			# return None,None
	
	gene_isoforms = {}
	for isoform,candidate in comb.items():
		if comb_best[isoform] == comb_best_AS:
			for i in candidate:
				if i[2] >= comb_best_AS - options.AS_diff:
					try:
						if transcripts_info[isoform]['ENSG'] not in gene_isoforms:
							gene_isoforms[transcripts_info[isoform]['ENSG']] = {}
						if isoform not in gene_isoforms[transcripts_info[isoform]['ENSG']]:
							gene_isoforms[transcripts_info[isoform]['ENSG']][isoform] = []
						gene_isoforms[transcripts_info[isoform]['ENSG']][isoform].append(i)
					except KeyError:
						if isoform not in gene_isoforms:
							gene_isoforms[isoform] = {}
						if isoform not in gene_isoforms[isoform]:
							gene_isoforms[isoform][isoform] = []
						gene_isoforms[isoform][isoform].append(i)

	gene_isoforms_len = len(gene_isoforms.keys())

	if gene_isoforms_len == 0:
		# AS_socre_too_low += 1
		return None,None,"AS_socre_too_low"
	elif gene_isoforms_len == 1:
		gene,isoforms = next(gene_isoforms.iteritems())
		if len(isoforms.keys()) == 1:
			isoform_name,isoform_values = next(isoforms.iteritems())
			if len(isoform_values) == 1:
				# unique_isoform_num += 1
				return isoform_values[0][0],isoform_values[0][1],"unique_isoform"
			else:
				# multi_location += 1
				return None,None,"multi_location"
		else:
			if options.no_longest == False:
				longest = -999
				longest_isoform = None
				for isoform in isoforms.keys():
					try:
						isoform_len = transcripts_info[isoform]['length']
					except KeyError:
						isoform_len = 0
					if isoform_len > longest:
						longest = transcripts_info[isoform]['length']
						longest_isoform = isoform
					elif isoform_len == longest:
						# diff_gene += 1
						return None,None,"diff_gene"
				isoform_values = isoforms[longest_isoform]
				if len(isoform_values) == 1:
					# longest_isoform_num += 1
					return isoform_values[0][0],isoform_values[0][1],"longest_isoform"
				else:
					# diff_isoform += 1
					return None,None,"diff_isoform"
			else:
				return None,None,"diff_isoform"
	else:
		if options.no_coding == False:
			coding = []
			for gene in gene_isoforms.keys():
				try:
					if gene_biotype[gene] == "protein_coding":
						coding.append(gene)
				except KeyError:
					pass
			if len(coding) == 1:
				gene = coding[0]
				isoforms = gene_isoforms[gene]
				# gene,isoforms = next(gene_isoforms.iteritems())
				if len(isoforms.keys()) == 1:
					isoform_name,isoform_values = next(isoforms.iteritems())
					if len(isoform_values) == 1:
						# coding_isoform_num += 1
						return isoform_values[0][0],isoform_values[0][1],"coding"
					else:
						# diff_isoform += 1
						return None,None,"diff_isoform"
				else:
					if options.no_longest == False:
						longest = -999
						longest_isoform = None
						for isoform in isoforms.keys():
							try:
								isoform_len = transcripts_info[isoform]['length']
							except KeyError:
								isoform_len = 0
							if isoform_len > longest:
								longest = transcripts_info[isoform]['length']
								longest_isoform = isoform
							elif isoform_len == longest:
								# diff_gene += 1
								return None,None,"diff_gene"
						isoform_values = isoforms[longest_isoform]
						if len(isoform_values) == 1:
							# coding_isoform_longest_num += 1
							return isoform_values[0][0],isoform_values[0][1],"coding_isoform_longest"
						else:
							# diff_isoform += 1
							return None,None,"diff_isoform"
					else:
						return None,None,"diff_isoform"
			else:
				# diff_gene += 1
				return None,None,"diff_gene"
		else:
			return None,None,"diff_gene"

def status_handler(mode,status):
	if status == "diff_gene":
		stats["diff_gene"] += 1
	elif status == "diff_isoform":
		stats["diff_isoform"] += 1
	elif status == "multi_location":
		stats["multi_location"] += 1
	elif status == "AS_socre_too_low":
		stats["AS_socre_too_low"] += 1
		
	elif status == "unique_isoform":
		if mode == "PE":
			stats["unique_isoform_PE"] += 1
		elif mode == "SE":
			stats["unique_isoform_SE"] += 1
	elif status == "coding":
		if mode == "PE":
			stats["coding_PE"] += 1
		elif mode == "SE":
			stats["coding_SE"] += 1
	elif status == "longest_isoform":
		if mode == "PE":
			stats["longest_isoform_PE"] += 1
		elif mode == "SE":
			stats["longest_isoform_SE"] += 1
	elif status == "coding_isoform_longest":
		if mode == "PE":
			stats["coding_isoform_longest_PE"] += 1
		elif mode == "SE":
			stats["coding_isoform_longest_SE"] += 1
	else:
		raise KeyError(status)

def mapping_result(unique_C2T,unmapped_C2T,multi_C2T,discordant_C2T,\
                   output,unal_read1,unal_read2,multi_output,mode,three_prime=False):
	global three_prime_num,PE_num,SE_num,multimapper_num,C2T_num,unmapped_num

	if mode == "PE" and len(unique_C2T) == 2:
		stats["fully_unique_PE"] += 1
	elif mode == "SE" and len(unique_C2T) == 1:
		stats["fully_unique_SE"] += 1
	elif multi_C2T:
		unique_1,unique_2 = None,None
		# try:
		unique_1,unique_2,status = classifer(multi_C2T,mode)
		status_handler(mode,status)
		# except TypeError:
			# print mode
			# for i in multi_C2T:
				# print i
			# sys.exit()
		if mode == "PE" and unique_1 and unique_2:
			unique_C2T.append(unique_1)
			unique_C2T.append(unique_2)
			multi_C2T = []
		elif mode == "SE" and unique_1:
			unique_C2T.append(unique_1)
			multi_C2T = []
	if multi_C2T:
		to_bam_transcriptome(multi_C2T,multi_output)
		if options.no_multi_fastq == False:
			to_unmapped_discordant(multi_C2T+unique_C2T+unmapped_C2T+discordant_C2T,unal_read1,unal_read2,mode)
		multimapper_num += 1
	elif unique_C2T:
		to_bam_transcriptome(unique_C2T,output,unique=True)
		if mode == "PE":
			PE_num += 1
		elif mode == "SE":
			if three_prime == True:
				three_prime_num += 1
				stats["three_prime"] += 1
			SE_num += 1
		C2T_num += 1
	elif mode == "PE" and len(unmapped_C2T) == 2:
		to_unmapped(unmapped_C2T[0],unal_read1,unal_read2)
		to_unmapped(unmapped_C2T[1],unal_read1,unal_read2)
		unmapped_num += 1
	elif mode == "SE" and len(unmapped_C2T) == 1:
		to_unmapped(unmapped_C2T[0],unal_read1,unal_read2)
		unmapped_num += 1
	else:
		to_unmapped_discordant(discordant_C2T+unique_C2T+unmapped_C2T+multi_C2T,unal_read1,unal_read2,mode)
		unmapped_num += 1

def check_list(segment,discordant_C2T):
	if segment.is_unmapped:
		return 0,discordant_C2T
	segment = set_tag(segment)
	status, segment = check_mapping(segment)
	if status == True:
		return 1,discordant_C2T
	elif status == False:
		discordant_C2T.append(segment)
		return 2,discordant_C2T
	else:
		raise Warning("Status error")

def set_discordant(read1_tmp,read2_tmp,unique_C2T,unmapped_C2T,multi_C2T,discordant_C2T):
	# three_prime_len = None
	discordant_C2T = read1_tmp + read2_tmp + unique_C2T + unmapped_C2T + multi_C2T + discordant_C2T
	# unique_C2T = []
	# multi_C2T = []
	# unmapped_C2T = []
	return None,[],[],[],discordant_C2T
	#three_prime_len,unique_C2T,multi_C2T,unmapped_C2T,discordant_C2T = set_discordant(read1_tmp,read2_tmp,unique_C2T,unmapped_C2T,multi_C2T,discordant_C2T)
	
def check_read1_and_read2(read1_tmp,read2_tmp,discordant_C2T,\
                          read1_unmapped, read2_unmapped,read1_len,read2_len,\
                          output,unal_read1,unal_read2,multi_output,mode):
	unique_C2T = []
	unmapped_C2T = []
	multi_C2T = []
	three_prime_len = None
	
	if read1_unmapped:
		unmapped_C2T.append(read1_unmapped)
	if read2_unmapped:
		unmapped_C2T.append(read2_unmapped)
	
	if mode == "PE":
		if not read1_unmapped and not read2_unmapped:
			if read1_len == read2_len == 1:
				unique_C2T = read1_tmp + read2_tmp
				multi_C2T = []
			elif read1_len > 0 and read2_len > 0:
				multi_C2T = read1_tmp + read2_tmp
				discordant_C2T = []
			elif read2_len == 0 and read1_len > 0:
				# if read1_len != 0 and read1_tmp[0].query_name == "E00552:46:HL7MJALXX:8:1104:23825:22440":
					# print "A"
				try:
					polyA = read2_unmapped.query_sequence.count("A") / (len(read2_unmapped.query_sequence) + 0.0)
				except AttributeError:
					tmp = unique_C2T + discordant_C2T + multi_C2T
					for segment in tmp:
						if segment.is_read2:
							read2_unmapped = segment
							seq = read2_unmapped.query_sequence
							if read2_unmapped.is_reverse:
								seq = reverse_complement(seq)
							break
					polyA = seq.count("A")/(len(seq)+0.0)
				if polyA >= 0.8:
					three_prime_reads = []
					three_prime_len = 0
					for read1 in read1_tmp:
						try:
							ref_len = transcripts_info[read1.reference_name]["length"]
							if ref_len - read1.reference_end < options.length:
								three_prime_reads.append(read1)
								three_prime_len += 1
							else:
								discordant_C2T.append(read1)
						except KeyError:
							discordant_C2T.append(read1)
					if options.no_three_prime == False:
						if three_prime_len == 1:
							unique_C2T = three_prime_reads
						elif three_prime_len > 1:
							multi_C2T = three_prime_reads
						else:
							three_prime_len,unique_C2T,multi_C2T,unmapped_C2T,discordant_C2T = set_discordant(read1_tmp,read2_tmp,unique_C2T,unmapped_C2T,multi_C2T,discordant_C2T)
					else:
						three_prime_len,unique_C2T,multi_C2T,unmapped_C2T,discordant_C2T = set_discordant(read1_tmp,read2_tmp,unique_C2T,unmapped_C2T,multi_C2T,discordant_C2T)
				else:
					three_prime_len,unique_C2T,multi_C2T,unmapped_C2T,discordant_C2T = set_discordant(read1_tmp,read2_tmp,unique_C2T,unmapped_C2T,multi_C2T,discordant_C2T)
			else:
				three_prime_len,unique_C2T,multi_C2T,unmapped_C2T,discordant_C2T = set_discordant(read1_tmp,read2_tmp,unique_C2T,unmapped_C2T,multi_C2T,discordant_C2T)
					
		elif read1_tmp and not read2_tmp:
			try:
				polyA = read2_unmapped.query_sequence.count("A") / (len(read2_unmapped.query_sequence) + 0.0)
			except AttributeError:
				for segment in discordant_C2T:
					if segment.is_read2:
						read2_unmapped = segment
						seq = read2_unmapped.query_sequence
						if read2_unmapped.is_reverse:
							seq = reverse_complement(seq)
						break
				polyA = seq.count("A") / (len(read2_unmapped.query_sequence) + 0.0)
			if polyA >= 0.8:
				three_prime_reads = []
				three_prime_len = 0
				for read1 in read1_tmp:
					try:
						ref_len = transcripts_info[read1.reference_name]["length"]
						if ref_len - read1.reference_end < options.length:
							three_prime_reads.append(read1)
							three_prime_len += 1
						else:
							discordant_C2T.append(read1)
					except KeyError:
						discordant_C2T.append(read1)

				if options.no_three_prime == False:
					if three_prime_len == 1:
						unique_C2T = three_prime_reads
					elif three_prime_len > 1:
						multi_C2T = three_prime_reads
					else:
						three_prime_len,unique_C2T,multi_C2T,unmapped_C2T,discordant_C2T = set_discordant(read1_tmp,read2_tmp,unique_C2T,unmapped_C2T,multi_C2T,discordant_C2T)
				else:
					three_prime_len,unique_C2T,multi_C2T,unmapped_C2T,discordant_C2T = set_discordant(read1_tmp,read2_tmp,unique_C2T,unmapped_C2T,multi_C2T,discordant_C2T)
			else:
				three_prime_len,unique_C2T,multi_C2T,unmapped_C2T,discordant_C2T = set_discordant(read1_tmp,read2_tmp,unique_C2T,unmapped_C2T,multi_C2T,discordant_C2T)
		elif not read2_unmapped and read1_unmapped:
			three_prime_len,unique_C2T,multi_C2T,unmapped_C2T,discordant_C2T = set_discordant(read1_tmp,read2_tmp,unique_C2T,unmapped_C2T,multi_C2T,discordant_C2T)
		else:
			three_prime_len,unique_C2T,multi_C2T,unmapped_C2T,discordant_C2T = set_discordant(read1_tmp,read2_tmp,unique_C2T,unmapped_C2T,multi_C2T,discordant_C2T)

		if three_prime_len is None:
			mapping_result(unique_C2T,unmapped_C2T,multi_C2T,discordant_C2T,\
						   output,unal_read1,unal_read2,multi_output,"PE")
		else:
			mapping_result(unique_C2T,unmapped_C2T,multi_C2T,discordant_C2T,\
						   output,unal_read1,unal_read2,multi_output,"SE",three_prime=True)
	elif mode == "SE":
		if read1_len == 1:
			unique_C2T = read1_tmp
		elif read1_len > 1:
			multi_C2T = read1_tmp
		# else:
			# print unmapped_C2T
			# raise Warning("SE error")
		mapping_result(unique_C2T,unmapped_C2T,multi_C2T,discordant_C2T,\
				   output,unal_read1,unal_read2,multi_output,"SE")
	
def itering_groups(C2T_iter,C2T_next,\
				   output,unal_read1,unal_read2,multi_output,mode):

	discordant_C2T = []
	#first read
	try:
		if not C2T_next:
			C2T_read = C2T_iter.__next__()
		else:
			C2T_read = C2T_next
	except StopIteration:
		return None,None
	
	last_read = C2T_read.query_name
	read1_unmapped = None
	read2_unmapped = None
	read1_tmp = []
	read2_tmp = []
	read1_len = 0
	read2_len = 0
	
	try:
		while C2T_read.query_name == last_read:
			passed,discordant_C2T = check_list(C2T_read,discordant_C2T)
			if passed == 0:
				if C2T_read.is_read1:
					read1_unmapped = C2T_read
				elif C2T_read.is_read2:
					read2_unmapped = C2T_read
				else:
					read1_unmapped = C2T_read
			elif passed == 1:
				if C2T_read.is_read1:
					read1_tmp.append(C2T_read)
					read1_len += 1
				elif C2T_read.is_read2:
					read2_tmp.append(C2T_read)
					read2_len += 1
				else:
					read1_tmp.append(C2T_read)
					read1_len += 1
					
			# if C2T_read.query_name == "E00552:46:HL7MJALXX:8:1104:23825:22440":
				# print read1_tmp,read2_tmp,discordant_C2T,read1_unmapped, read2_unmapped
					
			C2T_read = C2T_iter.__next__()
		
		check_read1_and_read2(read1_tmp,read2_tmp,discordant_C2T,\
		                      read1_unmapped, read2_unmapped,read1_len,read2_len,\
		                      output,unal_read1,unal_read2,multi_output,mode)
		C2T_next = C2T_read
		return C2T_iter,C2T_next
	except StopIteration:
		# if read1_tmp or read2_tmp:
			# discordant_C2T = []
		check_read1_and_read2(read1_tmp,read2_tmp,discordant_C2T,\
		                      read1_unmapped, read2_unmapped,read1_len,read2_len,\
		                      output,unal_read1,unal_read2,multi_output,mode)
		C2T_next = None
		return C2T_iter,C2T_next
		
def sam_handler(bowtie2_name,mode="PE"):
	unal_1 = options.fwd.replace(".fastq",".unmapped.fastq")
	if mode == "PE":
		unal_2 = options.rev.replace(".fastq",".unmapped.fastq")
	else:
		unal_2 = "single_end.ignore"
		
	with pysam.AlignmentFile(bowtie2_name,"r") as C2T_SAM, \
		 pysam.AlignmentFile(options.output+".bam","wb",template=C2T_SAM) as output, \
		 file(unal_1,'w') as unal_read1, \
		 file(unal_2,'w') as unal_read2, \
		 pysam.AlignmentFile(options.output+".multimappers.bam","wb",template=C2T_SAM) as multi_output:
		C2T_iter = C2T_SAM.fetch(until_eof=True)
		C2T_next = None
		while (C2T_iter):
			C2T_iter,C2T_next = itering_groups(C2T_iter,C2T_next,output,unal_read1,unal_read2,multi_output,mode)
	if mode == "SE":
		os.remove(unal_2)
		
def signal_handler(sig,frame):
	try:
		pool.terminate()
	except NameError:
		pass
	sys.exit()
	
if __name__ == "__main__":
	parser = argparse.ArgumentParser(prog="m5C_mapper",fromfile_prefix_chars='@',description=__doc__,formatter_class=argparse.RawTextHelpFormatter)
	#Require
	group_required = parser.add_argument_group("Required")
	group_required.add_argument("--fwd","-F",dest="fwd",required=True,help="Forward read (C less read)")
	group_required.add_argument("--rev","-R",dest="rev",required=False,help="Reverse read (G less read), if not given, use single end mode")
	group_required.add_argument("--index","-I",dest="index",required=False,help="Path to index")
	group_required.add_argument("--gene","-g",dest="genelist",required=True,help=".genelist")
	group_required.add_argument("--output","-o",dest="output",required=True,help="Output prefix")
	
	group_bowtie2 = parser.add_argument_group("Bowtie2")
	group_bowtie2.add_argument("--bowtie2-path",dest="bowtie2_path",help="Path to bowtie2 bin")
	group_bowtie2.add_argument("--bowtie2-param",dest="bowtie2_param",help="Index building parameters")
	
	group_filter = parser.add_argument_group("Filter")
	group_filter.add_argument("--library-length",dest="length",type=int,default=500,help="Possible library length, default=500")
	group_filter.add_argument("--AS-SE",dest="AS_min_SE",type=int,default=-15,help="Minimal AS for single end mapping, default=-8, quote it to pass if minus")
	group_filter.add_argument("--AS-PE",dest="AS_min_PE",type=int,default=-20,help="Minimal AS sum for paired end mapping, default=-15, quote it to pass if minus")
	group_filter.add_argument("--AS-diff",dest="AS_diff",type=int,default=5,help="If AS[primary mapping] - AS[secondary mapping] < this value, classify them as multimaper , default=5")
	group_filter.add_argument("--no-longest",dest="no_longest",default=False,action="store_true",help="Do not select the longest isoform in multimapper filtering")
	group_filter.add_argument("--no-coding",dest="no_coding",default=False,action="store_true",help="Do not select the coding genes in multimapper filtering")
	group_filter.add_argument("--no-three-prime",dest="no_three_prime",default=False,action="store_true",help="Do not select 3' end single end reads in paired end mapping")
	
	group_steps = parser.add_argument_group("Steps")
	group_steps.add_argument("--continue",dest="continue_file",help="If given, process the sam file. If not given, bowtie2_[pid].sam")
	group_steps.add_argument("--skip-convert",dest="skip_convert",default=False,action="store_true",help="Skip converting fastq")
	group_steps.add_argument("--skip-mapping",dest="skip_mapping",default=False,action="store_true",help="Skip mapping")
	group_steps.add_argument("--sorted",dest="sorted_bam",default=False,action="store_true",help="Sort bam")
	group_steps.add_argument("--no-index",dest="no_sorted_bam_index",default=False,action="store_true",help="Do not index bam")
	
	group_files = parser.add_argument_group("Files")
	group_files.add_argument("--del-convert",dest="del_convert",default=False,action="store_true",help="Do ot delete the converted fastq files")
	group_files.add_argument("--del-sam",dest="del_sam",default=False,action="store_true",help="Delete the sam files")
	group_files.add_argument("--del-bam",dest="del_bam",default=False,action="store_true",help="Delete the unsorted bam")
	group_files.add_argument("--no-multi-fastq",dest="no_multi_fastq",default=False,action="store_true",help="No multimapping in unmapped fastq")

	group_other = parser.add_argument_group("Other")
	group_other.add_argument("--version",action="version",version="%(prog)s 1.0")
	options = parser.parse_args()
	
	signal.signal(signal.SIGINT,signal_handler)
	
	if options.continue_file:
		bowtie2_name = options.continue_file
	else:
		bowtie2_name = "bowtie2_" + str(os.getpid()) + ".sam"

	#Statistics
	stats = { \
	'fully_unique_PE':0,
	'unique_isoform_PE':0,
	'longest_isoform_PE':0,
	'coding_PE':0,
	'coding_isoform_longest_PE':0,
	'three_prime' : 0,
	'three_prime_fully_unique' : 0,
	'fully_unique_SE':0,
	'unique_isoform_SE':0,
	'longest_isoform_SE':0,
	'coding_SE':0,
	'coding_isoform_longest_SE':0,
	'diff_gene':0,
	'diff_isoform':0,
	'multi_location':0,
	'AS_socre_too_low':0,
	'cannot_distinguish':0,
	}
	#overall
	multimapper_num = 0
	unmapped_num = 0
	unique_num = 0
	C2T_num = 0
	PE_num = 0
	SE_num = 0
	three_prime_num = 0
	#Statistics
	fwd_c2t = None
	rev_g2a = None
	if options.fwd:
		if options.fwd.endswith(".fastq") == False:
			raise ValueError("Forward file do not end with .fastq")
		fwd_c2t = options.fwd.replace(".fastq",".c2t.fastq")
	else:
		raise Warning("Please provide forward read.")
	if options.rev:
		if options.rev.endswith(".fastq") == False:
			raise ValueError("Reverse file do not end with .fastq")
		rev_g2a = options.rev.replace(".fastq",".g2a.fastq")
	else:
		rev_g2a = None
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
		sys.stderr.write("[%s]Mapping with bowtie2, TEMP name: %s\n" % (strftime("%Y-%m-%d %H:%M:%S", time.localtime()),bowtie2_name))
		bowtie2 = options.bowtie2_path
		bowtie2_parameters = read_parameters(options.bowtie2_param)
		bowtie2_parameters = mapping_string(bowtie2_parameters,fwd_c2t,rev_g2a)
		bowtie2_index = options.index
		
		if options.fwd:
			mapping(bowtie2_index,"C2T",bowtie2,bowtie2_parameters)
	else:
		if not options.continue_file:
			raise Warning("Please provide a SAM.")

	if options.skip_convert:
		#Step get C-contained reads
		sys.stderr.write("[%s]Finding non-converted C containing reads...\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
		fwd_read_dict,rev_read_dict,total_reads = read_C_position(options.fwd,options.rev)
	
	if options.del_convert == True:
		os.remove(fwd_c2t)
		if rev_g2a:
			os.remove(rev_g2a)
	
	transcripts_info,gene_biotype = read_transcripts_info(options.genelist)
	
	#Step Get unique mapping and report unmapped and multimappers
	sys.stderr.write("[%s]Handling SAM outputs...\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
	if options.fwd and options.rev:
		sam_handler(bowtie2_name,mode="PE")
	elif options.fwd and not options.rev:
		sam_handler(bowtie2_name,mode="SE")
	
	unique_num = C2T_num
	total_reads = unique_num + multimapper_num + unmapped_num
	sys.stderr.write("[%s]Completed successfully:\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
	sys.stderr.write(" Total reads: %d\n" % total_reads)
	sys.stderr.write(" Unique mapping: %d (%.3f%%), in which\n" % (unique_num,100*unique_num/(total_reads+0.0)))
	sys.stderr.write("   Paired end: %d (%.3f%%)\n" % (PE_num,100*PE_num/(total_reads+0.0)))
	sys.stderr.write("     Fully unique:                       %d (%.3f%%)\n" % (stats["fully_unique_PE"],100*stats["fully_unique_PE"]/(total_reads+0.0)))
	sys.stderr.write("     Unqiue gene (unique isoform):       %d (%.3f%%)\n" % (stats["unique_isoform_PE"],100*stats["unique_isoform_PE"]/(total_reads+0.0)))
	sys.stderr.write("     Unique gene (longest isoform):      %d (%.3f%%)\n" % (stats["longest_isoform_PE"],100*stats["longest_isoform_PE"]/(total_reads+0.0)))
	sys.stderr.write("     Coding selected (unique isoform):   %d (%.3f%%)\n" % (stats["coding_PE"],100*stats["coding_PE"]/(total_reads+0.0)))
	sys.stderr.write("     Coding selected (longest isoform):  %d (%.3f%%)\n" % (stats["coding_isoform_longest_PE"],100*stats["coding_isoform_longest_PE"]/(total_reads+0.0)))
	sys.stderr.write("     3' alignment (see SE):              %d (%.3f%%)\n" % (stats["three_prime"],100*stats["three_prime"]/(total_reads+0.0)))
	sys.stderr.write("   Single end: %d (%.3f%%)\n" % (SE_num,100*SE_num/(total_reads+0.0)))
	sys.stderr.write("     Fully unique:                       %d (%.3f%%)\n" % (stats["fully_unique_SE"],100*stats["fully_unique_SE"]/(total_reads+0.0)))
	sys.stderr.write("     Unqiue gene (unique isoform):       %d (%.3f%%)\n" % (stats["unique_isoform_SE"],100*stats["unique_isoform_SE"]/(total_reads+0.0)))
	sys.stderr.write("     Unique gene (longest isoform):      %d (%.3f%%)\n" % (stats["longest_isoform_SE"],100*stats["longest_isoform_SE"]/(total_reads+0.0)))
	sys.stderr.write("     Coding selected (unique isoform):   %d (%.3f%%)\n" % (stats["coding_SE"],100*stats["coding_SE"]/(total_reads+0.0)))
	sys.stderr.write("     Coding selected (longest isoform):  %d (%.3f%%)\n" % (stats["coding_isoform_longest_SE"],100*stats["coding_isoform_longest_SE"]/(total_reads+0.0)))
	sys.stderr.write(" Multiple mapping: %d (%.3f%%)\n" % (multimapper_num,100*multimapper_num/(total_reads+0.0)))
	sys.stderr.write("   Low scores:        %d (%.3f%%)\n" % (stats["AS_socre_too_low"],100*stats["AS_socre_too_low"]/(total_reads+0.0)))
	sys.stderr.write("   Different gene:    %d (%.3f%%)\n" % (stats["diff_gene"],100*stats["diff_gene"]/(total_reads+0.0)))
	sys.stderr.write("   Different isoform: %d (%.3f%%)\n" % (stats['diff_isoform'],100*stats['diff_isoform']/(total_reads+0.0)))
	sys.stderr.write("   Same isoform:      %d (%.3f%%)\n" % (stats['multi_location'],100*stats['multi_location']/(total_reads+0.0)))
	sys.stderr.write(" Unmapped: %d (%.3f%%)\n\n" % (unmapped_num,100*unmapped_num/(total_reads+0.0)))
	
	if options.sorted_bam == True:
		sys.stderr.write("[%s]Sorting bam...\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
		pysam.sort("-o", options.output+".sorted.bam", options.output+".bam")
		if options.no_sorted_bam_index == False:
			sys.stderr.write("[%s]Indexing bam...\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
			pysam.index(options.output+".sorted.bam")
	if options.del_sam == True:
		os.remove(bowtie2_name)
	if options.del_bam == True:
		os.remove(options.output+".bam")