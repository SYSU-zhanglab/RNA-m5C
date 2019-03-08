#!bin/usr/env python

#Jianheng Liu @ Zhanglab, SYSU
#Feb, 2018
#Email: liujh26@mail2.sysu.edu.cn
#Usage: This powerful program can convert BAMs in other reference to genomic coordinates, as long as the annotation table is given.
#Input: [.bam] [.anno] [.size]

import os,sys,time,signal,argparse,pysam,multiprocessing
from Bio import SeqIO
from Bio.Seq import reverse_complement
from time import gmtime, strftime
from collections import defaultdict
from collections import OrderedDict
from pysam import qualities_to_qualitystring
import subprocess
import copy

def read_anno(fn):
	output = defaultdict(dict)
	with open(fn,'r') as input:
		line = input.readline()
		while (line):
			line = line.strip().split("\t")
			trans_id = line[1]
			chr = line[2]
			dir = line[3]
			exonCounts = int(line[8])
			#All to 0-based
			exonStarts = map(lambda x:int(x),line[9].split(",")[:-1]) #Starts are 0-based 
			exonEnds = map(lambda x:int(x)-1, line[10].split(",")[:-1]) #Ends are 1-based
			gene_id = line[12]
			bins = zip(exonStarts,exonEnds)
			enst_start = 0
			
			output[trans_id]['dir'] = dir
			output[trans_id]['ensg'] = gene_id
			output[trans_id]['chr'] = chr
			output[trans_id]['introns'] = OrderedDict()
			output[trans_id]['exons'] = OrderedDict()
			
			last_end = None
			last_start = None
			enst_start = -1 #0-based
			if dir == "+":
				for start,end in bins:
					enst_start += 1
					enst_end = enst_start + end - start
					output[trans_id]['exons'][(enst_start, enst_end)] = (start,end)
					enst_start = enst_end
					last_end = end
				
			elif dir == "-":
				bins = bins[::-1]
				for start,end in bins:
					enst_start += 1
					enst_end = enst_start + end - start
					output[trans_id]['exons'][(enst_end, enst_start)] = (start,end) #Noted that '-' strand is reverse
					enst_start = enst_end
					last_start = start
			line = input.readline()
	return output
	
def cal(cigar):
	c = {0:0,1:0,2:0,3:0,4:0}
	for a,b in cigar:
		c[a] += b
	return c[0],"M"+str(c[0])+"N"+str(c[3])+"I"+str(c[1])+"D"+str(c[2])
	
def generate_new_cigar(all_bins,start,end,old_cigar,trans_dir):
	''' order: small --> big corrdinate '''
	new_cigar_tmp = [] #no del and insert, with intron
	if trans_dir == "-":
		old_cigar = old_cigar[::-1]
		all_bins = all_bins[::-1]
		start,end = end,start
	# print all_bins
	all_bins_iter = iter(all_bins)
	while (1):
		try:
			x,y = next(all_bins_iter)
			# print x,y,start,end
			if x <= start <= y <end:
				new_cigar_tmp.append([0,y-start+1])
				exon_edge = y
			elif x <= start <= end <= y:
				new_cigar_tmp.append([0,end - start +1])
				break
			elif start < x <= y < end:
				if x-exon_edge-1 > 0:
					new_cigar_tmp.append([3,x-exon_edge-1])
				new_cigar_tmp.append([0,y-x+1])
				exon_edge = y
			elif start < x <= end <= y:
				if x-exon_edge-1 > 0:
					new_cigar_tmp.append([3,x-exon_edge-1])
				new_cigar_tmp.append([0,end-x+1])
				break
		except StopIteration:
			print x,y,start,end
			sys.exit()
	
	new_cigar_tmp_tmp = []

	new_cigar_tmp_iter = iter(new_cigar_tmp)
	cigar_type,number = next(new_cigar_tmp_iter)
	while (1):
		try:
			cigar_type_1,number_1 = next(new_cigar_tmp_iter)
			if cigar_type == cigar_type_1:
				number = number + number_1
			else:
				new_cigar_tmp_tmp.append([cigar_type,number])
				cigar_type,number = cigar_type_1,number_1
		except StopIteration:
			new_cigar_tmp_tmp.append([cigar_type,number])
			break
	new_cigar_tmp = new_cigar_tmp_tmp
	# print new_cigar_tmp
	new_cigar = []
	#debug
	
	# print "*" * 30
	# print trans_dir
	# old_M, old = cal(old_cigar)
	# print old_cigar,old_M,old
	# NT_M, NT = cal(new_cigar_tmp)
	# print new_cigar_tmp,NT_M,NT
	
	new_cigar_tmp_iter = iter(new_cigar_tmp)
	block = next(new_cigar_tmp_iter)
	for cigar_type,num in old_cigar:
		try:
			if block[0] == 3:
				new_cigar.append((block[0],block[1]))
				block = next(new_cigar_tmp_iter)
			if cigar_type == 0: #matched
				if num < block[1]: #samller than the original block
					new_cigar.append((0,num))
					block[1] = block[1] - num
				elif num == block[1]: #remove a block
					new_cigar.append((0,num))
					block = next(new_cigar_tmp_iter)
					if block[0] == 3: #intron
						new_cigar.append((block[0],block[1]))
						block = next(new_cigar_tmp_iter)
				else:
					while num > block[1]:
						new_cigar.append((0,block[1]))
						num = num - block[1]
						block = next(new_cigar_tmp_iter)
						new_cigar.append((block[0],block[1])) #intron
						block = next(new_cigar_tmp_iter) #until the last exon
					if num == block[1]:
						new_cigar.append((0,num))
						block = next(new_cigar_tmp_iter)
					elif num < block[1]:
						block[1] = block[1] - num
						new_cigar.append((0,num))
					if block[1] < 0:
						raise Warning("Block start <0")
			elif cigar_type == 1: #insert
				new_cigar.append((1,num))
			elif cigar_type == 2: #del
				if num < block[1]:
					new_cigar.append((2,num))
					block[1] = block[1] - num
				elif num == block[1]:
					new_cigar.append((2,num))
					block = next(new_cigar_tmp_iter)
					if block[0] == 3:
						new_cigar.append((block[0],block[1]))
						block = next(new_cigar_tmp_iter)
				else:
					while num > block[1]:
						new_cigar.append((2,block[1]))
						num = num - block[1]
						block = next(new_cigar_tmp_iter)
						new_cigar.append((block[0],block[1])) #intron
						block = next(new_cigar_tmp_iter) #until the last exon
					if num == block[1]:
						new_cigar.append((2,num))
						block = next(new_cigar_tmp_iter)
					elif num < block[1]:
						block[1] = block[1] - num
						new_cigar.append((2,num))
					if block[1] < 0:
						raise Warning("Block start <0")
			elif cigar_type == 3:
				new_cigar.append((3,num))
			elif cigar_type == 4:
				new_cigar.append((4,num))
			elif cigar_type == 5:
				new_cigar.append((5,num))
			elif cigar_type == 6:
				new_cigar.append((6,num))
		except StopIteration:
			continue
	# new_M,new = cal(new_cigar)
	# print new_cigar,new_M,new
	# if new_M != old_M:
		# print False
		# print start,end,all_bins
	
	#debug
	
	return new_cigar

def map_to_genome(segment):
	global UNLIFT,total,lifted,unlifted
	total += 1
	try:
		genome_info = annotation.get(segment.reference_name)
		new_ref_id = header_dict.get(genome_info['chr'])
		trans_dir = genome_info['dir']
	except TypeError:
		genome_info = None
		new_ref_id = None
	# if genome_info['chr'] not in header_dict:
		# print genome_info['chr']
	# if not genome_info or not new_ref_id:
		# print segment.reference_name,genome_info['chr'],header_dict.get(genome_info['chr'])
	if genome_info and new_ref_id is not None:
		lifted += 1
		old_start = segment.reference_start #0-based
		old_end = segment.reference_end - 1 #1-based
		new_start = None
		new_end = None
		# print "*"*20
		# print segment.reference_name,old_start,old_end
		if trans_dir == "+":
			genome_info_iter = iter(genome_info["exons"].items())
		elif trans_dir == "-":
			genome_info_iter = iter(genome_info["exons"].items()[::-1])

		while new_start is None or new_end is None:
			key,values = genome_info_iter.next()
			start,end = key
			geno_start,geno_end = values
			# print start,end,geno_start,geno_end
			if trans_dir == "+":
				# print start,end,geno_start,geno_end
				if start <= old_start <= end:
					new_start = geno_start + old_start - start
				if start <= old_end <= end:
					new_end = geno_start + old_end - start
			elif trans_dir == "-":
				start,end = end,start
				# print start,end,geno_start,geno_end
				# if start <= old_end <= end:
					# new_end = geno_end - (old_end - start) 
				# if start <= old_start <= end:
					# new_start = geno_end - (old_start - start)
				if start <= old_end <= end:
					new_end = geno_start + (end - old_end)
					# print old_end,new_end
					# print old_start,new_start,geno_end,start,end
				if start <= old_start <= end:
					new_start = geno_start + (end - old_start)
					# print old_start,new_start
			else:
				raise Warning("Transcription direction loss.")

		new_cigar = generate_new_cigar(genome_info["exons"].values(),new_start,new_end,segment.cigar,genome_info['dir'])
		
		qual = segment.query_qualities
		mpq = segment.mapping_quality
		seq = segment.query_sequence
		
		segment_output = pysam.AlignedSegment()
		segment_output.tags = segment.tags
		if trans_dir == "-":
			new_start,new_end = new_end,new_start
			qual = qual[::-1]
			seq = reverse_complement(segment.query_sequence)
			if segment.is_reverse:
				segment.is_reverse = False
				segment.mate_is_reverse = True
			else:
				segment.is_reverse = True
				segment.mate_is_reverse = False
			segment_output.set_tag("TS","-")
			# segment_output.set_tag("YG","G2A")
		else:
			segment_output.set_tag("TS","+")
				
		segment_output.query_name = segment.query_name
		segment_output.flag = segment.flag
		segment_output.reference_id = new_ref_id
		segment_output.reference_start = new_start
		segment_output.cigar = new_cigar
		segment_output.query_sequence = seq
		segment_output.query_qualities = segment.query_qualities
		segment_output.mapping_quality = mpq
		segment_output.set_tag("GN",genome_info["ensg"])
		segment_output.set_tag("TN",segment.reference_name)
		segment_output.set_tag("TP",segment.reference_start+1) #1-based
		if segment_output:
			# print segment.reference_start,segment.reference_end
			# print segment_output.reference_start,segment_output.reference_end
			if options.verify == True:
				cigar_length = 0
				for cigar_type,num in segment_output.cigar:
					if cigar_type == 0 or cigar_type == 1:
						cigar_length += num
				if cigar_length != len(segment_output.query_sequence):
					print segment
					print segment_output
					raise ValueError("Cigar != sequence")
			return segment_output
		else:
			print segment.reference_name,segment.cigarstring,segment.reference_start,segment.reference_end,genome_info["chr"],genome_info["dir"],segment_output.cigarstring,segment_output.reference_start,segment_output.reference_end
	
	else:
		unlifted += 1
		if options.no_unlift == False:
			UNLIFT.write(segment)
		
if __name__ == "__main__":
	description = """
	"""
	parser = argparse.ArgumentParser(prog="m5C_mapper",fromfile_prefix_chars='@',description=description,formatter_class=argparse.RawTextHelpFormatter)
	#Require
	group_required = parser.add_argument_group("Required")
	group_required.add_argument("--input","-i",dest="input",required=True,help="input SAM/BAM, default is BAM")
	group_required.add_argument("--output","-o",dest="output",required=True,help="output SAM/BAM, default is BAM")
	group_required.add_argument("--anno","-a",dest="anno",required=True,help="UCSC-all table like annotation")
	#Filter
	group_optional = parser.add_argument_group("Optional")
	group_required.add_argument("--ref","-r",dest="ref",required=False,help="Reference for new header")
	group_required.add_argument("--fasta","-f",dest="fasta",required=False,help="Reference for new header")
	group_required.add_argument("--dict","-d",dest="dict",required=False,help="Reference dict for new header")
	group_required.add_argument("--header","-H",dest="header",required=False,help="Header file for new header")
	group_optional.add_argument("--no-unlift",dest="no_unlift",default=False,action="store_true",help="Do not report unlifted sequences")
	# group_optional.add_argument("--unlift-to-unmapped",dest="unlift_to_unmapped",default=False,action="store_true",help="Set unlift reads unmapped")
	group_optional.add_argument("--sort",dest="sort",default=False,action="store_true",help="Sort bam (and delete unsort)")
	group_optional.add_argument("--no-del-bam",dest="no_del_bam",default=False,action="store_true",help="Do not del bam file after sorting")
	group_optional.add_argument("--index",dest="index",default=False,action="store_true",help="Index sorted bam")
	group_optional.add_argument("--verify",dest="verify",default=False,action="store_true",help="Check output cigar and length")
	group_other = parser.add_argument_group("Other")
	group_other.add_argument("--version",action="version",version="%(prog)s 1.0")
	options = parser.parse_args()
	
	sys.stderr.write("[%s]Loading annotations...\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
	annotation = read_anno(options.anno)
	
	total = 0
	lifted = 0
	unlifted = 0
	
	in_mode = "rb"
	out_mode = "wb"

	sys.stderr.write("[%s]Buidling genome header...\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
	header = {}
	header['HD'] = {'SO': 'unsorted', 'VN': '1.0'}
	header['SQ'] = []
	if options.fasta:
		with file("ref.dict",'w') as output:
			for seq in SeqIO.parse(options.fasta,'fasta'):
				length = len(seq.seq)
				header['SQ'].append({"SN":seq.id,"LN":length})
				output.write(seq.id+"\t"+str(length)+"\n")
	elif options.dict:
		with open(options.dict) as input:
			for line in input.readlines():
				line = line.strip().split("\t")
				header['SQ'].append({"SN":line[0],"LN":int(line[1])})
	elif options.ref:
		with pysam.AlignmentFile(options.ref, 'rb') as input:
			header = input.header
	else:
		raise Warning("Please provide a reference or a header")

	sys.stderr.write("[%s]Coordinate conversion starts\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
	
	if options.no_unlift == False:
		unlift_fn = options.input.replace(".bam",".unlift.bam")
		with pysam.AlignmentFile(options.input, in_mode) as INPUT, pysam.AlignmentFile(options.output, out_mode, header = header) as OUTPUT, pysam.AlignmentFile(unlift_fn, 'wb', template=INPUT) as UNLIFT:
			header_dict = {}
			n_header = 0
			for header in OUTPUT.header["SQ"]:
				header_dict[header['SN']] = n_header
				n_header += 1
			# print len(header_dict)
			# sys.exit()
			for segment in INPUT:
				segment_output = map_to_genome(segment)
				if segment_output:
					OUTPUT.write(segment_output)
	else:
		with pysam.AlignmentFile(options.input, in_mode) as INPUT, pysam.AlignmentFile(options.output, out_mode, header = header) as OUTPUT:
			header_dict = {}
			n_header = 0
			for header in OUTPUT.header["SQ"]:
				header_dict[header['SN']] = n_header
				n_header += 1
			# print len(header_dict)
			for segment in INPUT:
				segment_output = map_to_genome(segment)
				if segment_output:
					OUTPUT.write(segment_output)
	sys.stderr.write("[%s]Finished.\n  Total: %d\n  Lifted: %d\n  Unlifted: %d\n\n" % (strftime("%Y-%m-%d %H:%M:%S", time.localtime()),total,lifted,unlifted))
	
	if options.sort == True:
		sys.stderr.write("[%s]Sorting bam...\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
		pysam.sort("-o", options.output.replace(".bam",".sorted.bam"), options.output)
		if options.index == True:
			sys.stderr.write("[%s]Indexing bam...\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
			pysam.index(options.output.replace(".bam",".sorted.bam"))
		if options.no_del_bam == False:
			os.remove(options.output)

