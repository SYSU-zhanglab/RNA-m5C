#!bin/usr/env python

#Jianheng Liu @ Zhanglab, SYSU
#Feb, 2018
#Email: liujh26@mail2.sysu.edu.cn
#Usage: Create indexes (both C2T and G2A) for HISAT2
#Input: [.fa] [.gtf] [path to HISAT2]

__doc__ = """
Basic usage:
python BS_hisat2_index.py \
-i references.fa \
-o output_dir \
--gtf references.gtf \
--control-fasta control.fa \
--hisat2-path /path/to/hisat2/

Build index for hisat2 in genome_trans mode (splice sites 
and exons are extracted from gtf annotation). Output prefix
should be HISAT2.

Note that if spike ins were added, this script will append
them to the C2T genome as default, since mostly they were
treated as transcript references.

If you have chromosomes that should be mapped at certain 
directions, you can also treat them as controls. But please
make sure that the order of common references between C2T
and G2A index should be the same, or a chaos will happen 
in handling hisat2 mapping results.
"""
import os,sys,time,signal,argparse,pysam,multiprocessing
from Bio import SeqIO
from Bio.Seq import reverse_complement
from time import gmtime, strftime
from collections import defaultdict
from pysam import qualities_to_qualitystring
import subprocess

def fasta_c2t(output):
	with file(output,'w') as output:
		for seq_record in SeqIO.parse(options.input,"fasta"):
			output.write("".join([">",seq_record.description,"\n"]))
			sequence = str(seq_record.seq)
			sequence = sequence.replace("C","T")
			sequence = sequence.replace("c","t")
			output.write("".join([sequence,"\n"]))
		if options.control:
			for seq_record in SeqIO.parse(options.control,"fasta"):
				output.write("".join([">",seq_record.description,"\n"]))
				sequence = str(seq_record.seq)
				sequence = sequence.replace("C","T")
				sequence = sequence.replace("c","t")
				output.write("".join([sequence,"\n"]))
		
def fasta_g2a(output):
	with file(output,'w') as output:
		for seq_record in SeqIO.parse(options.input,"fasta"):
			output.write("".join([">",seq_record.description,"\n"]))
			sequence = str(seq_record.seq)
			sequence = sequence.replace("G","A")
			sequence = sequence.replace("g","a")
			output.write("".join([sequence,"\n"]))
		if options.control and options.control_G2A == True:
			for seq_record in SeqIO.parse(options.control,"fasta"):
				output.write("".join([">",seq_record.description,"\n"]))
				sequence = str(seq_record.seq)
				sequence = sequence.replace("G","A")
				sequence = sequence.replace("g","a")
				output.write("".join([sequence,"\n"]))

def index(mode,hisat2_build,genome_name,ss,exon):
	if mode == "C2T":
		out_dir = options.output+"/C2T"
		fasta_name = options.output+"/C2T/HISAT2_C2T.fa"
		index_base = options.output+"/C2T/HISAT2_C2T"
		fasta_c2t(fasta_name)
	elif mode == "G2A":
		out_dir = options.output+"/G2A/"
		fasta_name = options.output+"/G2A/HISAT2_G2A.fa"
		index_base = options.output+"/G2A/HISAT2_G2A"
		fasta_g2a(fasta_name)
	cmd = [hisat2_build,"--ss",ss,"--exon",exon,"-p",str(options.threads),fasta_name,index_base]
	hisat2_out = subprocess.Popen(cmd,\
								  stdin=subprocess.PIPE,\
								  stdout=subprocess.PIPE,\
								  stderr=subprocess.PIPE,\
								  shell=False)
	stdout,stderr = hisat2_out.communicate()
	# print stdout
	sys.stderr.write("[%s]%s report:\n" % (strftime("%Y-%m-%d %H:%M:%S", time.localtime()),mode))
	sys.stderr.write(stdout)
	sys.stderr.write(stderr)
	sys.stderr.write("\n")
	
def signal_handler(sig,frame):
	pool.terminate()
	sys.exit()

if __name__ == "__main__":
	description = __doc__
	parser = argparse.ArgumentParser(prog="BS_hisat2_index",fromfile_prefix_chars='@',description=description,formatter_class=argparse.RawTextHelpFormatter)
	#Require
	group_required = parser.add_argument_group("Required")
	group_required.add_argument("-i","--input",dest="input",required=True,help="Input fasta")
	group_required.add_argument("-o","--output",dest="output",required=True,help="Output dir")
	group_required.add_argument("--gtf",dest="gtf",required=False,help="Gtf file, not needed if skip ss and exon extraction")
	group_optional = parser.add_argument_group("Optional")
	group_optional.add_argument("--control-fasta",dest="control",help="Control fasta file, only append to C2T index")
	group_optional.add_argument("--control-G2A",dest="control_G2A",default=False,action="store_true",help="Create G2A control index as well")
	group_optional.add_argument("-p","--threads",dest="threads",type=int,default=4,help="Threads for SINGLE hisat2-build, default=4")
	group_optional.add_argument("--I-have-big-memory",dest="sim",default=False,action="store_true",help="I have memory big enough to run two hisa2-build at the same time LOL!!!\n[Warning] This option needs >250G memory for hg19")
	group_optional.add_argument("--C2T-only",dest="C2T_only",default=False,action="store_true",help="Only create C2T index")
	group_optional.add_argument("--G2A-only",dest="G2A_only",default=False,action="store_true",help="Only create G2A index")
	group_optional.add_argument("--skip-ss",dest="skip_ss",default=False,action="store_true",help="Skip generate .ss file")
	group_optional.add_argument("--skip-exon",dest="skip_exon",default=False,action="store_true",help="Skip generate .exon file")
	group_optional.add_argument("--hisat2-path",dest="hisat2_path",help="Path to hisat2 bin")

	
	group_other = parser.add_argument_group("Other")
	group_other.add_argument("--version",action="version",version="%(prog)s 1.0")
	options = parser.parse_args()
	
	sys.stderr.write("[%s]Build index, cmd: %s\n" % (strftime("%Y-%m-%d %H:%M:%S", time.localtime())," ".join(sys.argv)))
	
	hisat2_splice = options.hisat2_path + "/extract_splice_sites.py"
	hisat2_exons = options.hisat2_path + "/extract_exons.py"
	hisat2_build = options.hisat2_path + "/hisat2-build"
	genome_name = options.input.split("/")[-1].strip(".fa")
	exon_file = options.output+"/"+genome_name+".exon"
	ss_file = options.output+"/"+genome_name+".ss"
	

	
	if not os.path.exists(options.output):
		os.mkdir(options.output)

	if not os.path.exists(options.output+"/C2T") and options.G2A_only == False:
		os.mkdir(options.output+"/C2T")
	if not os.path.exists(options.output+"/G2A") and options.C2T_only == False:
		os.mkdir(options.output+"/G2A")
	
	if options.skip_exon == False:
		sys.stderr.write("[%s]Extract exons\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
		cmd = "python %s %s > %s" % (hisat2_exons,options.gtf,exon_file)
		print cmd
		os.system(cmd)
	if options.skip_ss == False:
		sys.stderr.write("[%s]Extract splicing sites\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
		cmd = "python %s %s > %s" % (hisat2_splice,options.gtf,ss_file)
		print cmd
		os.system(cmd)
	if options.sim:
		signal.signal(signal.SIGINT,signal_handler)
		pool = multiprocessing.Pool(2)
		try:
			pool.apply_async(index, args=("C2T",hisat2_build,genome_name,ss_file,exon_file,))
			pool.apply_async(index, args=("G2A",hisat2_build,genome_name,ss_file,exon_file,))
			pool.close()
			pool.join()
		finally:
			pool.terminate()
	#NOT ENOUGH MEMORY	
	else:
		if options.G2A_only == False:
			index("C2T",hisat2_build,genome_name,ss_file,exon_file)
		if options.C2T_only == False:
			index("G2A",hisat2_build,genome_name,ss_file,exon_file)
	sys.stderr.write("[%s]Index built!\n\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))