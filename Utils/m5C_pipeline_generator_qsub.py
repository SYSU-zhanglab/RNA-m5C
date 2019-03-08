import sys,os

# input format:
# sample	read1	read2

if not sys.argv[1]:
	print "Usage: python %s <input.txt>" % (sys.argv[0])

#metadata
ref_genome_index = check_env("</path/to/hisat2_index/>",is_path=True)
ref_genome = check_env("</path/to/genome_fasta>")
ref_mRNA = check_env("</path/to/mRNA_fasta>")
ref_c2t = check_env("</path/to/c2t_mRNA_fasta>")
bowtie2_index = check_env("</path/to/bowtie2_index>",is_prefix=True)
anno = check_env("</path/to/anno file>")
gtf = check_env("</path/to/gtf>")
reg = check_env("</path/to/reg file>")
genelist = check_env('</path/to/genelist file>')
db = check_env( "</path/to/non-redundance base file>")
size = check_env( "</path/to/size file>")

#environment
python = check_env('</path/to/python>')
bowtie2 = check_env('</path/to/bowtie2>')
bowtie2_path = check_env('/path/to/bowtie2/',is_path=True)
samtools = check_env('/path/to/samtools')
cutadapt = check_env('/path/to/cutadapt')
trimmomatic = check_env('/path/to/trimmomatic')
java = check_env('/path/to/java')
hisat2_path = check_env('/path/to/hisat2/',is_path=True)

PATH = "./"

with open(sys.argv[1],'r') as input:
	for line in input.readlines():
		name,read1,read2 = line.strip().split("\t")
		
		qsub_file = name + ".qsub_file"
		with file(qsub_file,"w") as output:
			workpath = PATH + "/" + name + "/"
			
			# headers
			output.write("-N %s \n" % name)
			output.write("-d %s \n" % workpath)
			output.write("-o %s \n" % workpath)

			if os.path.isdir(workpath) == False:
					os.mkdir(workpath)
			
			# cutadapt
			output.write("{cutadapt} -a NNNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A NNNNNNAGATCGGAAGAGCGTCGTGTAGGGAAAGAG -e 0.25 -q 25 --trim-n -o read1.cutadapt.fastq -p read2.cutadapt.fastq ../{read1} ../{read2} \n".format(cutadapt=cutadapt,name=name,read1=read1,read2=read2))
			
			# trimmomatic
			output.write("{java} -jar {trimmomatic} PE -threads 2 -phred33 read1.cutadapt.fastq read2.cutadapt.fastq rev.fastq rev.UP.fastq fwd.fastq fwd.UP.fastq HEADCROP:10 SLIDINGWINDOW:4:22 AVGQUAL:25 MINLEN:40 \n".format(java=java,trimmomatic=trimmomatic))
			
			# remove cutadapt fastq
			output.write("rm read1.cutadapt.fastq read2.cutadapt.fastq \n")
			
			# hisat2 to genome
			output.write("{python} /share/public1/data/liujh/script/m5c/step_by_step/BSpipe/BS_hisat2.py -F fwd.fastq -R rev.fastq -o hisat2_genome -I {ref_genome_index} --index-prefix HISAT2 --hisat2-path {hisat2_path} --del-convert --del-sam \n".format(python=python,hisat2_path=hisat2_path,ref_genome_index=ref_genome_index))

			# bowtie2 to transcriptome
			output.write("{python} /share/public1/data/liujh/script/m5c/step_by_step/BSpipe/BS_bowtie2.py -F fwd.unmapped.fastq -R rev.unmapped.fastq -o bowtie2_trans --bowtie2-path {bowtie2_path} -I {bowtie2_index} -g {genelist} --del-convert --del-sam \n".format(python=python,bowtie2_path=bowtie2_path,bowtie2_index=bowtie2_index,genelist=genelist))

			# transcriptome to genome
			output.write("{python} /share/public1/data/liujh/script/m5c/step_by_step/BSpipe/Bam_transcriptome_to_genome_v1.0.py -i bowtie2_trans.bam -o bowtie2_genome.bam -a {anno} --dict {size} \n".format(python=python,anno=anno,size=size))
			
			# merge BAM
			output.write("{python} /share/public1/data/liujh/script/m5c/step_by_step/BSpipe/concat_bam.py -t 4 -o genome_merged.bam --sort --index -i hisat2_genome.bam bowtie2_genome.bam \n".format(python=python))
			
			# statistics
			output.write("{samtools} flagstat genome_merged.sorted.bam \n".format(samtools=samtools))

			# pileup
			output.write("{python} /share/public1/data/liujh/script/m5c/step_by_step/5_m5C_step-by-step_pileup/pileup_genome_multiprocessing_v1.4.py -P 8 -i genome_merged.sorted.bam -o {name}.pileups.tmp -f {ref_genome} \n".format(python=python,name=name,ref_genome=ref_genome))

			# format pileup
			output.write("{python} /share/public1/data/liujh/script/m5c/step_by_step/5_m5C_step-by-step_pileup/m5C_caller_temp_filter.in_memory.mm.py -i {name}.pileups.tmp -o {name}.pileups.txt --CR {name}.pileups.CR --db {db} \n".format(python=python,name=name,db=db))

		#os.system("qsub %s" % qsub_file)