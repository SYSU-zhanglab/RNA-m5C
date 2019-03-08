import sys,os
from sjm_tools import job,check_env

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

		SJM = name + ".sjm"
		workpath = PATH + "/" + name + "/"

		if os.path.isdir(workpath) == False:
				os.mkdir(workpath)

		JOB = job(workpath,SJM)
		SJM = name + ".sjm"
		workpath = PATH + "/" + name + "/"
		
		if os.path.isdir(workpath) == False:
			os.mkdir(workpath)
			
		JOB = job(workpath,SJM) 
		
		JOB.step_start(step_name="QC",memory="50G")
		JOB.add_process("{cutadapt} -a NNNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A NNNNNNAGATCGGAAGAGCGTCGTGTAGGGAAAGAG -e 0.25 -q 25 --trim-n -o read1.cutadapt.fastq -p read2.cutadapt.fastq ../{read1} ../{read2}".format(cutadapt=cutadapt,name=name,read1=read1,read2=read2))
		JOB.add_process("{java} -jar {trimmomatic} PE -threads 2 -phred33 read1.cutadapt.fastq read2.cutadapt.fastq rev.fastq rev.UP.fastq fwd.fastq fwd.UP.fastq HEADCROP:10 SLIDINGWINDOW:4:22 AVGQUAL:25 MINLEN:40".format(java=java,trimmomatic=trimmomatic))
		JOB.add_process("rm read1.cutadapt.fastq read2.cutadapt.fastq")
		JOB.step_end()
		
		
		JOB.step_start(step_name="hisat2",memory="80G")
		JOB.add_process("{python} /share/public1/data/liujh/script/m5c/step_by_step/BSpipe/BS_hisat2.py -F fwd.fastq -R rev.fastq -o hisat2_genome -I {ref_genome_index} --index-prefix HISAT2 --hisat2-path {hisat2_path} --del-convert --del-sam".format(python=python,hisat2_path=hisat2_path,ref_genome_index=ref_genome_index))
		JOB.step_end()
		
		JOB.step_start(step_name="bowtie2",memory="80G")
		JOB.add_process("{python} /share/public1/data/liujh/script/m5c/step_by_step/BSpipe/BS_bowtie2.py -F fwd.unmapped.fastq -R rev.unmapped.fastq -o bowtie2_trans --bowtie2-path {bowtie2_path} -I {bowtie2_index} -g {genelist} --del-convert --del-sam".format(python=python,bowtie2_path=bowtie2_path,bowtie2_index=bowtie2_index,genelist=genelist))
		JOB.step_end()
		
		JOB.step_start(step_name="LiftOver",memory="80G")
		JOB.add_process("{python} /share/public1/data/liujh/script/m5c/step_by_step/BSpipe/Bam_transcriptome_to_genome_v1.0.py -i bowtie2_trans.bam -o bowtie2_genome.bam -a {anno} --dict {size}".format(python=python,anno=anno,size=size))
		JOB.step_end()
		
		JOB.step_start(step_name="Merge",memory="80G")
		JOB.add_process("{python} /share/public1/data/liujh/script/m5c/step_by_step/BSpipe/concat_bam.py -t 4 -o genome_merged.bam --sort --index -i hisat2_genome.bam bowtie2_genome.bam".format(python=python))
		JOB.add_process("{samtools} flagstat genome_merged.sorted.bam".format(samtools=samtools))
		JOB.step_end()

		JOB.step_start(step_name="Pileup",memory="80G")
		JOB.add_process("{python} /share/public1/data/liujh/script/m5c/step_by_step/5_m5C_step-by-step_pileup/pileup_genome_multiprocessing_v1.4.py -P 8 -i genome_merged.sorted.bam -o {name}.pileups.tmp -f {ref_genome}".format(python=python,name=name,ref_genome=ref_genome))
		JOB.step_end()

		JOB.step_start(step_name="Call",memory="80G")
		JOB.add_process("{python} /share/public1/data/liujh/script/m5c/step_by_step/5_m5C_step-by-step_pileup/m5C_caller_temp_filter.in_memory.mm.py -i {name}.pileups.tmp -o {name}.pileups.txt --CR {name}.pileups.CR --db {db}".format(python=python,name=name,db=db))
		JOB.step_end()
		
		JOB.job_finish()
		#JOB.job_submit()