import sys
from Bio import SeqIO

for seq in SeqIO.parse(sys.argv[1],'fasta'):
	print ">"+seq.id
	print seq.seq
