import gzip
import sys
from Bio import SeqIO

fasta = gzip.open(sys.argv[1], 'rU')
out_file = sys.argv[2]

o=open(out_file, 'w')

for i in SeqIO.parse(fasta, 'fasta'):
	name = i.id
	o.write(str(name)+'\n')

fasta.close()
o.close()

#usage example: python get_ids.py calJac3.fa.masked.gz calJac3.fa.ids.txt