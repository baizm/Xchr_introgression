import sys
import re

in_file = sys.argv[1]
out_file = sys.argv[2]

f = open(in_file, 'r')
o = open(out_file, 'w')

for i in f:
	if not 'chrUn_' in i:
		o.write(i)

f.close()
o.close()

#usage example: python filter_chrUn.py calJac3.fa.ids.txt target_ids.txt
#makes file of ids that do not include unassigned contigs

#then you can make a copy of the genome without the chrUn's
#module load samtools/1.5
#gunzip calJac3.fa.masked.gz
#samtools faidx calJac3.fa.masked
#xargs samtools faidx calJac3.fa.masked < target_ids.txt >> marmoset.fa