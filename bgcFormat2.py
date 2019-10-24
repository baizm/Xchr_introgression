import subprocess
import re
import pandas as pd

fem_ind = '/scratch/lsa_flux/baizm/bgc_ch3/female_ids.txt'
in_vcf = '/scratch/lsa_flux/baizm/snpCalling/filtered_final_4.vcf'
contigs = '/scratch/lsa_flux/baizm/bgc_ch3/all_contigs.txt'
api_ind = '/scratch/lsa_flux/baizm/bgc_ch3/female_ids_api.txt'
apm_ind = '/scratch/lsa_flux/baizm/bgc_ch3/female_ids_apm.txt'
hz_ind = '/scratch/lsa_flux/baizm/bgc_ch3/female_ids_hz.txt'
out_dir = '/scratch/lsa_flux/baizm/bgc_ch3/'

#filter vcf to keep only female individuals
subprocess.call("bcftools view -S %s %s > %sfemales.vcf" % (fem_ind, in_vcf, out_dir), shell=True)

#filter to keep only biallelic SNPs
subprocess.call("bcftools view -m2 -M2 -v snps %sfemales.vcf > %sbiSNPs.vcf" % (out_dir, out_dir), shell=True)

#filter vcf to keep only SNPs on X and autosomal contigs
subprocess.call("vcftools --bed %s --vcf %sbiSNPs.vcf --recode --out %scontigs" % (contigs, out_dir, out_dir), shell=True)

#filter vcf to get rid of monomorphic sites by setting minor allele frequency low
subprocess.call("vcftools --vcf %scontigs.recode.vcf --maf 0.05 --recode --out %smaf" % (out_dir, out_dir), shell=True)

#filter vcf to keep only sites with minimum mean depth of 12 across individuals
subprocess.call("vcftools --vcf %smaf.recode.vcf --min-meanDP 12 --recode --out %sDP12" % (out_dir, out_dir), shell=True)

#subset vcf by populations
subprocess.call("bcftools view -S %s %sDP12.recode.vcf > %shz.vcf" % (hz_ind, out_dir, out_dir), shell=True)
subprocess.call("bcftools view -S %s %sDP12.recode.vcf > %sapi.vcf" % (api_ind, out_dir, out_dir), shell=True)
subprocess.call("bcftools view -S %s %sDP12.recode.vcf > %sapm.vcf" % (apm_ind, out_dir, out_dir), shell=True)

#make vcfs with sites present in at least 14 individuals in each parental pop
subprocess.call("vcftools --vcf %sapm.vcf --max-missing-count 3 --recode --out %sapm_28" % (out_dir, out_dir), shell=True) #exclude sites missing in 3 or more individuals (6 alleles)
subprocess.call("vcftools --vcf %sapi.vcf --max-missing-count 9 --recode --out %sapi_28" % (out_dir, out_dir), shell=True) #exclude sites missing in 9 or more indivdiuals (18 alleles)

#make list of sites in each file so they can be merged dropping sites not present in both
f='%sapm_28.recode.vcf' % (out_dir)
o='%sapm_28_sites.txt' % (out_dir)

f=open(f, 'r')
o=open(o, 'w')

l1=[]

#for i in f:
	if '#' not in i:
		split=re.split('\t', i.strip())
		l1.append(map(str,split[0:2]))

print('N sites in apm_28.recode.vcf: '+str(len(l1)))

f.close()
o.close()

f='%sapi_28.recode.vcf' % (out_dir)
o='%sapi_28_sites.txt' % (out_dir)

f=open(f, 'r')
o=open(o, 'w')

l2=[]

#for i in f:
	if '#' not in i:
                split=re.split('\t', i.strip())
                l2.append(map(str,split[0:2]))

print('N sites in api_28.recode.vcf: '+str(len(l2)))

f.close()
o.close()

l3=[i for i in l1 if i in l2]
print('N sites in both files: '+str(len(l3)))

f='%ssites_both.txt' % (out_dir)

with open(f, 'w') as file:
	for i in l3:
		contig,pos=i[0],i[1]
		file.write(str(contig)+'\t'+str(pos)+'\n')
file.close()

#make file with one random SNP per contig
data=pd.read_csv('sites_both.txt',sep='\t',header=None)
data.columns=['contig','pos']
random=data.groupby('contig',group_keys=False).apply(lambda data: data.sample(1))
random.sort_index(inplace=True)
random.to_csv('random_snps.txt', sep='\t',index=False,header=None)

#filter population vcfs to include only random SNPs in common between api_28 and apm_28
subprocess.call("vcftools --vcf %shz.vcf --positions %srandom_snps.txt --recode --out %shz_filtered" % (out_dir, out_dir, out_dir), shell=True)
subprocess.call("vcftools --vcf %sapi.vcf --positions %srandom_snps.txt --recode --out %sapi_filtered" % (out_dir, out_dir, out_dir), shell=True)
subprocess.call("vcftools --vcf %sapm.vcf --positions %srandom_snps.txt --recode --out %sapm_filtered" % (out_dir, out_dir, out_dir), shell=True)


###get data in format for bgc genotype uncertainty model###
#for hz pop
f='%shz_filtered.recode.vcf' % (out_dir)
o='%sbgc_hz_filtered.txt' % (out_dir)

f=open(f, 'r')
o=open(o, 'w')

n_loci=0

for i in f:
	if '#' not in i:
		n_loci += 1
		o.write('locus '+str(n_loci)+'\npop 0\n')
		split=re.split('\t', i.strip())
		t=split[9:]
		for i in t:
			GT,PL,DP,DPR=i.split(':')
			o.write(DPR.replace(',',' ')+'\n')
			
f.close()
o.close()

##for api pop
f='%sapi_filtered.recode.vcf' % (out_dir)
o='%sbgc_api_filtered.txt' % (out_dir)

f=open(f, 'r')
o=open(o, 'w')

n_loci=0

for i in f:
	if '#' not in i:
		n_loci += 1
		o.write('locus '+str(n_loci)+'\n')
		split=re.split('\t', i.strip())
		t=split[9:]
		for i in t:
			GT,PL,DP,DPR=i.split(':')
			o.write(DPR.replace(',',' ')+'\n')
			
f.close()
o.close()

##for apm pop
f='%sapm_filtered.recode.vcf' % (out_dir)
o='%sbgc_apm_filtered.txt' % (out_dir)

f=open(f, 'r')
o=open(o, 'w')

n_loci=0

for i in f:
	if '#' not in i:
		n_loci += 1
		o.write('locus '+str(n_loci)+'\n')
		split=re.split('\t', i.strip())
		t=split[9:]
		for i in t:
			GT,PL,DP,DPR=i.split(':')
			o.write(DPR.replace(',',' ')+'\n')
			
f.close()
o.close()

#usage: python bgcFormat.py
