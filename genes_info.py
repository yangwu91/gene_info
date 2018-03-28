#!/usr/bin/env python

'''Usage:
python gene_info.py <GFF3>
'''

import sys
from interval import interval

try:
	upstream_len = int(sys.argv[2])
except:
	upstream = 1000

def gene_info(intervals):
	# input format: [[gene], [exon1], [exon2], ...]
	# [[1,1000],[1,100],[300,1000]]
	gene = intervals[0]
	exons = interval()
	for i in intervals[1:]:
		exons = exons | interval(i)
	num_intron = len(exons) - 1
	if gene[0] < gene[-1]:
		gene5 = min(gene)
		gene3 = max(gene)
		exons5 = min(min(exons))
		exons3 = max(max(exons))
		if exons5 - gene5 > 1:
			num_intron += 1
		if gene3 - exons3 > 1:
			num_intron += 1
		upstream = [gene5-upstream_len, gene5]
		strand = '+'
	elif gene[0] > gene[-1]:
		gene3 = min(gene)
		gene5 = max(gene)
		exons3 = min(min(exons))
		exons5 = max(max(exons))
		if exons3 - gene3 > 1:
			num_intron += 1
		if gene5 - exons5 > 1:
			num_intron += 1
		upstream = [gene5+1, gene5+upstream_len+1]
		strand = '-'
	if upstream[0] < 0:
		upstream[0] = 0
	gene_len = max(gene) - min(gene)
	return map(str, (num_intron, gene_len, upstream[0], upstream[-1], strand))

def gff_info(gff_data):
	gff_data = gff_data.strip().split(';')
	gff_data_dict = {}
	for g in gff_data:
		g = g.split('=')
		gff_data_dict[g[0]] = g[-1]
	return gff_data_dict

if __name__ == '__main__':
	try:
		gff3 = sys.argv[1]
	except IndexError:
		print __doc__
		sys.exit(1)
	with open(gff3, 'r') as inf:
		gene_region = []
		outf1 = open('%s.gene_structure.tsv' % gff3, 'w')
		outf1.write('Scaffold\tTranscript\t#intron\tLength\n')
		outf2 = open('%s.upstream%s.bed' % (gff3, upstream_len), 'w')
		for line in inf:
			if line[0] != '#':
				line = line.strip().split('\t')
				feature = line[2]
				if feature == 'mRNA':
					if len(gene_region) != 0:
						gene_structure = gene_info(gene_region)
						outf1.write('%s\t%s\t%s\n' % (scfd, gene_id, '\t'.join(gene_structure[:2])))
						if int(gene_structure[3]) - int(gene_structure[2]) >= upstream_len:
							outf2.write('%s\t%s\t%s\t%s\t%s\n' % (scfd, gene_structure[2], gene_structure[3], gene_id, gene_structure[-1]))
					scfd = line[0]
					gene_id = gff_info(line[8])['ID']
					strand = line[6]
					if strand == '+':
						gene_region = [[int(line[3]), int(line[4])], ]
					elif strand == '-':
						gene_region = [[int(line[4]), int(line[3])], ]
				elif feature == 'exon' and gff_info(line[8])['Parent'] == gene_id:
					gene_region.append([int(line[3]), int(line[4])])
		outf1.write('%s\t%s\t%s\n' % (scfd, gene_id, '\t'.join(gene_structure)))
		if int(gene_structure[3]) - int(gene_structure[2]) >= upstream_len:
			outf2.write('%s\t%s\t%s\t%s\t%s\n' % (scfd, gene_structure[2], gene_structure[3], gene_id, gene_structure[-1]))
	outf1.close()
	outf2.close()
