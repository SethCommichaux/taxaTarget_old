import argparse
import sys
import math
import pandas as pd
import numpy as np
from collections import Counter
from scipy import stats

def read_file_info():
	for h,i in enumerate(open('read_file_info.txt')):
		if h == 0: number_reads = float(i.strip())
		elif h == 1: read_len = float(i.strip())
	return number_reads,read_len

def taxa_mappings(): # function for parsing the NCBI fullnamelineage.dmp file
  taxaID2Lineage = {}
  for i in open(args.dir+'/fullnamelineage.dmp'):
      tmp = i.strip().upper().split('|')
      taxaID = tmp[0].strip()
      taxaName = tmp[1].lower().strip()
      lineage = tmp[2].lower().strip().replace('; ',';')+taxaName
      if 'eukaryota' in lineage:
          taxaID2Lineage[taxaID] = lineage
  return taxaID2Lineage

def used_classifiers():
	c = {}
	for i in open(args.d):
		tmp = i.strip().split('\t')
		MG,start,end = tmp[1],int(tmp[8]),int(tmp[9])
		region = str(int(20 * math.floor(min(start,end)/20)))
		MG_region = MG+'_'+region
		c[MG_region] = {}
	return c

def collect_thresholds(classifiers):
	for i in open(args.dir+'/strict_classifiers_filtered.txt'):
		tmp = i.strip().split('\t')
		region = tmp[0]
		if region in classifiers:
			thresholds = []
			for x in tmp[1:]:
				if x == 'na': thresholds.append('na')
				else: thresholds.append(float(x))
			classifiers[region] = thresholds
	return classifiers

def collect_classifier_metadata(taxaID2Lineage):
	metadata,mg_per_species,prot2len = {},{},{}
	
	for i in open(args.dir+'/marker_gene_metadata.txt'):
		tmp = i.strip().split('\t')
		prot,MG,prot_len,species,genus,family,lineage = tmp[0],tmp[1],float(tmp[2]),tmp[3],tmp[4],tmp[5],tmp[6]
		metadata[prot] = [MG,species,genus,family]
		species_lineage = taxaID2Lineage[species]
		prot2len[prot] = prot_len
		if species_lineage not in mg_per_species: mg_per_species[species_lineage] = set([MG])
		else: mg_per_species[species_lineage].update([MG])
	return metadata,mg_per_species,prot2len

def LCA(lineages): # function for finding the lowest common ancestor of a set of lineages
    lineages = [x.split(';') for x in lineages]
    new_lineage = [i[0].strip() for i in zip(*lineages) if len(set(i)) == 1 if '' not in i]
    return new_lineage

parser = argparse.ArgumentParser()
parser.add_argument("-d", help="diamond output file")
parser.add_argument("-dir", help="path to data directory of taxaTarget")
parser.add_argument("-p", help="proportion of range(smax - 0) to add to thresholds for species, genus and family",default=0.05,type=float)
args = parser.parse_args()

# extract number of reads in sample and mean read length
number_reads,read_len = read_file_info()
print('Extracted number of reads and mean read length for sample')

taxaID2Lineage = taxa_mappings()
print('collected taxonomy dictionary')

# maps number of marker genes per family, genus, and species per taxaID
MGs_per_taxa = {taxaID2Lineage[i.strip().split('\t')[0]]:i.strip().split('\t')[2:] for i in open(args.dir+'/MGs_per_species_genus_family.txt')}

classifiers = used_classifiers()
print('diamond file first pass')

classifiers = collect_thresholds(classifiers)
print('collected relevant classifiers')

metadata,MG_per_species,prot2len = collect_classifier_metadata(taxaID2Lineage)
print('collected relevant metadata for marker genes')

read_classifications = {}
multimapped = Counter([i.strip().split('\t')[0] for i in open(args.d)])
genus_species = {}

with open('classified_reads.txt','w') as out:
	for i in open(args.d):
		tmp = i.strip().split('\t')
		read,MG,pident,aln_len,start,end,bitscore = tmp[0],tmp[1],float(tmp[2]),float(tmp[3]),int(tmp[8]),int(tmp[9]),float(tmp[11])
		mean_bitscore = bitscore/aln_len
		if aln_len >= 30:
			region = str(int(20 * math.floor(min(start,end)/20)))
			MG_region = MG+'_'+region
			busco,species,genus,family = metadata[MG]
			if (multimapped[read] > 1) and (pident == 100):
				if read not in read_classifications:
					read_classifications[read] = [[taxaID2Lineage[species]],set([busco]),mean_bitscore,[prot2len[MG]]]
				else:
					read_classifications[read][0].append(taxaID2Lineage[species])
					read_classifications[read][1].update([busco])
					read_classifications[read][3].append(prot2len[MG])
			elif classifiers[MG_region] != {}:
				nmax,smin,smax,gmin,gmax,fmin,fmax = classifiers[MG_region]
				padding = args.p*smax
				lineage = ''
				if mean_bitscore >= smin + padding:
					lineage = taxaID2Lineage[species]
					out.write(read+'\t'+'species'+'\t'+busco+'\t'+MG_region+'\t'+str(mean_bitscore)+'\t'+lineage+'\n')
					if genus not in genus_species:
						genus_species[genus] = set([lineage])
					else:
						genus_species[genus].update([lineage])
				elif gmin != 'na':
					if mean_bitscore >= gmin + padding:
						lineage = taxaID2Lineage[genus]
						out.write(read+'\t'+'genus'+'\t'+busco+'\t'+MG_region+'\t'+str(mean_bitscore)+'\t'+lineage+'\n')
				elif fmin != 'na':
					if mean_bitscore >= fmin + padding:
						lineage = taxaID2Lineage[family]
						out.write(read+'\t'+'family'+'\t'+busco+'\t'+MG_region+'\t'+str(mean_bitscore)+'\t'+lineage+'\n')
				if lineage != '':
					if read not in read_classifications:
						read_classifications[read] = [[lineage],set([busco]),mean_bitscore,[prot2len[MG]]]
					else:
						read_classifications[read][0].append(lineage)
						read_classifications[read][1].update([busco])
						read_classifications[read][3].append(prot2len[MG])

mg_names = [i.strip() for i in open(args.dir+'/255_MGs.txt')] + ['Total_reads']
sample_taxa = {'names':mg_names}
final_read_classifications = {}

for read,v in read_classifications.items():
    lineages,buscos,mean_bitscore,average_prot_len = v[0],v[1],v[2],sum(v[3])/float(len(v[3]))
    lca = LCA(lineages)
    line = ';'.join(lca)
    if line not in sample_taxa:
        sample_taxa[line] = np.zeros(len(mg_names))
    if line not in final_read_classifications:
        final_read_classifications[line] = [read]
    else:
        final_read_classifications[line].append(read)
    for i in buscos:
        sample_taxa[line][mg_names.index(i)] += 1 # read_len/average_prot_len
        sample_taxa[line][-1] += 1 #read_len/average_prot_len

del_keys = []

for k,v in sample_taxa.items():
    if k == 'names': continue
    lineage_MGS = MGs_per_taxa.get(k,0)
    if lineage_MGS == 0: continue
    tmp = []
    non_zero = 0
    for i in lineage_MGS:
        rd_cnt = sample_taxa[line][mg_names.index(i)]
        tmp.append(rd_cnt)
        if rd_cnt > 0: non_zero += 1
    total_rd_cnt = sum(tmp)
    expected_number_MGs = 0.96270 * total_rd_cnt + -0.00128 * total_rd_cnt**2 + 1.89510
    if non_zero < 20:
        if non_zero/float(len(tmp)) < 0.25:
            if non_zero < expected_number_MGs*0.5:
                del_keys.append(k)
    # if np.std(tmp)/np.mean(tmp) >= 2.7: del_keys.append(k)

for i in del_keys:
	del sample_taxa[i]


'''
experimental was never used

    if lineage_total_MGS == 0: continue
    if median_depth_coverage >= 1: expected_MGs_mapped = lineage_total_MGS
    else: expected_MGs_mapped = median_depth_coverage * lineage_total_MGS
    expected_median_coverage = sample_taxa[line][-1] * read_len / expected_MGs_mapped # total_reads_mapped * read_len / expected_MGs_mapped
    print(line,sample_taxa[line][-1],median_depth_coverage,lineage_total_MGS,expected_MGs_mapped,expected_median_coverage,median_depth_coverage/expected_median_coverage)
    break
'''

df = pd.DataFrame(sample_taxa)
t = df.transpose()

t.to_csv('marker_gene_read_counts_per_taxa.txt',sep='\t')

print('marker gene counts collected')

for k,v in genus_species.items():
	s = {}
	for line in v:
		if line in sample_taxa:
			s[line] = sample_taxa[line]
	
	try:
		primary,primary_counts = max(s.items(),key=lambda x: x[1][-1])
		
		for z,w in s.items():
			if z != primary:
				shared,different = [0.1,0.1],[0.0,0.0] # [number_genes,number_reads]; shared are set to 0.1 to avoid division by zero errors
				for j in range(214):
					if (mg_names[j] in MG_per_species[primary]) and (mg_names[j] in MG_per_species[z]):
						shared[0] += 1
						shared[1] += w[j]
					elif mg_names[j] in MG_per_species[z]:
						different[0] += 1
						different[1] += w[j]
				if (different[1]/shared[1]) > ((different[0]/shared[0])+0.1*different[0]/shared[0]): # a false positive would be expected to have a lower number of reads mapping to shared genes than to not-shared genes; here within an epsilon of 10% is allowed
					del sample_taxa[z]
	except: continue

aggregate_taxa = {}

for k,v in sample_taxa.items():
    if k == 'names': continue
    counts = v[-1] # number of mapped reads
    buscos = {mg_names[h] for h,i in enumerate(v) if i > 0} # number of mapped marker genes
    if counts < 4:
        if k in final_read_classifications: del final_read_classifications[k]
        continue
    if len(buscos) < 2:
        if k in final_read_classifications: del final_read_classifications[k]
        continue
    '''if len(buscos) < 6:
        if len(buscos) < expected_number_MGs*0.5: 
            if k in final_read_classifications: del final_read_classifications[k]
            continue
    if len(buscos) < 6:
        if counts > 1.7**len(buscos):
            if k in final_read_classifications: del final_read_classifications[k]
            continue'''
    line = ''
    for j in k.split(';'):
        line += j+';'
        if line not in aggregate_taxa:
            aggregate_taxa[line] = [counts,buscos]
        else:
            aggregate_taxa[line][0] += counts
            aggregate_taxa[line][1].update(buscos)

with open('final_read_classifications.txt','w') as out:
	for k,v in final_read_classifications.items():
		for i in v:
			out.write(i+'\t'+k+'\n')

print('aggregated results. writing out to taxonomic report')

with open('Taxonomic_report.txt','w') as out:
	for k,v in sorted(aggregate_taxa.items()):
		counts,buscos = v[0],v[1]
		out.write(k+'\t'+str(counts)+'\t'+str(len(buscos))+'\n')





