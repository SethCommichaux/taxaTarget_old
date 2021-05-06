import argparse
import ast
import math

parser = argparse.ArgumentParser()
parser.add_argument("-d", help="diamond output file")
parser.add_argument("-m", help="marker gene mapping file")
parser.add_argument("-c", help="marker gene thresholds file")
args = parser.parse_args()

classifiers_to_use = {}

for i in open(args.d):
        tmp = i.strip().split('\t')
        MG,start,end = tmp[1],int(tmp[8]),int(tmp[9])
        region = str(int(20 * math.floor(min(start,end)/20)))
        MG_region = MG+'_'+region
	classifiers_to_use[MG_region] = 0

classifiers = {}

for i in open(args.c):
	tmp = i.strip().split('\t')
	region = tmp[0]
	if region in classifiers_to_use:
		thresholds = ast.literal_eval(tmp[1])
		classifiers[region] = thresholds

for i in open(args.d):
	tmp = i.strip().split('\t')
	read,MG,aln_len,start,end,bitscore = tmp[0],tmp[1],float(tmp[3]),int(tmp[8]),int(tmp[9]),float(tmp[11])
	region = str(int(20 * math.floor(min(start,end)/20)))
	MG_region = MG+'_'+region
	if MG_region in classifiers:
		mean_bitscore = bitscore/aln_len
		if mean_bitscore >= classifiers[MG_region]['s']: print read,'s',mean_bitscore,classifiers[MG_region]['s']
		elif mean_bitscore >= classifiers[MG_region].get('g',100): print read,'g',mean_bitscore,classifiers[MG_region].get('g',100)
		elif mean_bitscore >= classifiers[MG_region].get('f',100): print read,'f',mean_bitscore,classifiers[MG_region].get('f',100)
