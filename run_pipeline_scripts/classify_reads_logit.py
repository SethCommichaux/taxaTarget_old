import argparse
import ast
import math

def logistic_prob(slope,intercept,x):
	if 'na' in [slope,intercept]: return 0
	return 1/(1+(math.exp(-1*(intercept+slope*x))))

parser = argparse.ArgumentParser()
parser.add_argument("-d", help="diamond output file")
parser.add_argument("-m", help="marker gene mapping file")
parser.add_argument("-c", help="marker gene thresholds file")
parser.add_argument("-t", help="probability threshold",default=0.8,type=float)
args = parser.parse_args()

classifiers = {}

for i in open(args.d):
	tmp = i.strip().split('\t')
	MG,start,end = tmp[1],int(tmp[8]),int(tmp[9])
	region = str(int(20 * math.floor(min(start,end)/20)))
	MG_region = MG+'_'+region
	classifiers[MG_region] = {}

for i in open(args.c):
	tmp = i.strip().split('\t')
	region = tmp[0]
	if region in classifiers:
		thresholds = ast.literal_eval(tmp[1])
		classifiers[region] = thresholds

for i in open(args.d):
	tmp = i.strip().split('\t')
	read,MG,aln_len,start,end,bitscore = tmp[0],tmp[1],float(tmp[3]),int(tmp[8]),int(tmp[9]),float(tmp[11])
	region = str(int(20 * math.floor(min(start,end)/20)))
	MG_region = MG+'_'+region
	if classifiers[MG_region] != {}:
		mean_bitscore = bitscore/aln_len
		sslope,sint = classifiers[MG_region].get('s',['na','na'])
		sprob = logistic_prob(sslope,sint,mean_bitscore)
		gslope,gint = classifiers[MG_region].get('g',['na','na'])
		gprob = logistic_prob(gslope,gint,mean_bitscore)
		fslope,fint = classifiers[MG_region].get('f',['na','na'])
		fprob = logistic_prob(fslope,fint,mean_bitscore)
		probs = [fprob,gprob,sprob]
		ranks = 'fgs'
		rank = 'n'
		for x,prob in enumerate(probs):
			if prob >= args.t:
				rank = ranks[x]
			else:
				break
		if rank != 'n':
			print(read,rank,probs,classifiers[MG_region],mean_bitscore)





'''


		elif mean_bitscore >= classifiers[MG_region]['s']: print(read,'s',mean_bitscore,classifiers[MG_region]['s'])
		elif mean_bitscore >= classifiers[MG_region].get('g',100): print(read,'g',mean_bitscore,classifiers[MG_region].get('g',100))
		elif mean_bitscore >= classifiers[MG_region].get('f',100): print(read,'f',mean_bitscore,classifiers[MG_region].get('f',100))
'''
