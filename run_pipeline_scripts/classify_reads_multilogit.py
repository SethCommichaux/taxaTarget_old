import pickle
import argparse
import ast
import math

parser = argparse.ArgumentParser()
parser.add_argument("-d", help="diamond output file")
parser.add_argument("-m", help="marker gene mapping file")
parser.add_argument("-t", help="probability threshold",default=0.8,type=float)
parser.add_argument("-p", help="pickles file")
args = parser.parse_args()

with open(args.p,'rb') as f:
	classifiers = pickle.load(f)

print(list(classifiers.items())[:10])

print('pickle loaded!')

for i in open(args.d):
	tmp = i.strip().split('\t')
	read,MG,aln_len,start,end,bitscore = tmp[0],tmp[1],float(tmp[3]),int(tmp[8]),int(tmp[9]),float(tmp[11])
	region = str(int(20 * math.floor(min(start,end)/20)))
	MG_region = MG+'_'+region
	if MG_region not in classifiers: continue
	mean_bitscore = bitscore/aln_len
	probs = classifiers[MG_region].predict_proba([[mean_bitscore]])[0]
	classes = classifiers[MG_region].classes_
	m = max(probs)
	next_best = 0
	rank = []
	for h,j in enumerate(probs):
		if j == m:
			rank.append(classes[h])
		else:
			next_best = max(next_best,j)
	if 'n' not in rank:
		if m - next_best >= 1./len(classes)/2.:
			print(read,rank,next_best,probs,classes,mean_bitscore)




