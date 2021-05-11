import sys
import numpy as np
from sklearn.linear_model import LogisticRegression
import argparse

def train_logistic_regression(x):
	x = np.array(x)
	labels = [1,1,0,0]
	X = x[:, np.newaxis]
	clf = LogisticRegression(C=1e5)
	clf.fit(X, labels)
	return clf.coef_[0][0],clf.intercept_[0]

parser = argparse.ArgumentParser()
parser.add_argument("-s", help="strict classifier file")
parser.add_argument("-o", help="desired output file")
args = parser.parse_args()

f = args.s   # 'strict_classifiers.txt'

with open(args.o,'w') as out:
	for h,i in enumerate(open(f)):
		if h%10000 == 0: print(h)
		line = i.strip().split('\t')
		n,smin,smax,gmin,gmax,fmin,fmax = float(line[1]),float(line[2]),float(line[3]),float(line[4]),float(line[5]),float(line[6]),float(line[7])
		if smin == 100: continue # skip regions with no species-specific read mappings
		if gmin == 100: gmin = 0
		if fmin == 100: fmin = 0
		if smax > max([n,gmin,gmax,fmin,fmax]):
			tmp = {}
			SMIN = max(smin,gmax,gmin,fmax,fmin,n)
			x = [smax,SMIN,max(n,SMIN/3.),0]
			slope,intercept = train_logistic_regression(x)		
			tmp['s'] = [slope,intercept]
			if fmax > n:
				FMIN = max(fmin,n)
				x = [fmax,FMIN,max(n,FMIN/3.),0]
				slope,intercept = train_logistic_regression(x)
				tmp['f'] = [slope,intercept]
			if gmax > max(n,fmax):
				GMIN = max(gmin,n,fmax)
				x = [gmax,GMIN,max(n,GMIN/3.),0]
				slope,intercept = train_logistic_regression(x)
				tmp['g'] = [slope,intercept]
			out.write(line[0]+'\t'+str(tmp)+'\n')	
