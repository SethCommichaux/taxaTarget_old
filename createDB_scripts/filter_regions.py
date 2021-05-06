f = 'strict_classifiers.txt'

with open('strict_classifiers_filtered.txt','w') as out:
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
			tmp['s'] = SMIN
			tmp['n'] = n
			if fmax > n:
				FMIN = max(fmin,n)
				tmp['f'] = FMIN
			if gmax > max(n,fmax):
				GMIN = max(gmin,n,fmax)
				tmp['g'] = GMIN
			out.write(line[0]+'\t'+str(tmp)+'\n')	
