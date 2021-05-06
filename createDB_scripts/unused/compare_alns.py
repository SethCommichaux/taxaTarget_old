import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-d", help="diamond alignment file")
parser.add_argument("-f", help="protist functions file")
parser.add_argument("-i", help="msa global coordinates file")
parser.add_argument("-names", help="NCBI names file")
parser.add_argument("-nodes", help="NCBI nodes file")
parser.add_argument("-o", help="output file")
args = parser.parse_args()


def create_dicts(names,nodes):
	uniref2lineage = {}

	for h,i in enumerate(open(args.f)):
		if h == 0: continue
		tmp = i.strip().split('\t')
		uniref = tmp[0]
		lineage = map(lambda x:x.strip(),tmp[3].split(';'))
		uniref2lineage[uniref] = {}
		for j in lineage:
			rank = nodes[names[j]]
			if rank in ['family','genus','species']:
				uniref2lineage[uniref][rank] = names[j]

	print 'Built dictionaries...'
	return uniref2lineage

def global_coordinates():
	indices = {}
	for i in open(args.i):
		tmp = i.strip().split('\t')
		uniref = tmp[0]
		f_name = tmp[1]
		f_gaps = map(int,tmp[2].split(','))
		f_gap_lens = map(int,tmp[3].split(','))
		g_name = tmp[4]
		g_gaps = map(int,tmp[5].split(','))
		g_gap_lens = map(int,tmp[6].split(','))
		s_name = tmp[7]
		s_gaps = map(int,tmp[8].split(','))
		s_gap_lens = map(int,tmp[9].split(','))
		indices[uniref] = {'F':[f_gaps,f_gap_lens,f_name],'G':[g_gaps,g_gap_lens,g_name],'S':[s_gaps,s_gap_lens,s_name]}
	print "Extracted global MSA coordinates!!!"
	return indices


def get_taxa_rank(x,y,r):
	x_rank = uniref2lineage.get(x,0)
	y_rank = uniref2lineage.get(y,0)
	if y_rank == 0:
		return "NA"
	if x_rank == 0: return 'N'
	rank = ''
	for i in ['species','genus','family']:
		if x_rank.get(i,0) == y_rank.get(i,1):
			rank = i
			return r[rank]
	if rank == '': return 'N'


def find_global_coordinates(gaps,start):
	gap_starts = gaps[0]
	gap_lens = gaps[1]
	taxaName = gaps[2]
	if gap_starts == ['']: return start,taxaName
	else:
		for i in range(len(gap_starts)):
			if start >= gap_starts[i]: start += gap_lens[i]
	return start,taxaName


names = {i.strip().split('\t')[1].upper():i.strip().split('\t')[0] for i in open(args.names)}
nodes = {i.strip().split('\t')[0]:i.strip().split('\t')[1] for i in open(args.nodes)}
indices = global_coordinates()
uniref2lineage = create_dicts(names,nodes)

r = {'species':'S','genus':'G','family':'F'}

# format of diamonod output
# neg     UniRef100_A0A1Y2HDR0    47.8    67      35      0       2       68      52      118     4.3e-14 77.8
# UniRef100_E5A3J1_10     UniRef100_E5A3J1	100.0   70      0       0       1       70      11      80      1.5e-32 139.0

with open(args.o,'w') as out:
	for i in open(args.d):
		tmp = i.strip().split('\t')
		query = tmp[0]
		subject = tmp[1]
		sstart = int(tmp[8])
		aln_len = float(tmp[3])
		bitscore = float(tmp[11])
		slope = bitscore/aln_len
		if query == 'neg':
			relation = 'N'
		else:
			query_id = '_'.join(query.split(' ')[0].split('_')[:2])
			relation = get_taxa_rank(query_id,subject,r)
		if subject in indices:
			if relation == 'S':
				global_msa_start,taxaName = find_global_coordinates(indices[subject]['S'],sstart)
				if taxaName != '0':
					out.write(taxaName+'\t'+str(global_msa_start)+'\t'+str(slope)+'\t'+relation+'\n')
			elif relation == 'G':
				global_msa_start,taxaName = find_global_coordinates(indices[subject]['G'],sstart)
				if taxaName != '0':
					out.write(taxaName+'\t'+str(global_msa_start)+'\t'+str(slope)+'\t'+relation+'\n')
				global_msa_start,taxaName = find_global_coordinates(indices[subject]['S'],sstart)
				if taxaName != '0':
					out.write(taxaName+'\t'+str(global_msa_start)+'\t'+str(slope)+'\t'+'NS'+'\n')
			elif relation == 'F':
				global_msa_start,taxaName = find_global_coordinates(indices[subject]['F'],sstart)
				if taxaName != '0':
					out.write(taxaName+'\t'+str(global_msa_start)+'\t'+str(slope)+'\t'+relation+'\n')
				global_msa_start,taxaName = find_global_coordinates(indices[subject]['S'],sstart)
				if taxaName != '0':
					out.write(taxaName+'\t'+str(global_msa_start)+'\t'+str(slope)+'\t'+'NS'+'\n')
				global_msa_start,taxaName = find_global_coordinates(indices[subject]['G'],sstart)
				if taxaName != '0':
					out.write(taxaName+'\t'+str(global_msa_start)+'\t'+str(slope)+'\t'+'NG'+'\n')
			elif relation == 'N':
				for j in 'FGS':
					global_msa_start,taxaName = find_global_coordinates(indices[subject][j],sstart)
					if taxaName != '0':
						out.write(taxaName+'\t'+str(global_msa_start)+'\t'+str(slope)+'\t'+relation+j+'\n')




