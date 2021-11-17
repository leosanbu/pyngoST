import os
import glob
import subprocess
import pickle
import ahocorasick
import pandas as pd
from Bio import SeqIO, AlignIO
from pyfaidx import Fasta
from Bio.Align.Applications import MuscleCommandline

## Classes ##

class Allele:
	def __init__(self, gene, allele):
		self.gene = gene
		self.allele = allele

## Functions ##

def getInput(flist, fread):
	if fread:
		filelist = []
		with open(fread, 'r') as paths:
			for genomepath in paths:
				filelist.append(genomepath.rstrip())
	else:
		filelist = flist
	return filelist

def readAlleles(path):
	allelesDB = {}
	for files in glob.glob(path+'/'+'*.fas'):
		gene = files.replace('.fas', '')
		gene = gene.split('/').pop()
		with open(files, 'r') as fasta:
			for record in SeqIO.parse(fasta, 'fasta'):
				name = record.name
				name = name.split('_')[1]
				if 'penA' not in record.name:
					name = name.split('.')[0]
				seq = record.seq
				rc = revComp(seq)
				allelesDB[str(seq)] = Allele(gene, name)
				allelesDB[str(rc)] = Allele(gene, name)
	return allelesDB

def readProfiles(path):
	pathsfile = glob.glob(path+'/'+'*profiles.tab')
	paths = []
	for p in pathsfile:
		tmp = p.split('/').pop().replace('_profiles.tab', '')
		paths.append(tmp)
	profilesdic = {}
	for n, p in enumerate(pathsfile):
		profilesDB = {}
		with open(p, 'r') as tab:
			for line in tab:
				linesplit = line.rstrip().split('\t')
				st = linesplit[0]
				profile = '\t'.join(linesplit[1:8])
				profilesDB[profile] = st
		profilesdic[paths[n]] = profilesDB
	return profilesdic

def readNGSTARCCs(file):
	ngstarCCs = {}
	with open(file, 'r') as tab:
		for line in tab:
			if not line.startswith('Sequence Type'):
				linesplit = line.rstrip().split(',')
				ngstarCCs[linesplit[0]] = linesplit[1]
	return ngstarCCs

def revComp(seq):
	basesDic = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
	seqlist = list(seq)
	revseq = list(reversed(seqlist))
	comp = [basesDic[x] for x in revseq]
	return ''.join(comp)

def reportProfile(order, results):
	profile = ''
	for i in order:
		# check if multiple alleles (i.e. 23S)
		check_copies = results[i]
		if len(check_copies)>1:
			unique_copies = '_'.join(list(pd.unique(results[i])))
			if profile is '':
				profile = unique_copies
			else:
				profile += '\t'+unique_copies
		else:
			if profile is '':
				profile = check_copies[0]
			else:
				profile += '\t'+check_copies[0]
	return profile

def blastNewAlleles(query, subject, path):
	# run makeblastdb on subject (genome)
	subject_name = subject.split('/').pop()
	makeblastdb_cmd = ['makeblastdb', '-in', subject, '-dbtype', 'nucl', '-out', 'tmp/'+subject_name]
	subprocess.call(makeblastdb_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
	# run blastn
	blast_cmd = ['blastn', '-out', 'tmp/'+subject_name+'.blastn', '-outfmt', '6', '-query', path+'/'+query, '-db', 'tmp/'+subject_name, '-evalue', '0.001']
	subprocess.call(blast_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
	ret = ['-']
	# read tabular output with pandas if file exists
	if os.path.getsize('tmp/'+subject_name+'.blastn') > 0:
		bresult = pd.read_csv('tmp/'+subject_name+'.blastn', sep='\t', header=None)
		# get value of max bitscore
		max_bitscore = bresult[11].max()
		# get alleles with same bitscore
		subresult = bresult.loc[bresult[11] == max_bitscore]
		sublist = subresult[[0]][0].tolist()
		# get coordinates of queried alleles
		coords = subresult.iloc[0,8:10].tolist()
		# get contig name of hit
		contigloc = subresult.iloc[0,1]
		# get number from closest allele names
		clean_alleles = []
		for x in sublist:
			newx = x.split('_')[1]
			if 'penA' not in x:
				newx = newx.split('.')[0]
			clean_alleles.append(newx)
		ret = ['|'.join(clean_alleles[0:3])+'*', coords, contigloc]
	return ret

def printNewAlleleSeqs(gene, coords, contigloc, fasta, allout, path):
	strand = 1
	startcoord, endcoord = coords
	if startcoord > endcoord:
		strand = -1
		endcoord, startcoord = coords
	# read genome and extract locus
	genome = Fasta(fasta)
	locus = genome[contigloc][(startcoord-1):endcoord]
	if strand == -1:
		locus.seq = revComp(locus.seq)
	# print to file
	outfile = fasta.split('/').pop()+'.'+gene+'.fasta'
	if allout:
		with open(path+'/'+outfile, 'w') as out:
			out.write('>'+gene+'_'+locus.name+'_'+str(startcoord)+':'+str(endcoord)+'\n')
			out.write(locus.seq+'\n')
	return locus
	
def updatePKL(path, allelesDB):
	allelesAC = ahocorasick.Automaton()
	for idx,key in enumerate(allelesDB):
		allelesAC.add_word(key, (idx, key))
	allelesAC.make_automaton()
	pickle.dump((allelesDB, allelesAC), open(path+'/'+'ST_alleles_AC.pkl', 'wb'))

def AC_fast(seq, order, allelesDB, allelesAC):
	results = {}
	for i in order:
		results[i] = '-'
	for end_index, (insert_order, original_value) in allelesAC.iter(seq):
		if allelesDB[original_value].gene in order:
			if type(results[allelesDB[original_value].gene]) is not list:
				results[allelesDB[original_value].gene] = [allelesDB[original_value].allele]
			else:
				results[allelesDB[original_value].gene].append(allelesDB[original_value].allele)
	return results

def run_ngmaster(filelist, genogroups):
	if genogroups:
		filelist_str = 'ngmaster --printseq ngmaster_printseq.fas '+' '.join(filelist)
	else:
		filelist_str = 'ngmaster '+' '.join(filelist)
	ngmaster_output = subprocess.run(filelist_str, shell=True, capture_output=True).stdout.decode('utf-8')
	ngmaster_output = ngmaster_output.split('\n')[1:-1]
	ngmastDic = {}
	ngmastClusters = {}
	for i in ngmaster_output:
		tmp = i.split('\t')
		ngmastDic[tmp[0]] = '\t'.join(tmp[1:])
		if genogroups:
			if not tmp[1] in ngmastClusters:
				ngmastClusters[tmp[1]] = []
				#ngmastClusters[tmp[1]] = [tmp[0], tmp[2], tmp[3]]
			ngmastClusters[tmp[1]].append([tmp[0], tmp[2], tmp[3]])
	return [ngmastDic, ngmastClusters]

def split_loci_sequences(file, out_path):
	check_files = os.listdir(out_path)
	if 'POR.fas' in check_files:
		subprocess.run(['rm POR.fas'], shell=True)
	if 'TBPB.fas' in check_files:
		subprocess.run(['rm TBPB.fas'], shell=True)
	por_fas = open('POR.fas', 'a')
	tbpb_fas = open('TBPB.fas', 'a')
	with open(out_path+'/'+file, 'r') as fasta:
		for record in SeqIO.parse(fasta, 'fasta'):
			if 'POR' in record.description:
				SeqIO.write(record, por_fas, 'fasta')
			elif 'TBPB' in record.description:
				SeqIO.write(record, tbpb_fas, 'fasta')
	por_fas.close()
	tbpb_fas.close()

def align_sequences(infilename, outfilename): #use maxiters=2 for speed
	cline = MuscleCommandline("muscle", input=infilename, out=outfilename) #maxiters=maxiters
	try:
		stdout, stderr = cline()
	except OSError:
		print("Alignment failed:")
		print(cline)
		sys.exit()

def calculate_distances(aln):
	align = AlignIO.read(aln, "fasta")
	dist={}
	for recorda in align:
		seqa=str(recorda.seq).upper()
		dist[recorda.id]={}
		for recordb in align:
			seqb=str(recordb.seq).upper()
			count=0
			for x in range(len(seqa)):
				if seqa[x]!=seqb[x] and seqa[x] in ["A","C","G","T"] and seqb[x] in ["A","C","G","T"]:
					count+=1
			dist[recorda.id][recordb.id]=count
	return dist

def prepare_files_for_genogroups(out_path):
	split_loci_sequences('ngmaster_printseq.fas', out_path)
	align_sequences('POR.fas', 'POR.aln')
	align_sequences('TBPB.fas', 'TBPB.aln')

def calculateGenogroups(out_path, ngmastClusters):
	por_dist = calculate_distances('POR.aln')
	tbpb_dist = calculate_distances('TBPB.aln')
	# calculate clusters
	clusters=[]
	for st in ngmastClusters:
		if st is not '-':
			s=ngmastClusters[st][0]+[st]
			hit_clusters=set([])
			for x, cluster in enumerate(clusters):
				for c in cluster:
					hit=False
					if s[0] in por_dist and s[0] in tbpb_dist:
						if c[0] in por_dist[s[0]] and c[0] in tbpb_dist[s[0]]:
							if por_dist[s[0]][c[0]]<=5 and tbpb_dist[s[0]][c[0]]==0:
								hit=True
							elif tbpb_dist[s[0]][c[0]]<=4 and por_dist[s[0]][c[0]]==0:
								hit=True
							elif tbpb_dist[s[0]][c[0]]+por_dist[s[0]][c[0]]<=8:
								hit=True
							if hit:
								hit_clusters.add(x)
								continue
			hit_clusters=list(hit_clusters)
			if len(hit_clusters)==0:
				clusters.append([s])
			elif len(hit_clusters)==1:
				clusters[hit_clusters[0]].append(s)
			else:
				hit_clusters.sort()
				clusters[hit_clusters[0]].append(s)
				for c in hit_clusters[1:]:
					clusters[hit_clusters[0]]=clusters[hit_clusters[0]]+clusters[c]
				hit_clusters.reverse()
				for c in hit_clusters[:-1]:
					clusters.pop(c)
	# get genogroups
	genogroups={}
	output=open("STs_in_genogroups.tsv", "w")
	for cluster in clusters:
		g=set([])
		for c in cluster:
			g.add(int(c[-1]))
		sts=[]
		for c in g:
			sts.append([len(ngmastClusters[str(c)]), str(c)])
		sts.sort()
		sts.reverse()
		genogroup_name="G"+str(sts[0][1])
		stssorted = []
		for s in sts:
			genogroups[s[1]]=genogroup_name
			stssorted.append(s[1]+'('+str(s[0])+')')
		print(genogroup_name, '; '.join(map(str,stssorted)), file=output, sep='\t')
	output.close()
	genopergenome = {}
	for n in ngmastClusters:
		for sample in ngmastClusters[n]:
			if n is not '-':
				genopergenome[sample[0]] = genogroups[n]
			else:
				genopergenome[sample[0]] = '-'
	return genopergenome

def processFiles(fpath, order, MLSTorder, NGSTARorder, profilesDB, allelesDB, allelesAC, schemes, 
	ngstarCCs, ngmastDic, genopergenome, db_path, allout, out_path, outfile, outfilehandle):
	fname = fpath.split('/').pop()
	with open(fpath, 'r') as fasta:
		concat = ''
		for record in SeqIO.parse(fasta, 'fasta'):
			concat += record.seq
		#AC algorithm:
		resultsDB = AC_fast(str(concat), order, allelesDB, allelesAC)
		# run Blast to get closest alleles to missing loci
		save_new_alleles = {}
		for x in resultsDB:
			if resultsDB[x] == '-':
				query_file = x+'.fas'
				subject_file = fpath
				blastout = blastNewAlleles(query_file, subject_file, db_path)
				if len(blastout)>1:
					need_blast = 1
					closest_allele, coords, contigloc = blastout
					resultsDB[x] = closest_allele
					save_new_alleles[x] = printNewAlleleSeqs(x, coords, contigloc, fpath, allout, out_path)
	# get profile and assign ST
	defaultST = 'NEW'
	st_list = {}
	for s in schemes:
		if s=='MLST':
			prof = reportProfile(MLSTorder, resultsDB)
			if '-' in prof:
				st = '-'
			else:
				st = profilesDB['MLST'].get(prof, defaultST)
			st_list['MLST'] = st+'\t'+prof
		elif s=='NG-STAR':
			prof = reportProfile(NGSTARorder, resultsDB)
			if '-' in prof:
				st = '-'
			else:
				st = profilesDB['NGSTAR'].get(prof, defaultST)
			if ngstarCCs:
				if st=='NEW' or st=='-':
					st_list['NG-STAR'] = st+'\t'+prof+'\t'+'-'
				else:
					st_list['NG-STAR'] = st+'\t'+prof+'\tCC'+ngstarCCs[st]
			else:		
				st_list['NG-STAR'] = st+'\t'+prof
		elif s=='NG-MAST':
			ngmast = ngmastDic[fpath]
			if genopergenome:
				st_list['NG-MAST'] = ngmast+'\t'+genopergenome[fpath]
			else:
				st_list['NG-MAST'] = ngmast
	outline = ''
	for s in st_list:
		outline += '\t'+st_list[s]
	if outfile:	
		outfilehandle.write(fname+outline+'\n')
	else:
		print(fname+outline)








