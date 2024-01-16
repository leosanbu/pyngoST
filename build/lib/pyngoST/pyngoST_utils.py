import os
import re
import sys
import glob
import subprocess
import openpyxl
import pickle
import requests
import ahocorasick
import pandas as pd
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pyfaidx import Fasta
from Bio.Align.Applications import MuscleCommandline

## Classes ##

class Allele:
	def __init__(self, gene, allele, revcomp):
		self.gene = gene
		self.allele = allele
		self.revcomp = revcomp

## Functions ##

def get_args(arg):
	# Retrieve and return the command-line arguments as a dictionary
	args = {
		'db_path': arg.path if arg.path else False,
		'db_name': arg.db_name if arg.db_name else 'allelesDB',
		'out_path': arg.out_path if arg.out_path else os.getcwd(),
		'schemes': arg.schemes.split(',') if arg.schemes else ['NG-STAR'],
		'ngstarccs': arg.ngstarccs if arg.ngstarccs else False,
		'mosaic_pena': arg.mosaic_pena if arg.mosaic_pena else False,
		'genogroups': True if arg.genogroups else False,
		'blast_new_alleles': True if arg.blast_new_alleles else False,
		'allout': True if arg.alleles_out else False,
		'update': True if arg.update else False,
		'download': True if arg.download_db else False,
		'ccsfile': arg.ngstarccsfile if arg.ngstarccsfile else None,
		'outfile': arg.out_filename if arg.out_filename else None,
		'num_threads': arg.num_threads if arg.num_threads else 1,
		'only_assignccs': int(arg.only_assignccs) if arg.only_assignccs else None,
		'only_assignsts': arg.only_assignsts if arg.only_assignsts else None
    }
	if args['ngstarccs'] is True and 'NG-STAR' not in args['schemes']:
		print('## NG-STAR CCs requested but not NG-STAR. Please include NG-STAR in -s\n')
		sys.exit()
	if args['genogroups'] is True and 'NG-MAST' not in args['schemes']:
		print('## NG-MAST genogroups requested but not NG-MAST. Please include NG-MAST in -s\n')
		sys.exit()
	return args

def get_input(flist, fread):
	if fread:
		filelist = []
		with open(fread, 'r') as paths:
			for genomepath in paths:
				filelist.append(genomepath.rstrip())
	else:
		filelist = flist
	return filelist

def download_db(db_path, db_name, ccsfile):
	if db_name and db_name not in os.listdir(os.getcwd()):
		print('## Creating updated database in', db_name)
		subprocess.run(['mkdir '+db_name], shell=True)
		db_path = db_name
	else:
		print('## Folder', db_name, 'already exists. Exiting...\n')
		sys.exit()
	download_updated_dbs(db_path, ccsfile)
	allelesDB = read_alleles(db_path)
	make_ACautomaton(db_path, allelesDB)
	sys.exit()

def update_db(db_path, ccsfile):
	if not db_path:
		print('## Please, specify the path to an existing database with -p. Exiting...\n')
		sys.exit()
	else:
		print('## Updating database in', db_path)
	if ccsfile:
		print('## Integrating NG-STAR CCs from', ccsfile)
		integrate_ngstar_ccs(db_path, 'NGSTAR_profiles.tab', ccsfile)
	allelesDB = read_alleles(db_path)
	make_ACautomaton(db_path, allelesDB)
	sys.exit()

def load_db(db_path):
	# Load pickled database - dictionary and automaton    
	print('## Loading databases...')
	allelesDB, allelesAC = pickle.load(open(db_path+'/'+'ST_alleles_AC.pkl', 'rb'))
	return allelesDB, allelesAC

def download_updated_dbs(db_path, ccsfile):
	if ccsfile:
		ngstarccsfile = read_ngstar_ccs(ccsfile)
	# MLST and NG-MAST v2 from PubMLST
	print('## Downloading MLST and NG-MAST v2 alleles and Sequence Type profiles from https://pubmlst.org')
	loci = ['abcZ', 'adk', 'aroE', 'fumC', 'gdh', 'pdhC', 'pgm', 'NG-MAST_porB', 'NG-MAST_tbpB']
	pubmlstURLloci = 'https://pubmlst.org/bigsdb?db=pubmlst_neisseria_seqdef&page=downloadAlleles&locus='
	pubmlstMLSTs = 'https://pubmlst.org/bigsdb?db=pubmlst_neisseria_seqdef&page=downloadProfiles&scheme_id=1'
	pubmlstNGMASTs = 'https://pubmlst.org/bigsdb?db=pubmlst_neisseria_seqdef&page=downloadProfiles&scheme_id=71'
	for g in loci:
		r = requests.get(pubmlstURLloci+g, allow_redirects=True)
		open(db_path+'/'+g+'.fas', 'wb').write(r.content)
	m = requests.get(pubmlstMLSTs, allow_redirects=True)
	open(db_path+'/MLST_profiles.tab', 'wb').write(m.content)
	n = requests.get(pubmlstNGMASTs, allow_redirects=True)
	open(db_path+'/NGMAST_profiles.tab', 'wb').write(n.content)
	# Formatting NG-MAST files
	por_fas = open(db_path+'/POR.fas', 'w+')
	with open(db_path+'/NG-MAST_porB.fas', 'r') as fasta:
		for record in SeqIO.parse(fasta, 'fasta'):
			record.id = record.id.replace('NG-MAST_porB_', 'POR_')
			record.description = ''
			SeqIO.write(record, por_fas, 'fasta')
	if 'NG-MAST_porB.fas' in os.listdir(db_path):
		subprocess.call(['rm', db_path+'/NG-MAST_porB.fas'])
	tbpb_fas = open(db_path+'/TBPB.fas', 'w+')
	with open(db_path+'/NG-MAST_tbpB.fas', 'r') as fasta:
		for record in SeqIO.parse(fasta, 'fasta'):
			record.id = record.id.replace('NG-MAST_tbpB_', 'TBPB_')
			record.description = ''
			SeqIO.write(record, tbpb_fas, 'fasta')
	if 'NG-MAST_tbpB.fas' in os.listdir(db_path):
		subprocess.call(['rm', db_path+'/NG-MAST_tbpB.fas'])
	# NG-STAR from ngstar.canada.ca
	print('## Downloading NG-STAR alleles and Sequence Type profiles from https://ngstar.canada.ca')
	if ccsfile:
		print('## Integrating NG-STAR CCs from', ccsfile)
	locistar = ['penA', 'mtrR', 'porB', 'ponA', 'gyrA', 'parC', '23S']
	ngstarURLloci = 'https://ngstar.canada.ca/alleles/download?lang=en&loci_name='
	ngstarNGSTARs = 'https://ngstar.canada.ca/sequence_types/download?lang=en'
	ngstarMOSAIC = 'https://ngstar.canada.ca/alleles/download_metadata?lang=en&loci_name=penA'
	# Download allele fasta files
	for g in locistar:
		r = requests.get(ngstarURLloci+g, allow_redirects=True, verify=False)
		open(db_path+'/'+g+'.fas', 'wb').write(r.content)
	# Download and process profiles
	s = requests.get(ngstarNGSTARs, allow_redirects=True, verify=False)
	open(db_path+'/NGSTAR_profiles.xlsx', 'wb').write(s.content)
	df = pd.read_excel(db_path+'/NGSTAR_profiles.xlsx', sheet_name='Sheet 1')
	with open(db_path+'/NGSTAR_profiles.tmp.tab', 'w') as outfile:
		df.to_string(outfile, index=None)
	outfile = open(db_path+'/NGSTAR_profiles.tab', 'w+')
	with open(db_path+'/NGSTAR_profiles.tmp.tab', 'r') as profiles:
		for line in profiles:
			if 'Sequence Type':
				line=line.replace('Sequence Type', 'ST')
			linesplit = line.strip().split()
			if ccsfile:
				if linesplit[0] in ngstarccsfile:
					linesplit.append(ngstarccsfile[linesplit[0]])
				else:
					if linesplit[0] == 'ST':
						linesplit.append('NG-STAR_CC')
					else:
						linesplit.append('-')
			outfile.write('\t'.join(linesplit)+'\n')
	if 'NGSTAR_profiles.tmp.tab' in os.listdir(db_path):
		subprocess.call(['rm', db_path+'/NGSTAR_profiles.tmp.tab'])
	# Download penA alleles metadata to extract mosaic information
	r = requests.get(ngstarMOSAIC, allow_redirects=True, verify=False)
	open(db_path+'/penA_alleles_metadata.xlsx', 'wb').write(r.content)
	df2 = pd.read_excel(db_path+'/penA_alleles_metadata.xlsx', sheet_name='Sheet 1')
	with open(db_path+'/penA_metadata.tmp.tab', 'w') as outfile2:
		df2.to_string(outfile2, index=None)
	outmeta = open(db_path+'/penA_mosaics.tab', 'w+')
	outmeta.write('penA_allele'+'\t'+'mosaic_type'+'\n')
	with open(db_path+'/penA_metadata.tmp.tab', 'r') as metadata:
		for line in metadata:
			linesplit=re.split(' |;|,|:', line)
			linesplit2=[x for x in linesplit if x!='']
			if linesplit2[0] != 'Allele':
				found = 0
				for x in linesplit2:
					if 'Mosaic' in x:
						mosaicline = linesplit2[0]+'\t'+x
						outmeta.write(mosaicline+'\n')
						found = 1
					elif 'mosaic' in x:
						mosaicline = linesplit2[0]+'\t'+x
						outmeta.write(mosaicline+'\n')
						found = 1
				if found == 0:
					mosaicline = linesplit2[0]+'\t-'
					outmeta.write(mosaicline+'\n')
	if 'penA_metadata.tmp.tab' in os.listdir(db_path):
		subprocess.call(['rm', db_path+'/penA_metadata.tmp.tab'])

def read_alleles(path):
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
				rc = revcomp(seq)
				allelesDB[str(seq)] = Allele(gene, name, 0)
				allelesDB[str(rc)] = Allele(gene, name, 1)
	return allelesDB

def read_profiles(path, ngstarccs, mosaic_pena):
	pathsfile = glob.glob(path+'/'+'*profiles.tab')
	paths = []
	for p in pathsfile:
		tmp = p.split('/').pop().replace('_profiles.tab', '')
		paths.append(tmp)
	profilesdic = {}
	ngstardic = {}
	for n, p in enumerate(pathsfile):
		profilesDB = {}
		with open(p, 'r') as tab:
			for line in tab:
				linesplit = line.rstrip().split('\t')
				st = linesplit[0]
				profile = '\t'.join(linesplit[1:8])
				if ngstarccs:
					if 'NGSTAR_profiles.tab' in p:
						if st != 'ST':
							ngstardic[st] = linesplit[8]
						else:
							if 'NG-STAR_CC' not in linesplit:
								print('## NG-STAR CCs not in the database. Update the database with -u including the CCs file with -cc. Exiting...')
								sys.exit()
				profilesDB[profile] = st
		profilesdic[paths[n]] = profilesDB
		mosaicpenadic = {}
		if mosaic_pena:
			with open(path+'/penA_mosaics.tab', 'r') as mosaics:
				for line in mosaics:
					if 'penA_allele' not in line:
						linesplit=line.strip().split()
						mosaicpenadic[linesplit[0]] = linesplit[1]
	return profilesdic, ngstardic, mosaicpenadic

def read_ngstar_ccs(file):
	ngstarCCs = {}
	with open(file, 'r') as tab:
		for line in tab:
			if not line.startswith('Sequence Type'):
				linesplit = line.rstrip().split(',')
				ngstarCCs[linesplit[0]] = linesplit[1]
	return ngstarCCs

def integrate_ngstar_ccs(db_path, ngstarfile, ccsfile):
	ngstarccsfile = read_ngstar_ccs(ccsfile)
	ngstarprofiles = [] 
	with open(db_path+'/'+ngstarfile, 'r') as profiles:
		for line in profiles:
			editedline = line.rstrip()
			if 'ST' in editedline:
				if 'NG-STAR_CC' in editedline:
					print('## NG-STAR CCs already in the database. Exiting...')
					sys.exit()
				else:
					editedline+='\tNG-STAR_CC'
					ngstarprofiles.append(editedline)
			else:
				linesplit = editedline.split()
				if linesplit[0] in ngstarccsfile:
					editedline+='\t'+ngstarccsfile[linesplit[0]]
				else:
					editedline+='\t-'
				ngstarprofiles.append(editedline)
		subprocess.run(['rm '+db_path+'/NGSTAR_profiles.tab'], shell=True)
		outfile = open(db_path+'/NGSTAR_profiles.tab', 'w+')
		for line in ngstarprofiles:
			outfile.write(line+'\n')
		outfile.close()

def assign_ccs_only(filename, outfilename, ngstarccs, ngstarCCsdic, column):
	if ngstarccs:
		ngstarccsvec = []
		with open(filename, 'r') as tab:
			header = next(tab)
			if '>' in header:
				print('## Input is a fasta file, a table with NG-STAR STs is required to assign CCs. Exiting...')
				sys.exit()
			sep = '\t' if '\t' in header else ','
			ngstarccsvec.append(header.rstrip()+sep+'NG-STAR_CC')
			for line in tab:
				linesplit = line.rstrip().split(sep)
				st = linesplit[int(column)-1]
				if st in ngstarCCsdic:
					cc = ngstarCCsdic[st]
				else:
					cc = '-'
				ngstarccsvec.append(sep.join(linesplit)+sep+'CC'+cc)
		# Print results
		if outfilename:
			outfilehandle = open(outfilename, 'w+')
			for v in ngstarccsvec: 
				outfilehandle.write(v+'\n')
			outfilehandle.close()
		else:
			for v in ngstarccsvec:
				print(v)
	else:
		print('## Please, include -c to get the NG-STAR CCs')
	sys.exit()

def assign_sts_only(filename, outfilename, schemes, profilesDB, MLSTorder, NGSTARorder, NGMASTorder, ngstarccs, ngstarCCsdic, mosaic_pena, penAmosaicsdic):
	resultsvec = []
	with open(filename, 'r') as tab:
		header = next(tab)
		if '>' in header:
			print('## Input is a fasta file, a table with typing profiles is required to assign sequence types. Exiting...')
			sys.exit()
		sep = '\t' if '\t' in header else ','
		headersplit = header.rstrip().split(sep)
		all_loci = headersplit[1:] # asuming the first column contains strain names
		# Generate new header
		newheader = [headersplit[0]]
		for s in schemes:
			if s=='MLST':
				newheader.append('MLST')
				newheader.append(sep.join(MLSTorder))
			elif s=='NG-MAST':
				newheader.append('NG-MAST')
				newheader.append(sep.join(NGMASTorder))
			elif s=='NG-STAR':
				newheader.append('NG-STAR')
				newheader.append(sep.join(NGSTARorder))
				if ngstarccs:
					newheader.append('NG-STAR_CC')
				if mosaic_pena:
					newheader.append('penA_mosaic_type')
		resultsvec.append(sep.join(newheader))
		# Get all loci and assign STs	
		save_all = {}
		for line in tab:
			wholelinesplit = line.rstrip().split(sep)
			linesplit = wholelinesplit[1:]
			for index, locus in enumerate(all_loci):
				save_all[locus] = linesplit[index]
			MLSTprof = []
			NGSTARprof = []
			NGMASTprof = []
			for l in MLSTorder:
				if l in save_all:
					MLSTprof.append(save_all[l])
			for l in NGSTARorder:
				if l in save_all:
					NGSTARprof.append(save_all[l])
			for l in NGMASTorder:
				if l in save_all:
					NGMASTprof.append(save_all[l])
			mlst = profilesDB['MLST'].get('\t'.join(MLSTprof), '-')
			ngmast = profilesDB['NGMAST'].get('\t'.join(NGMASTprof), '-')
			ngstar = profilesDB['NGSTAR'].get('\t'.join(NGSTARprof), '-')
			to_join = [wholelinesplit[0]]
			for s in schemes:
				if s=='MLST':
					to_join.append(mlst)
					to_join.append(sep.join(MLSTprof))
				elif s=='NG-MAST':
					to_join.append(ngmast)
					to_join.append(sep.join(NGMASTprof))
				elif s=='NG-STAR':
					to_join.append(ngstar)
					to_join.append(sep.join(NGSTARprof))
					if ngstarccs:
						if ngstar in ngstarCCsdic:
							cc = 'CC'+ngstarCCsdic[ngstar]
						else:
							cc = '-'
						to_join.append(cc)
					if mosaic_pena:
						penA_al = NGSTARprof[0]
						if penA_al in penAmosaicsdic:
							penmos = penAmosaicsdic[penA_al]
						else:
							penmos = '-'
						to_join.append(penmos)
			resultsvec.append(sep.join(to_join))
		# Print results
		if outfilename:
			outfilehandle = open(outfilename, 'w+')
			for r in resultsvec: 
				outfilehandle.write(r+'\n')
			outfilehandle.close()
		else:
			for r in resultsvec:
				print(r)
	sys.exit()

def revcomp(seq):
	basesDic = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
	seqlist = list(seq)
	revseq = list(reversed(seqlist))
	comp = [basesDic[x] for x in revseq]
	return ''.join(comp)

def report_profile(order, results):
	profile = ''
	for i in order:
		# check if multiple alleles (i.e. 23S)
		check_copies = results[i]
		if len(check_copies)>1:
			unique_copies = '_'.join(list(pd.unique(results[i])))
			if profile == '':
				profile = unique_copies
			else:
				profile += '\t'+unique_copies
		else:
			if profile == '':
				profile = check_copies[0]
			else:
				profile += '\t'+check_copies[0]
	return profile

def blast_newalleles(query, subject, path):
	# run makeblastdb on subject (genome)
	subject_name = subject.split('/').pop()
	makeblastdb_cmd = ['makeblastdb', '-in', subject, '-dbtype', 'nucl', '-out', 'tmp/'+subject_name]
	subprocess.call(makeblastdb_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
	# run blastn
	blast_cmd = ['blastn', '-out', 'tmp/'+subject_name+'.blastn', '-outfmt', '6', '-query', path+'/'+query, '-db', 'tmp/'+subject_name, '-evalue', '0.001', '-max_target_seqs', '1']
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
		# Only first three closest alleles remain
		clean_alleles = []
		for x in sublist[0:3]:
			newx = x.split('_')[1]
			allele_name = newx
			if 'penA' not in x:
				allele_name = newx.split('.')[0]
			clean_alleles.append(allele_name)
		joint_allele = '|'.join(clean_alleles[0:3])
		ret = ['|'.join(clean_alleles[0:3]), coords, contigloc]
	return ret

def print_newallele_seqs(gene, coords, contigloc, fasta, allout, path):
	strand = 1
	startcoord, endcoord = coords
	if startcoord > endcoord:
		strand = -1
		endcoord, startcoord = coords
	# read genome and extract locus
	genome = Fasta(fasta)
	locus = genome[contigloc][(startcoord-1):endcoord]
	if strand == -1:
		locus.seq = revcomp(locus.seq)
	# print to file
	outfile = fasta.split('/').pop()+'.'+gene+'.fasta'
	if allout:
		with open(path+'/'+outfile, 'w') as out:
			out.write('>'+gene+'_'+locus.name+'_'+str(startcoord)+':'+str(endcoord)+'\n')
			out.write(locus.seq+'\n')
	return locus
	
def make_ACautomaton(path, allelesDB):
	allelesAC = ahocorasick.Automaton()
	for idx,key in enumerate(allelesDB):
		allelesAC.add_word(key, (idx, key))
	allelesAC.make_automaton()
	pickle.dump((allelesDB, allelesAC), open(path+'/'+'ST_alleles_AC.pkl', 'wb'))
 
def ac_fast(fname, seq, order, allelesDB, allelesAC, genogroups, out_path, PORout_results, TBPBout_results):
	results = {}
	for i in order:
		results[i] = ['-']
	for end_index, (insert_order, original_value) in allelesAC.iter(seq):
		gene = allelesDB[original_value].gene
		allele = allelesDB[original_value].allele
		direction = allelesDB[original_value].revcomp
		if gene in order:
			if type(gene) is not list:
				results[gene] = [allele]
			else:
				results[gene].append(allele)
			# save POR and TBPB sequences if genogroups are requested
			if genogroups:
				save_ngmast_genes(fname, original_value, gene, allele, direction, out_path, PORout_results, TBPBout_results)
	return results

def save_ngmast_genes(fname, seq, gene, allele, direction, out_path, PORout_results, TBPBout_results):
	if direction == 1: #revcomp
		sequence = revcomp(seq)
	else:
		sequence = seq
	if gene=='POR':
		PORrecord = SeqRecord(seq=Seq(sequence),id=fname,description=gene+'_'+allele) 
		PORout_results.append(PORrecord)
		#SeqIO.write(PORrecord, PORout, 'fasta')
	if gene=='TBPB':
		TBPBrecord = SeqRecord(seq=Seq(sequence),id=fname,description=gene+'_'+allele)
		TBPBout_results.append(TBPBrecord)
		#SeqIO.write(TBPBrecord, TBPBout, 'fasta')

def write_ngmast_genes(filename, results):
	for r in results:
		SeqIO.write(r, filename, 'fasta')

def align_sequences(infilename, outfilename):
	#use maxiters=2 for speed
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
				if seqa[x]!=seqb[x] and seqa[x] in ['A','C','G','T'] and seqb[x] in ['A','C','G','T']:
					count+=1
			dist[recorda.id][recordb.id]=count
	return dist

def prepare_files_for_genogroups(out_path):
	align_sequences('POR_out.fas', 'POR_out.aln')
	align_sequences('TBPB_out.fas', 'TBPB_out.aln')

def calculate_genogroups(out_path, PORout_results, TBPBout_results, ngmastClusters):
	# Write PORout and TBPBout results to output files
	PORout_file = open(out_path+'/POR_out.fas', 'w+')
	TBPBout_file = open(out_path+'/TBPB_out.fas', 'w+')
	write_ngmast_genes(PORout_file, PORout_results)
	write_ngmast_genes(TBPBout_file, TBPBout_results)
	PORout_file.close()
	TBPBout_file.close()
	print('## Calculating NG-MAST genogroups...')
	prepare_files_for_genogroups(out_path)
	# Calculate distances
	por_dist = calculate_distances('POR_out.aln')
	tbpb_dist = calculate_distances('TBPB_out.aln')
	# calculate clusters
	clusters=[]
	for st in ngmastClusters:
		if st != '-' and st != 'multiple':
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
	#print(clusters)
	for cluster in clusters:
		g=set([])
		for c in cluster:
			g.add(c[-1])
		sts=[]
		for c in g:
			sts.append([len(ngmastClusters[str(c)]), str(c)])
		sts.sort()
		sts.reverse()
		if sts[0][1]!='-':
			genogroup_name='G'+str(sts[0][1])
			stssorted = []
			for s in sts:
				genogroups[s[1]]=genogroup_name
				stssorted.append(s[1]+'('+str(s[0])+')')
			print(genogroup_name, '; '.join(map(str,stssorted)), file=output, sep='\t')
	output.close()
	genopergenome = {}
	for n in ngmastClusters:
		for sample in ngmastClusters[n]:
			if n != '-':
				genopergenome[sample[0]] = genogroups[n]
			else:
				genopergenome[sample[0]] = '-'
	return genopergenome

def create_header(schemes, genogroups, ngstarccs, mosaic_pena, MLSTorder, NGSTARorder, NGMASTorder):
	# Create and print output header
	header = ['strain']
	order = []
	for s in schemes:
		if s=='MLST':
			tmporder = ['MLST']+MLSTorder
			header+=tmporder
			order+=MLSTorder
		elif s=='NG-STAR':
			tmporder = ['NG-STAR']+NGSTARorder
			header+=tmporder
			order+=NGSTARorder
			if ngstarccs:
				header+=['NG-STAR_CC']
			if mosaic_pena:
				header+=['penA_mosaic_type']
		elif s=='NG-MAST':
			tmporder = ['NG-MAST']+NGMASTorder
			header+=tmporder
			order+=NGMASTorder
			if genogroups:
				header+=['Genogroup']
	return header, order

def process_files(args):
	fpath, order, MLSTorder, NGSTARorder, NGMASTorder, profilesDB, allelesDB, allelesAC, schemes, \
    ngstarCCsdic, penAmosaicsdic, genogroups, db_path, PORout_results, TBPBout_results, allout, blast_new_alleles, out_path, \
    new_alleles, indices_new_alleles = args
	fname = fpath.split('/').pop()
	blast_hit = 0
	if not fpath.endswith('.fasta') or fpath.endswith('.fas') or fpath.endswith('.fa'):
		print('## Input file is not recognised as a fasta file. Please, include fasta files ending in .fasta, .fas or .fa')
		sys.exit()
	with open(fpath, 'r') as fasta:
		concat = ''
		for record in SeqIO.parse(fasta, 'fasta'):
			concat += record.seq
	#AC algorithm:
	resultsDB = ac_fast(fname, str(concat), order, allelesDB, allelesAC, genogroups, out_path, PORout_results, TBPBout_results)
	# run Blast to get closest alleles to missing loci
	if blast_new_alleles:
		for x in resultsDB:
			if resultsDB[x] == ['-']:
				query_file = x+'.fas'
				subject_file = fpath
				blastout = blast_newalleles(query_file, subject_file, db_path)
				if len(blastout)>1:
					blast_hit = 1
					closest_allele, coords, contigloc = blastout
					resultsDB[x] = [closest_allele]
					extractedallele = print_newallele_seqs(x, coords, contigloc, fpath, allout, out_path)
					if str(extractedallele) in new_alleles:
						item_index = new_alleles[str(extractedallele)]
					else:
						callele = x+'_'+closest_allele
						if callele in indices_new_alleles:
							itemi = indices_new_alleles[callele]+1
							indices_new_alleles[callele] = itemi
							item_index = closest_allele+'-'+str(itemi)
							new_alleles[str(extractedallele)] = item_index
						else:
							indices_new_alleles[callele] = 1
							item_index = closest_allele+'-1'
							new_alleles[str(extractedallele)] = item_index
					resultsDB[x] = [item_index]
	else:
		blast_hit = 0
	# get profile and assign ST
	st_list = {}
	ngmastcl = {}
	ngmast = ''
	for s in schemes:
		if s=='MLST':
			prof = report_profile(MLSTorder, resultsDB)
			st = profilesDB['MLST'].get(prof, '-')
			st_list['MLST'] = st+'\t'+prof
		elif s=='NG-STAR':
			if '|' not in resultsDB['penA'][0] and resultsDB['penA'][0] != '-':
				checkpena = resultsDB['penA'][0].split('.')
				checkpena2 = checkpena[1].split('-')
				if len(checkpena2[0])==2:
					checkpena2[0] = checkpena2[0]+'0'
					resultsDB['penA'][0] = checkpena[0]+'.'+'-'.join(checkpena2)
			prof = report_profile(NGSTARorder, resultsDB)
			st = profilesDB['NGSTAR'].get(prof, '-')
			if ngstarCCsdic:
				if st=='-':
					st_list['NG-STAR'] = st+'\t'+prof+'\t'+'-'
				else:
					if ngstarCCsdic[st]=='Ungroupable':
						st_list['NG-STAR'] = st+'\t'+prof+'\tUngroupable'
					elif ngstarCCsdic[st]=='-':
						st_list['NG-STAR'] = st+'\t'+prof+'\t-'
					else:
						st_list['NG-STAR'] = st+'\t'+prof+'\tCC'+ngstarCCsdic[st]
			else:		
				st_list['NG-STAR'] = st+'\t'+prof
			if penAmosaicsdic:
				penA = prof.split('\t')[0]
				st_list['NG-STAR'] += '\t'+penAmosaicsdic[penA]
		elif s=='NG-MAST':
			prof = report_profile(NGMASTorder, resultsDB)
			st = profilesDB['NGMAST'].get(prof, '-')
			st_list['NG-MAST'] = st+'\t'+prof
			if genogroups:
				ngmastcl[st] = [fname, prof]
				ngmast = st
	return fname, [st_list, ngmast, ngmastcl, new_alleles, indices_new_alleles]

def access_results(results, genogroups):
	finalresults = {}
	ngmastClusters = {}
	for f, all_results in results:
		st_list, ngmast, ngmastcl, new_alleles, indices_new_alleles = all_results
		finalresults[f] = st_list
		# Handle genogroups if requested
		if genogroups:
			if ngmast not in ngmastClusters:
				ngmastClusters[ngmast] = [ngmastcl[ngmast]]
			else:
				ngmastClusters[ngmast].append(ngmastcl[ngmast])
	return finalresults, ngmastClusters

def print_results(header, finalresults, genogroups, genopergenome, out_path, outfile):
	# Print header
	if outfile:
		outfilehandle = open(out_path+'/'+outfile, 'w+')
		outfilehandle.write('\t'.join(header)+'\n')
	else:
		outfilehandle = None
		print('\n'+'\t'.join(header))
	# Print by line
	for f in finalresults:
		outline = ''
		for s in finalresults[f]:
			if genogroups and s=='NG-MAST':
				finalresults[f][s]+='\t'+genopergenome[f]
			outline += '\t'+finalresults[f][s]
		if outfile:	
			outfilehandle.write(f+outline+'\n')
		else:
			print(f+outline)

def clean_tmp_files(out_path):
	subprocess.call(['rm', '-r', out_path+'/tmp'])
	# Clean *.fai files generated in the same directory as the input files
	for r in glob.glob(out_path+"/*.fai"):
		subprocess.call(['rm', r])


