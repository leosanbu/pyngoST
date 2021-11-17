import os
import sys
import pickle
import subprocess
import argparse as arg
import concurrent.futures
from pyngoST_functions import * 

parser = arg.ArgumentParser(description='pyngoST: high throughput molecular typing of N. gonorrhoeae genome assemblies', usage = '%(prog)s [options]')
parser.add_argument('-i', '--input', help='Assembly files in fasta format', required=False, nargs='+')
parser.add_argument('-r', '--read_file', help='File containing the paths to the input files', required=False)
parser.add_argument('-s', '--schemes', help='Typing schemes to query separated by commas (options: NG-STAR, MLST, NG-MAST) (default=NG-STAR)', required=False, default='NG-STAR')
parser.add_argument('-c', '--ngstarcc', help='File with the NG-STAR clonal complexes (NG-STAR CCs) database to include NG-STAR CCs in the results table', required=False)
parser.add_argument('-g', '--genogroups', help='Calculate NG-MAST genogroups from NG-MAST types (default=False)', required=False, action='store_true')
parser.add_argument('-p', '--path', help='Path to database containing MLST/NG-STAR alleles and profiles', required=True)
parser.add_argument('-a', '--alleles_out', help='Print fasta files with new alleles (optional, default: False)', required=False, action='store_true')
parser.add_argument('-o', '--out_filename', help='Name of file to print the results table to (optional, default: screen output)', required=False)
parser.add_argument('-u', '--update', help='Re-create .pkl file after updating the database', required=False, action='store_true')
arg = parser.parse_args()

#####################
### Get arguments ###
#####################

out_path = os.getcwd()
db_path = arg.path
schemes = arg.schemes.split(',') if arg.schemes else ['NG-STAR']
ccsfile = arg.ngstarcc if arg.ngstarcc else False
genogroups = True if arg.genogroups else False
allout = True if arg.alleles_out else False
outfile = arg.out_filename if arg.out_filename else False
update = True if arg.update else False
filelist = getInput(arg.input, arg.read_file)

###############################
### Update pickled database ###
###############################

if update:
	# Only update the pickle file if requested
	allelesDB = readAlleles(db_path)
	updatePKL(db_path, allelesDB)
	print('Database (pickle file) updated in '+db_path+'/'+'ST_alleles_AC.pkl')
	sys.exit()
else:
	# Load pickled database - dictionary and automaton    
	allelesDB, allelesAC = pickle.load(open(db_path+'/'+'ST_alleles_AC.pkl', 'rb'))

###################
### Main script ###
###################

# order of loci 
MLSTorder = ['abcZ', 'adk', 'aroE', 'fumC', 'gdh', 'pdhC', 'pgm']
NGSTARorder = ['penA', 'mtrR', 'porB', 'ponA', 'gyrA', 'parC', '23S']
NGMASTorder = ['POR', 'TBPB']

# read profiles file
profilesDB = readProfiles(db_path)

# read NG-STAR CCs
if ccsfile:
	ngstarCCs = readNGSTARCCs(ccsfile)
else:
	ngstarCCs = None

# create and print output header
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
		if ngstarCCs:
			header+=['NG-STAR_CC']
	elif s=='NG-MAST':
		tmporder = ['NG-MAST']+NGMASTorder
		header+=tmporder
		if genogroups:
			header+=['Genogroup']

# create tmp/ folder for intermediate operations
if not 'tmp' in os.listdir(out_path):
	subprocess.run(['mkdir tmp'], shell=True)

# run NG-MAST and genogroups if requested
if 'NG-MAST' in schemes:
	print('## Running NGMASTER...')
	ngmastDic, ngmastClusters = run_ngmaster(filelist, genogroups)
	if genogroups:
		print('## Splitting and aligning POR.fas and TBPB.fas files...')
		prepare_files_for_genogroups(out_path)
		print('## Calculating NG-MAST genogroups...')
		genopergenome = calculateGenogroups(out_path, ngmastClusters)
	else:
		genopergenome = None
else:
	ngmastDic = None
	genopergenome = None

# open ouput file and print header
if outfile:
	outfilehandle = open(out_path+'/'+outfile, 'w+')
	outfilehandle.write('\t'.join(header)+'\n')
else:
	outfilehandle = None
	print('\t'.join(header))

# process files to call MLST/NG-STAR and add info on NG-MAST/genogroups if requested
for f in filelist:
	processFiles(f, order, MLSTorder, NGSTARorder, profilesDB, allelesDB, allelesAC, schemes, 
		ngstarCCs, ngmastDic, genopergenome, db_path, allout, out_path, outfile, outfilehandle)

# clean output files
subprocess.call(['rm', '-r', 'tmp'])
	




