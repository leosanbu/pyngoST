import os
import sys
import pickle
import subprocess
import argparse as arg
import concurrent.futures
from pyngoST_functions import * 

parser = arg.ArgumentParser(description='pyngoST: high throughput retrieval of sequence types from N. gonorrhoeae assembly collections', usage = '%(prog)s [options]')
parser.add_argument('-i', '--input', help='Assembly files', required=False, nargs='+') # add multiple assemblies
parser.add_argument('-r', '--read_file', help='Read file containing the paths to the input files', required=False) # use if too many input assemblies or several input paths
parser.add_argument('-s', '--schemes', help='Typing schemes to query: NG-STAR, MLST, cgMLST (default=MLST,NG-STAR)', required=False, default='MLST,NG-STAR')
parser.add_argument('-c', '--ngstarcc', help='Include NG-STAR clonal complexes (NG-STAR CCs)', required=False)
parser.add_argument('-p', '--path', help='Path to database containing alleles and profiles', required=True)
parser.add_argument('-a', '--alleles_out', help='Print fasta files with new alleles (optional, default: False)', required=False, action='store_true')
parser.add_argument('-o', '--out_filename', help='Name of file to print output to (optional, default: screen output)', required=False)
parser.add_argument('-u', '--update', help='Re-create .pkl file after updating the database', required=False, action='store_true')
arg = parser.parse_args()

# Get arguments #
out_path = os.getcwd()
db_path = arg.path
schemes = arg.schemes.split(',') if arg.schemes else ['MLST', 'NG-STAR']
ccsfile = arg.ngstarcc if arg.ngstarcc else False
allout = True if arg.alleles_out else False
outfile = arg.out_filename if arg.out_filename else False
update = True if arg.update else False
filelist = getInput(arg.input, arg.read_file)

if update:
	# Only update the pickle file if requested
	allelesDB = readAlleles(db_path)
	updatePKL(db_path, allelesDB)
	print('Database (pickle file) updated in '+db_path+'/'+'ST_alleles_AC.pkl')
	sys.exit()
else:
	# Load pickled database - dictionary and automaton    
	allelesDB, allelesAC = pickle.load(open(db_path+'/'+'ST_alleles_AC.pkl', 'rb'))

## Main script ##

# order of loci 
MLSTorder = ['abcZ', 'adk', 'aroE', 'fumC', 'gdh', 'pdhC', 'pgm']
NGSTARorder = ['penA', 'mtrR', 'porB', 'ponA', 'gyrA', 'parC', '23S']

# read profiles file
#if 'MLST' and 'NG-STAR' in schemes:
profilesDB = readProfiles(db_path)
#elif schemes==['MLST']:
#
#elif schemes==['NG-STAR']:


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

if outfile:
	outfilehandle = open(out_path+'/'+outfile, 'w+')
	outfilehandle.write('\t'.join(header)+'\n')
else:
	outfilehandle = None
	print('\t'.join(header))

# create tmp/ folder for intermediate operations
subprocess.call(['mkdir', 'tmp'])

# process files
for f in filelist:
	processFiles(f, order, MLSTorder, NGSTARorder, profilesDB, allelesDB, allelesAC, schemes, ngstarCCs, db_path, allout, out_path, outfile, outfilehandle)

# clean output files
subprocess.call(['rm', '-r', 'tmp'])
	




