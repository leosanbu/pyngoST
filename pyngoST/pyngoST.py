#!/usr/bin/env python3
import os
import sys
import glob
import time
import pickle
import subprocess
import argparse as arg
from argparse import RawTextHelpFormatter
import concurrent.futures
from pyngoST_utils import * 

## Start the timer
start_time = time.time()

## Script parameters
parser = arg.ArgumentParser(prog="pyngoST",
	formatter_class=RawTextHelpFormatter,
	description='pyngoST: fast, simultaneous and accurate multiple sequence typing of Neisseria gonorrhoeae genome collections\n'
	'\nCitation:\n'
	'    Sanchez-Buso L, Sanchez-Serrano A, Golparian D and Unemo M.\n'
	'    pyngoST: fast, simultaneous and accurate multiple sequence typing of Neisseria gonorrhoeae genome collections.\n'
	'    Preprint: https://www.biorxiv.org/content/10.1101/2023.10.23.563537v2\n'
	'    GitHub: https://github.com/leosanbu/pyngoST\n',
	usage = '%(prog)s [options]')

parser.add_argument('-i', '--input', help='Input files (fasta or tab/csv)', required=False, nargs='+')
parser.add_argument('-r', '--read_file', help='File containing the paths to the input files', required=False)
parser.add_argument('-s', '--schemes', help='Typing schemes to query separated by commas (options: NG-STAR, MLST, NG-MAST) (default=NG-STAR)', required=False, default='NG-STAR')
parser.add_argument('-g', '--genogroups', help='Calculate NG-MAST genogroups from NG-MAST types (default=False)', required=False, action='store_true')
parser.add_argument('-c', '--ngstarccs', help='Include NG-STAR CCs in output table (default=False)', required=False, action='store_true')
parser.add_argument('-m', '--mosaic_pena', help='Report if the penA allele sequence is a mosaic or semimosaic (default=False)', required=False, action='store_true')
parser.add_argument('-b', '--blast_new_alleles', help='Use blastn to find the closest alleles to new ones (default=False)', required=False, action='store_true')
parser.add_argument('-a', '--alleles_out', help='Print fasta files with new alleles (optional, default=False)', required=False, action='store_true')
parser.add_argument('-q', '--out_path', help='Path used to save output files (default: current directory)', required=False)
parser.add_argument('-o', '--out_filename', help='Name of file to print the results table to (optional, default=screen output)', required=False)
parser.add_argument('-y', '--only_assignccs', help='Only assign CCs from a table with NG-STAR STs. Indicate as value the number of the column that contains the STs (optional, default=None)', required=False)
parser.add_argument('-z', '--only_assignsts', help='Only assign STs from a table with NG-STAR, MLST and/or NG-MAST profiles (optional, default=None)', required=False, action='store_true')
parser.add_argument('-t', '--num_threads', help='Number of processes to use for computation (optional, default=1)', required=False)
# Options for the database
parser.add_argument('-p', '--path', help='Path to database containing MLST/NG-STAR alleles and profiles. If not available, use -d to create an updated database', required=False)
parser.add_argument('-d', '--download_db', help='Download updated allele and profile files and create database', required=False, action='store_true')
parser.add_argument('-n', '--db_name', help='Name of the folder that will contain the database in case a download is requested with -d (default=allelesDB in working directory)', required=False)
parser.add_argument('-u', '--update', help='Only recreate the database pickle file', required=False, action='store_true')
parser.add_argument('-cc', '--ngstarccsfile', help='File with the NG-STAR clonal complexes (NG-STAR CCs) database (csv) to integrate to NG-STAR profiles', required=False)
if len(sys.argv)==1:
	parser.print_help(sys.stderr)
	sys.exit(0)
arg = parser.parse_args()

if __name__ == '__main__':

	## Print script header
	print("\n## pyngoST: multiple sequence typing of Neisseria gonorrhoeae for large assembly collections\n")

	## Get arguments
	args = get_args(arg)

	## Extract variables from the args dictionary
	db_path = args['db_path']
	db_name = args['db_name']
	out_path = args['out_path']
	schemes = args['schemes']
	ngstarccs = args['ngstarccs']
	genogroups = args['genogroups']
	mosaic_pena = args['mosaic_pena']
	blast_new_alleles = args['blast_new_alleles']
	allout = args['allout']
	update = args['update']
	download = args['download']
	ccsfile = args['ccsfile']
	outfile = args['outfile']
	num_threads = int(args['num_threads'])
	only_assignccs = args['only_assignccs']
	only_assignsts = args['only_assignsts']

	## Download/Update/Load pickled database or only assign CCs
	if download:
		download_db(db_path, db_name, ccsfile)
	elif update:
		update_db(db_path, ccsfile)
	else:
		## Load dbs
		allelesDB, allelesAC = load_db(db_path)
		profilesDB, ngstarCCsdic, penAmosaicsdic = read_profiles(db_path, ngstarccs, mosaic_pena)
		## Read input files
		filelist = get_input(arg.input, arg.read_file)
		## Define profiles and loci order 
		MLSTorder = ['abcZ', 'adk', 'aroE', 'fumC', 'gdh', 'pdhC', 'pgm']
		NGSTARorder = ['penA', 'mtrR', 'porB', 'ponA', 'gyrA', 'parC', '23S']
		NGMASTorder = ['POR', 'TBPB']
		## Create header and set order of schemes
		header, order = create_header(schemes, genogroups, ngstarccs, mosaic_pena, MLSTorder, NGSTARorder, NGMASTorder)
		if only_assignccs:
			input_file = filelist[0]
			column = only_assignccs
			assign_ccs_only(input_file, outfile, ngstarccs, ngstarCCsdic, column)
		if only_assignsts:
			if genogroups:
				print('## Genogroups can only be obtained from fasta files. Exiting...')
				sys.exit()
			if blast_new_alleles:
				print('## Blast to obtain close alleles can only be perform with fasta files. Exiting...')
				sys.exit()
			if allout:
				print('## New allele sequences can only be saved with fasta files. Exiting...')
				sys.exit()
			input_file = filelist[0]
			assign_sts_only(input_file, outfile, schemes, profilesDB, MLSTorder, NGSTARorder, NGMASTorder, ngstarccs, ngstarCCsdic, mosaic_pena, penAmosaicsdic)

	## Print schemes requested
	sc = ','.join(schemes)
	sc += ',NG-STAR CCs' if ngstarccs else ''
	sc += ',NG-MAST Genogroups' if genogroups else ''
	print("## Schemes requested:", sc)

	## Create tmp/ folder for intermediate operations
	if not 'tmp' in os.listdir(out_path):
		subprocess.run(['mkdir '+out_path+'/tmp'], shell=True)

	## Run the process_files function
	print('## Number of processes:', num_threads)
	print('## Finding matches...')
	new_alleles = {}
	indices_new_alleles = {}
	PORout_results = []
	TBPBout_results = []
	args_list = [(f, order, MLSTorder, NGSTARorder, NGMASTorder, profilesDB, allelesDB, allelesAC, schemes,
				ngstarCCsdic, penAmosaicsdic, genogroups, db_path, PORout_results, TBPBout_results, allout, blast_new_alleles, out_path,
				new_alleles, indices_new_alleles) for f in filelist]

	with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as executor:
		results = executor.map(process_files, args_list)
		
	## Access results
	finalresults, ngmastClusters = access_results(results, genogroups)

	## Calculate genogroups if requested
	genopergenome = calculate_genogroups(out_path, PORout_results, TBPBout_results, ngmastClusters) if genogroups else None

	## Print final results
	print_results(header, finalresults, genogroups, genopergenome, out_path, outfile)

	## Clean output files
	clean_tmp_files(out_path)

	## Calculate and print the elapsed time
	elapsed_time = time.time() - start_time
	elapsed_time_minutes = elapsed_time/60
	print("\n## Computation time: %.2f seconds, %.2f minutes" % (elapsed_time, elapsed_time_minutes))


	
