# pyngoST: multiple sequence typing of _Neisseria gonorrhoeae_ assembly collections

`pyngoST` unifies molecular typing for _Neisseria gonorrhoeae_. By providing genome assemblies in fasta format, `pyngoST` can extract up to three typing schemes (NG-STAR, MLST and NG-MAST) and two modifications (NG-STAR Clonal Complexes and NG-MAST Genogroups). These schemes are detailed below:
* **NG-STAR**: _N. gonorrhoeae_ Sequence Typing for Antimicrobial Resistance, hosted at Public Health Agency of Canada, National Microbiology Laboratory [NG-STAR Canada](https://ngstar.canada.ca/). This is a typing scheme published by [Demczuk et al. 2017](https://doi.org/10.1128/JCM.00100-17) that targets 7 genes associated to cephalosporin, macrolides and fluoroquinolones resistance: _**penA**, **mtrR**, **porB**, **ponA**, **gyrA**, **parC**_ and _**23S rRNA**_.
* **NG-STAR CCs** (Clonal Complexes), published by [Golparian et al, 2021](https://doi.org/10.1093/jac/dkaa552). This scheme groups NG-STAR sequence types (STs) into Clonal Complexes (CCs) for a better fit of the typing scheme with the population structure of antimicrobial resistant (AMR) lineages.
* **MLST** (Multi-Locus Sequence Typing), hosted at [PubMLST Neisseria](https://pubmlst.org/bigsdb?db=pubmlst_mlst_seqdef). This scheme was published by [Maiden et al, 1998](https://www.pnas.org/doi/full/10.1073/pnas.95.6.3140) and targets 7 housekeeping genes of the _Neisseria_ genus: _**abcZ**, **adk**, **aroE**, **fumC**, **gdh**, **pdhC**_ and _**pgm**_.
* **NG-MAST v2** (_N. gonorrhoeae_ Multi-Antigen Sequence Typing), hosted at [PubMLST Neisseria](https://pubmlst.org/bigsdb?db=pubmlst_mlst_seqdef). This scheme was published by [Martin et al. 2004](https://doi.org/10.1086/383047) and targets 2 genes encoding rapidly-evolving surface antigens: _**porB**_ and _**tbpbB**_.
* **NG-MAST Genogroups**, as described by [Harris et al, 2018](https://doi.org/10.1016/S1473-3099(18)30225-1). As NG-MAST loci are highly-variable, they poorly match the population structure of _N. gonorrhoeae_. Genogroups were described to group NG-MAST types, however they are calculated on the specific dataset under study and thus, are not comparable among studies.

`pyngoST` is written in python3 and uses the [_pyahocorasick_](https://github.com/WojciechMula/pyahocorasick) library for fast simultaneous multi-loci search in _N. gonorrhoeae_ genomes. 

First, a database is built that contains updated fasta files with all the known alleles for 7 NG-STAR, 7 MLST and 2 NG-MAST loci as well as the known NG-STAR, MLST and NG-MAST profiles up to the building of the database. Then, a simultaneous screening of all loci is performed on the input assembly files (in FASTA) format using the [Aho-Corasick string-searching algorithm](https://en.wikipedia.org/wiki/Aho%E2%80%93Corasick_algorithm).

## Installation

I recommend to create a python3 virtual environment and install pyngoST using `pip`.

```
## Create and activate a virtual environment (i.e. `venv` or any other name):
python3 -m venv venv
source venv/bin/activate

## Install pyngoST using pip install on the latest distribution
pip install pyngoST
```
Activate the virtual environment whenever you want to use `pyngoST`. Exit the virtual environment by running `deactivate`:
```
## Activate virtualenv
source venv/bin/activate

## Run pyngSTar

## Exit virtualenv
deactivate
```

## Bulding the database

You can download updated allele and profiles files for the three schemes using `-d`. NG-STAR alleles and profiles will be downloaded and adapted from [NG-STAR Canada](https://ngstar.canada.ca/)). MLST and NG-MAST alleles and profiles will be downloaded from [PubMLST Neisseria](https://pubmlst.org/bigsdb?db=pubmlst_mlst_seqdef). By default, it will create the database under the `allelesDB` folder if it does not exist in the working directory. You can provide a database name with `-n`. If a CSV file containing NG-STAR CCs is included with `-cc`, CCs will be integrated into the database.
```
## Create the database as 'allelesDB'
pyngoST.py -d

## Create the database with another name
pyngoST.py -d -n othernameDB

## Create the database with another name and include NG-STAR CCs
pyngoST.py -d -n othernameDB -cc NGSTAR_CCs.csv
```
During the creation of the database, a dictionary containing all alleles of all loci is saved as a pickle file, which will be loaded every time `pyngoST` is run.

## Updating the database

If you choose to manually modify any of the files in the database and/or you want to include NG-STAR CCs (not previously included when the database was downloaded), you need to update the existing pickle file in the database with `-u`. Use `-p` to specify the path to the existing database.
```
## Update the pickle file in the database
pyngoST.py -u -p /path/to/allelesDB

## Incorporate NG-STAR CCs and update the pickle file in the database
pyngoST.py -u -p /path/to/allelesDB -cc NGSTAR_CCs.csv
```

## Available options

```
usage: pyngoST.py [options]

pyngoST: fast, simultaneous and accurate multiple sequence typing of Neisseria gonorrhoeae genome collections

Citation:
    Sanchez-Buso L, Sanchez-Serrano A, Golparian D and Unemo M.
    pyngoST: fast, simultaneous and accurate multiple sequence typing of Neisseria gonorrhoeae genome collections.
    GitHub: https://github.com/leosanbu/pyngoST

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT [INPUT ...], --input INPUT [INPUT ...]
                        Input files (fasta or tab/csv)
  -r READ_FILE, --read_file READ_FILE
                        File containing the paths to the input files
  -s SCHEMES, --schemes SCHEMES
                        Typing schemes to query separated by commas (options: NG-STAR, MLST, NG-MAST) (default=NG-STAR)
  -g, --genogroups      Calculate NG-MAST genogroups from NG-MAST types (default=False)
  -c, --ngstarccs       Include NG-STAR CCs in output table (default=False)
  -b, --blast_new_alleles
                        Use blastn to find the closest alleles to new ones (default: False)
  -a, --alleles_out     Print fasta files with new alleles (optional, default: False)
  -q OUT_PATH, --out_path OUT_PATH
                        Path used to save output files (default: current directory)
  -o OUT_FILENAME, --out_filename OUT_FILENAME
                        Name of file to print the results table to (optional, default: screen output)
  -y ONLY_ASSIGNCCS, --only_assignccs ONLY_ASSIGNCCS
                        Only assign CCs from a table with NG-STAR STs. Indicate as value the number of the column that contains the STs (optional, default: None)
  -z, --only_assignsts  Only assign STs from a table with NG-STAR, MLST and/or NG-MAST profiles (optional, default: None)
  -t NUM_THREADS, --num_threads NUM_THREADS
                        Number of processes to use for computation (optional, default: 1)
  -p PATH, --path PATH  Path to database containing MLST/NG-STAR alleles and profiles. If not available, use -d to create an updated database
  -d, --download_db     Download updated allele and profile files and create database
  -n DB_NAME, --db_name DB_NAME
                        Name of the folder that will contain the database in case a download is requested with -d (default=allelesDB in working directory)
  -u, --update          Only recreate the database pickle file
  -cc NGSTARCCSFILE, --ngstarccsfile NGSTARCCSFILE
                        File with the NG-STAR clonal complexes (NG-STAR CCs) database (csv) to integrate to NG-STAR profiles
  
```
| Option | Description |
| :---: | :--- |
| -i | list of assembly files (* can be used to input multiple files in a directory). |
| -r | text file containing a list of paths to assembly files. |
| -s | schemes to query and report: NG-STAR,MLST,NG-STAR or any combination/order separated by commas. |
| -g | calculate NG-MAST genogroups from NG-MAST typing using the rules described in [Harris et al, 2018](https://doi.org/10.1016/S1473-3099(18)30225-1). |
| -c | report NG-STAR Clonal Complexes from NG-STAR STs in the final results table. |
| -b | blast new alleles to identify the closest among the known ones. |
| -a | save new alleles to output fasta files. |
| -o | name of output file where to save the results table. By default, results will print to the screen. |
| -q | path where to save results if different from the working directory. |
| -y | only call NG-STAR Clonal Complexes from a table of NG-STAR STs. Specify the number of column containing the STs. |
| -z | only assign STs from a table containing allelic profiles of any of NG-STAR, MLST and NG-STAR (specify with -s). |
| -t | number of threads for multithreading the main processing function from fasta files. |
| -p | path to the database file to be used by `pyngoST`. |
| -d | download allele and ST profile files and create the database. |
| -n | name of the database folder if different from `allelesDB`. |
| -u | update the database by recreating the pickle file. |
| -cc | input CSV file containing the NG-STAR Clonal Complexes database. |

## Running pyngoST from fasta files

The main feature of `pyngoST` is to extract molecular typing information from _N. gonorrhoeae_ assembly files. Input files can be indicated with `-i` directly or inside a text file, which can be read by `pyngoST`using `-r`. The path to the database is specified with `-p` and the typing schemes requested with `-s` (any of NG-STAR, NG-MAST and/or MLST in any order). Examples are shown below, including other options.

By default, if `-s` is not specified, only NG-STAR STs will be reported. Use `-r filelist.txt` with a text containing the paths to multiple fasta files instead of `-i *.fasta`.
```
pyngoST.py -i *.fasta -p /path/to/allelesDB 
```

Any of the NG-STAR, MLST and NG-MAST schemes can be requested with `-s` in any order. This order will be also used for the output table. NG-STAR CCs can be included into the reported table with `-c` and NG-MAST genogroups can be calculated with `-g`. Please note that genogroups calculation is computationally expensive for large datasets.
```
pyngoST.py -i *.fasta -p /path/to/allelesDB -s MLST,NG-MAST,NG-STAR -c -g
```

The output table for `-s MLST,NG-MAST,NG-STAR -c -g` looks like the following, including the columns `Genogroup` after the `NG-MAST` profile and `NG-STAR_CC` after the `NG-STAR` profile:
```
strain	MLST	abcZ	adk	aroE	fumC	gdh	pdhC	pgm	NG-MAST	POR	TBPB	Genogroup	NG-STAR	penA	mtrR	porB	ponA	gyrA	parC	23S	NG-STAR_CC
31919_4_101.fasta	8135	59	39	67	156	188	153	133	387	266	118	G387	3658	15.001	424	100	100	100	2	100	CC893
31919_4_12.fasta	11975	59	39	762	111	188	71	223	13571	7834	16	G2	84	22.001	15	100	100	100	7	100	CC352
31919_4_16.fasta	9363	126	39	67	238	148	153	133	12302	908	267	G12302	168	2.001	39	11	100	7	4	100	CC168
31919_4_17.fasta	11975	59	39	762	111	188	71	223	2	2	16	G2	84	22.001	15	100	100	100	7	100	CC352
```

Use `-o <outfilename>` to save output to file.

### Blasting new allele sequences

Blasting of new alleles to find the closest among the known ones can be requested with `-b`. These new alleles can be saved to fasta files with `-a`.
```
pyngoST.py -i *.fasta -p /path/to/allelesDB -s MLST,NG-MAST,NG-STAR -b -a
```

New alleles will be saved to individual .fasta files ending in `<input_file>.<locus_name>.fasta`. As an example, a new `abcZ` allele of the MLST scheme will be saved as `32386_1_293.fasta.abcZ.fasta` and contains the following sequence:

```
>abcZ_NODE_4_length_153129_cov_24.801849_9917:10349
TTTGATACTGTTGCCGAAGGTTTGGGCTAAATTCGCGATTTATTGCGCCGTTATCATCATGTCAGCCATGAGTTGGAAAATGGTTCGAGTGAGCTCTTATTGAAAGAGCTTAACGAATTGCAACTTGAAATCGAAGCGAAGGACGGTTGGAAACTGGATGCGGCAGTCAAGCAGACTTTGGGGGAACTCGGTTTGCCGGAAAACAAAAAAATCGGCAACCTTTCCGGCGGTCAGAAAAAACGCGTCGCCTTAGCGCAGGCTTGGGTACAGAAGCCCGATGTATTGTTACTGGACGAACCGACCAACCATTTGGATATCGACGCGATTATTTGGTTGGAAAACCTGCTCAAAGCGTTTGAAGGCAGCTTGGTCGTGATTACCCACGACCGCCGTTTTTTGGACAATATCGCCACGCGCATTGTCGAACTCGATC
```
The structure of the sequence name is as follows: 

`<locus_name>_<contig_name>_<start_coordinate>:<end_coordinate>`

where the coordinates indicate the exact position in the contig where this sequence was extracted. 

New allelic sequences and STs should be submitted to curators at [PubMLST Neisseria](https://pubmlst.org/bigsdb?db=pubmlst_mlst_seqdef) or [NG-STAR Canada](https://ngstar.canada.ca/) for assignment.

## Running pyngoST from csv/tab files

`pyngoST` has two extra functionalities that result from common procedures used in molecular epidemiology of _N. gonorrhoeae_:
* ST assignment from a table of allelic profiles of any of the NG-STAR, MLST and/or NG-MAST schemes.
* NG-STAR CCs assignment from a table containing a column with NG-STAR STs.
In these two cases, the input file is a CSV or tabular table containing the mentioned information and is specified with `-i`. Other necessary parameter is `-p` to specify the path to an updated database. The specific options to use these features are exemplified below:

### Get STs from allelic profiles

To run `pyngoST` from a table of allelic profiles use `-z`. Specify the typing schemes to look at with `-s` and include `-c` if NG-STAR CCs are also requested.
```
pyngoST.py -i testassignSTs.tsv -p /path/to/allelesDB -z -s MLST,NG-MAST,NG-STAR -c
```
The example `testassignSTs.tsv` file looks like:
```
strain	penA	mtrR	porB	ponA	gyrA	parC	23S	abcZ	adk	aroE	fumC	gdh	pdhC	pgm	POR	TBPB
31919_4_101.fasta	15.001	424	100	100	100	2	100	59	39	67	156	188	153	133	266	118
31919_4_12.fasta	22.001	15	100	100	100	7	100	59	39	762	111	188	71	223	7834	16
31919_4_16.fasta	2.001	39	11	100	7	4	100	126	39	67	238	148	153	133	908	267
31919_4_17.fasta	22.001	15	100	100	100	7	100	59	39	762	111	188	71	223	2	16
```

And the output file will contain a `MLST`, `NG-MAST`, `NG-STAR` and `NG-STAR_CCs` columns as requested:
```
strain	MLST	abcZ	adk	aroE	fumC	gdh	pdhC	pgm	NG-MAST	POR	TBPB	NG-STAR	penA	mtrR	porB	ponA	gyrA	parC	23S	NG-STAR_CC
31919_4_101.fasta	8135	59	39	67	156	188	153	133	387	266	118	3658	15.001	424	100	100	100	2	100	CC893
31919_4_12.fasta	11975	59	39	762	111	188	71	223	13571	7834	16	84	22.001	15	100	100	100	7	100	CC352
31919_4_16.fasta	9363	126	39	67	238	148	153	133	12302	908	267	168	2.001	39	11	100	7	4	100	CC168
31919_4_17.fasta	11975	59	39	762	111	188	71	223	2	2	16	84	22.001	15	100	100	100	7	100	CC352
```
### Get NG-STAR CCs from NG-STAR STs

To only call NG-STAR CCs from a table containing a column with NG-STAR STs use `-y <column_number>`, where `<column_number>` is the number of the column (starting from 1) in the table that contains the NG-STAR STs. The `-c` option is also requested to activate the search of NG-STAR CCs.
```
pyngoST.py -i testNGSTAR.tsv -p /path/to/allelesDB -c -y 2
```

The example `testNGSTAR.tsv` file looks like:
```
strain	NG-STAR	penA	mtrR	porB	ponA	gyrA
31919_4_101.fasta	3658	15.001	424	100	100	100
31919_4_12.fasta	84	22.001	15	100	100	100
31919_4_16.fasta	168	2.001	39	11	100	7
31919_4_17.fasta	84	22.001	15	100	100	100
```

And the output table will have an appended column containing the NG-STAR CCs:
```
strain	NG-STAR	penA	mtrR	porB	ponA	gyrA	NG-STAR_CC
31919_4_101.fasta	3658	15.001	424	100	100	100	CC893
31919_4_12.fasta	84	22.001	15	100	100	100	CC352
31919_4_16.fasta	168	2.001	39	11	100	7	CC168
31919_4_17.fasta	84	22.001	15	100	100	100	CC352
```

In both `-y` and `-z` modes, the input file can also be comma-separated (CSV). The output results will show using the same separator than the input files.

Use `-o <outfilename>` to save output to file.

## Multithreading

`pyngoST` was conceived to be a fast command line tool for molecular typing of _N. gonorrhoeae_ genome collections. When using it on a large number of fasta files, some tasks can be computationally expensive, such as blasting new alleles to existing ones or, especially, calculating NG-MAST genogroups, which compromise speed. To overcome this, `pyngoST` includes multithreading using the `concurrent.futures` module.

Multithreading is applied to the `process_files()` function, which, for each input fasta file, finds exact matches of existing alleles for all loci of the three schemes using `pyahocorasick`, blasts new alleles if requested, and assigns STs from allele profiles.

Please note that only computationally intensive tasks are worth multithreading, such as combining the blasting of new alleles with the calculation of NG-MAST genogroups. You can multithread using `-t <number of threads`, i.e:
```
## Run pyngoST on 400 FASTA files requesting the blasting of new alleles with -b and genogroups calculation with -g on 8 threads
pyngoST.py -i /path/to/400genomes/*.fasta -p /path/to/allelesDB -s MLST,NG-MAST,NG-STAR -c -g -b -t 8
```
The request of NG-STAR CCs with `-c` is not a calculation but only a query of the existing database and it does not add significantly more computation time compared to not requesting CCs.
