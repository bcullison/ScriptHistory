import sys, os
from glob import glob
from Bio import SeqIO
from pathlib import Path
from Bio.Blast.NCBIWWW import qblast

###################This script takes the top 3 sequences from each file in the input folder based on best coverage and blasts them.###################

database = "/Users/katzlab_admin/Desktop/pr2dbblast/pr2dbblast" #insert path to your database here

def grab_top3_cov_rRNA(file):

	#this takes each file in the input folder and and creates a new file with the top 3 sequences based on highest coverage. 

	temp_list = []
	for thing in SeqIO.parse(file, 'fasta'):

		if '18S' in thing.id:
			name = file.split(".")[0].split("/")[-1]
			coverage = float(thing.id.split("_cov_")[1].split(":")[0])
			temp_list.append((thing, coverage))
	
	temp_list.sort(key = lambda x: -x[1])	
	
	rRNAnames = [f'>{name}_Germ18SrRNA_{thing[0].id.split("::")[1].split(":")[0]}\n{thing[0].seq}\n' for thing in temp_list[:3]]
	#SeqIO is expecting a fasta where each sequence starts with a >
	
	top3 = (f'{directory}_Top3')
	
	Path(top3).mkdir(exist_ok=True)
	
	with open(f'{top3}/{name}_18S_rRNA.fasta','w+') as final:
		final.write(''.join(rRNAnames))

	return top3

def rDNA_many_files(files):

	#this appends those top 3 sequences to a list in preparation for blasting

	top3_dirs = []

	for file in files:
	
		top3_dirs.append(grab_top3_cov_rRNA(file))
	
	return top3_dirs


def run_blastn(top3):

	#this blasts your sequences and creates an output file for your top 5 results

	print('\nBLAST-ing Sequences...\n')
	output = glob(f'{top3}/**fasta')
	Path("blast_results").mkdir(exist_ok=True)
	for outputs in output:
		newname = outputs.split("/")[-1].split(".fas")[0]
		print(newname)
		blast_cmd = f'blastn -query {outputs} -db {database} -max_target_seqs 5 -outfmt 6 -out blast_results/{newname}_pr2dbblast.tsv'
		os.system(blast_cmd)


if __name__ == '__main__':
	if len(sys.argv[1:]) != 1:
		print('Usage:\n    python3 barrnap_BLASTv2.py [DIRECTORY-WITH-BARRNAPd-FASTAs]\n')
		sys.exit(1)
	else:
		# directory = "Chilo_barrnap_output"
		directory = sys.argv[1]
	
	files = glob(f'{directory}/LKH**/**fasta')
	
	top3_dirs = rDNA_many_files(files)
	for top_dir in list(set(top3_dirs)):
		run_blastn(top_dir)

		





