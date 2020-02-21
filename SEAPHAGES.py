#!/usr/bin/env python3

import subprocess
import os
import glob
from Bio import SeqIO
import pandas as pd
import argparse

#import sys
#from Bio import SeqIO
#import pandas as pd

def my_args():
	# Create argument parser
	parser = argparse.ArgumentParser(description = "Consenus Genecalling of a FASTA file")

	# Positional mandatory arguments
	parser.add_argument("--infile", "-i", dest = "infile", required = True, help = "Define the input FASTA file in nuclotide format")
	parser.add_argument("--outfile", "-o", dest = "outfile", required = True, help = "Define the output .tsv file with consenus gene predictions")

	# Optional arguments
	parser.add_argument('--score', "-s", dest = "score", default = 3, type = int, help = "The minimum score a gene needs to be included in the final DF. Default = 3")

	# Parse arguments
	args = parser.parse_args()

	return args

class bcolours:
	HEADER = '\033[95m'
	OKBLUE = '\033[94m'
	OKGREEN = '\033[92m'
	WARNING = '\033[93m'
	FAIL = '\033[91m'
	ENDC = '\033[0m'
	BOLD = '\033[1m'
	UNDERLINE = '\033[4m'

def runGeneCallers(infile):
	def runProdigal(infile):
		'''
		Run the prodigal gene caller on the input FASTA file and output the gene calls in nuclotide format
		as well as the gene statistics. Use meta as prodigal requires over 100K length sequences to be accurate.
		Also combine the new prodigal gene names with the gene statistics for later analysis.
		'''
		Prodigal = f'prodigal -i {infile} -d Prodigal.fnn -o Prodigal.lst -f gff -p meta'
		stdout, stderr = execute(Prodigal)

	def runGlimmer(infile):
		'''
		Run the genecaller glimmer on the input FASTA file and output the gene calls in nuclotide format
		as well as gene statistics. Also combine the renamed Glimmer genes to the gene statistics for later analysis.
		Glimmer Flags Explained
		-----------------------
		long-orf : Find long, non-overlapping orfs to use as a training set.
		extract : Extract the training sequences from the genome file.
		build-icm : Build the icm from the training sequences.
		glimmer3 : Run the first instance of Glimmer.
		tail : Get training coordinates from first predictions.
		upstream-coords.awk : Create a position weight matrix (PWM) from the regions.
		elph : upstream of the start locations in run3.coords.
		startuse : Determine the distribution of start-codon usage in run3.coords.
		glimmer3 : Run second Glimmer.
		extract : Extract the nucleotide sequences from the genome file.
		'''

		Glimmer = f"""long-orfs -n -t 1.15 {infile} run3.longorfs;
			extract -t {infile} run3.longorfs > run3.train;
			build-icm -r run3.icm < run3.train;
			glimmer3 -o50 -g110 -t30 {infile} run3.icm run3.run1;
			tail -n +2 run3.run1.predict > run3.coords;
			upstream-coords.awk 25 0 run3.coords | extract {infile} - > run3.upstream;
			elph run3.upstream LEN=6 | get-motif-counts.awk > run3.motif;
			startuse="$(start-codon-distrib -3 {infile} run3.coords)";
			glimmer3 -o50 -g110 -t30 -b run3.motif -P $startuse {infile} run3.icm FINAL;
			extract -t {infile} FINAL.predict > Glimmer.fnn;
			mv FINAL.predict Glimmer.lst"""
		stdout, stderr = execute(Glimmer)

	def runMGA(infile):
		'''
		Run the gene caller MGA on the input FASTA file and output the genecalls in nuclotide format
		as well as the gene statistics. Use the `extract` function from glimmer to produce the nuclotide
		gene calls from the output file. Also combine the renamed MGA genes to the gene statistics for later analysis.
		'''
		MGA = f'mga {infile} -s > MGA.out; extract {infile} MGA.out > MGA.fnn'
		stdout, stderr = execute(MGA)


	def runGeneMarkS(infile):
		'''
		Run the gene caller GeneMarkS on the input FASTA file to create the .mat file for GeneMark and
		output the gene calls in nuclotide format as well as the gene caller statistics. Also combine the renamed
		GeneMarkS genes to the gene statistics for later analysis.
		'''
		GeneMarkS = f'gmsn.pl --gm --clean --fnn {infile} --format GFF; mv {infile}.fnn GeneMarkS.fnn'
		stdout, stderr = execute(GeneMarkS)

	def runGeneMark(infile):
		'''
		Run the gene caller GeneMark on the input FASTA file and the .mat file from GeneMarkS and
		output the gene calls in nuclotide format as well as the gene caller statistics. Also combine the renamed
		GeneMark genes to the gene statistics for later analysis
		'''
		GeneMark = f'gm -onq -loq -m GeneMark_heuristic.mat -D {infile}; mv {infile}.orf GeneMark.fnn'
		stdout, stderr = execute(GeneMark)

	def runGeneMarkS2(infile):
		'''
		Run the gene caller GeneMarkS2 on the input FASTA file and
		output the gene calls in nuclotide format as well as the gene caller statistics. Also combine the renamed
		GeneMarkS2 genes to the gene statistics for later analysis.
		'''
		GeneMarkS2 = f'gms2.pl --seq {infile} --genome-type auto --output GeneMarkS2.lst --fnn GeneMarkS2.fnn --format gff'
		stdout, stderr = execute(GeneMarkS2)

	print(f"{bcolours.HEADER}Running GeneCallers {bcolours.ENDC}")
	runProdigal(infile)
	print(f"{bcolours.OKGREEN}Prodigal Complete! {bcolours.ENDC}")
	runGlimmer(infile)
	print(f"{bcolours.OKGREEN}Glimmer Complete! {bcolours.ENDC}")
	runMGA(infile)
	print(f"{bcolours.OKGREEN}MGA Complete! {bcolours.ENDC}")
	runGeneMarkS(infile)
	print(f"{bcolours.OKGREEN}GeneMarkS Complete! {bcolours.ENDC}")
	runGeneMark(infile)
	print(f"{bcolours.OKGREEN}GeneMark Complete! {bcolours.ENDC}")
	runGeneMarkS2(infile)
	print(f"{bcolours.OKGREEN}GeneMarkS2 Complete! {bcolours.ENDC}")

def renameGeneCalls():
	'''
	Rename the outputted nuclotide genecalls into something sensible so they can be analysed in the future and combined with different genecalls
	that use the same name for the first protein. For example GeneMarkS and GeneMarkS2 have the same name for their first genecalled protein.
	'''
	file_list = glob.glob("*.fnn") # this needs to be renamed specifically for this files ["GeneMarkS2.fnn", "GeneMark.fnn", "GeneMarkS.fnn", "Glimmer.fnn", "Prodigal.fnn", "MGA.fnn"] located in the CWD.
	for f in file_list:
		results = []
		mappings = []
		counter = 1
		for record in SeqIO.parse(f, 'fasta'):
			new_name = os.path.basename(f).split(".fnn")[0] + f'_{counter:011d}'
			mappings.append((f, record.id, new_name))
			record.id = new_name
			record.name = ""
			record.description = ""
			results.append(record)
			counter +=1

		SeqIO.write(results, os.path.basename(f).split(".fnn")[0] + "_reformatted.fa", 'fasta')
		df = pd.DataFrame(mappings, columns=['filename', 'old_id', 'new_id'])
		df.to_csv(os.path.basename(f).split(".fnn")[0] + "_mapping.tsv", sep='\t', index=False)
	print(f"{bcolours.HEADER}Renamed the called genes into something more sensible {bcolours.ENDC}")

def getGeneStats():
	'''
	Combine the renamed genes with their gene statistics for comparison between gene callers
	'''
	print(f"{bcolours.HEADER}Collating the Gene Statistics into a list  {bcolours.ENDC}")
	Stats = f"""# Prodigal;
		grep ^[^#] Prodigal.lst | cut -f1,4,5,6,7 > Prodigal.trm.lst;
		grep ">" Prodigal_reformatted.fa | sed "s/>//g" | paste -d "\t" - Prodigal.trm.lst > Prodigal_genecalls.tsv;
		rm Prodigal.trm.lst Prodigal.lst;
		# Glimmer;
		grep "^[^>]" Glimmer.lst | sed "s/-[1-3]/-/g; s/+[1-3]/+/g; s/ \+ /\t/g" > Glimmer.trm.lst # change strand number to + or - for standardisation vs other gene callers;
		grep ">" Glimmer_reformatted.fa | sed "s/>//g" | paste -d "\t" - Glimmer.trm.lst > Glimmer_genecalls.tsv;
		rm FINAL.detail run3.* Glimmer.lst Glimmer.trm.lst;
		# MGA;
		grep "^[^#]" MGA.out | cut -f1,2,3,4,7 > MGA.trm.lst
		grep ">" MGA_reformatted.fa | sed "s/>//g" | paste -d "\t" - MGA.trm.lst > MGA_genecalls.tsv
		rm MGA.out MGA.trm.lst;
		# GeneMarkS
		grep "^[^#]" HTVC010P_reformatted.fasta.gff | cut -f4,5,6,7,9 > GeneMarkS.trm.lst;
		grep ">" GeneMarkS_reformatted.fa | sed "s/>//g" | paste -d "\t" - GeneMarkS.trm.lst > GeneMarkS_genecalls.tsv;
		rm GeneMark_hmm_heuristic.mod GeneMarkS.trm.lst gms.log HTVC010P_reformatted.fasta.gff;
		# GeneMark
		grep ">" GeneMark.fnn | grep -o ',.*$' | sed 's/, //g; s/ - /\t/g; s/)//g' > GeneMark_sstart_send.tsv;
		grep ">" GeneMark_reformatted.fa | sed "s/>//g" | paste -d "\t" - GeneMark_sstart_send.tsv > GeneMark_genecalls.tsv;
		rm GeneMark_heuristic.mat GeneMark_sstart_send.tsv HTVC010P_reformatted.fasta.gdata HTVC010P_reformatted.fasta.ldata HTVC010P_reformatted.fasta.lst;
		# GeneMarkS2
		grep "^[^#]" GeneMarkS2.lst | cut -f4,5,6,7 > GeneMarkS2.trm.lst;
		grep ">" GeneMarkS2_reformatted.fa | sed "s/>//g" | paste -d "\t" - GeneMarkS2.trm.lst > GeneMarkS2_genecalls.tsv;
		rm GMS2.mod GeneMarkS2.trm.lst GeneMarkS2.lst log"""

	stdout, stderr = execute(Stats)

def runCD_hit():
	'''
	Run cd-hit, a protein clustering software to determine if genecalled proteins from different protein softwares are the infact the same gene.
	Clusters are defined as 0.95 percentage identity over 90% of the gene with word search size of 10 and take description of gene to first defline(a space)
	'''
	print(f"{bcolours.HEADER}Clustering simular genes into groups to create consenus gene calling{bcolours.ENDC}")
	Cluster = f"""cat *_reformatted.fa > Genecallers_combined.fa;
		cd-hit-est -i Genecallers_combined.fa -o cd-hit.out -c 0.95 -n 10 -s 0.9 -d 0;
		clstr2txt.pl cd-hit.out.clstr > Gene_clstr.txt;
		rm cd-hit.out Genecallers_combined.fa cd-hit.out.clstr"""

	stdout, stderr = execute(Cluster)

def combineGeneClusters():
	'''
	1) Import all the genecaller statistics and merge with the cd-hit protein clustering results.
	2) Normalise the genescores to allow for comparisons against other genecallers genescores.
		a) This is done by dividing all the genescores from that genecaller against the highest score that specific genecaller scored in this dataset.
		Therefore this number is based on the genecaller not the merged DF as a whole
	3) Calculate gene length. This allows for standardisation of how Gene Length is defined as Glimmer's is different (3bps shorter then GeneMark Family for instance).
	4) Create a penality system for short genes.
		Scoring
		-------
		If Gene length is between 200 and 150 bps, -1 point
		If Gene Length is between 150 and 120 bps, -2 points
		If Gene Length is between 120 and 90 bps, -3 points
		If Gene Length is between 90 and 75 bps, -4 points
		If Gene Length is smaller than 75 bps, -inf points (as no gene as been recorded to be smaller than 27 amino acids long)
	5) Therefore the scoring system for each cluster of genes determined by the Normalised Gene Score, Cluster size and length penality.
	The representative of each cluster is chosen based on the highest score.
	If there is a tie, the first gene in the cluster is chosen. (First = position listed in the dataframe).
	'''
	print(f"{bcolours.HEADER}Collating all genecalled genes from different algorithims together{bcolours.ENDC}")
	#import DFs
	GeneMarkS2 = pd.read_csv('GeneMarkS2_genecalls.tsv', sep = "\t", names = ["id", "Sstart", "Send", "Gene_score", "Direction"], index_col=False)
	GeneMark = pd.read_csv('GeneMark_genecalls.tsv', sep = "\t", names = ["id", "Sstart", "Send"], index_col=False)
	GeneMarkS = pd.read_csv('GeneMarkS_genecalls.tsv', sep = "\t", names = ["id", "Sstart", "Send", "Gene_score", "Direction", "old_name"], index_col=False)
	Glimmer = pd.read_csv('Glimmer_genecalls.tsv', sep = "\t", names = ["id", "old_name", "Sstart", "Send", "Direction", "Gene_score"], index_col=False)
	MGA = pd.read_csv('MGA_genecalls.tsv', sep = "\t", names = ["id", "old_name", "Sstart", "Send", "Direction", "Gene_score"], index_col=False)
	Prodigal = pd.read_csv('Prodigal_genecalls.tsv', sep = "\t", names = ["id", "Contig", "Sstart", "Send", "Gene_score","Direction"], index_col=False)
	df = pd.read_csv("Gene_clstr.txt", sep = "\t")

	# normalise the genescores
	print(f"{bcolours.HEADER}Normalising Gene Scores based on highest value{bcolours.ENDC}")
	GeneMarkS2["normalised_gene_score"] = GeneMarkS2["Gene_score"]/max(GeneMarkS2["Gene_score"])
	GeneMarkS["normalised_gene_score"] = GeneMarkS["Gene_score"]/max(GeneMarkS["Gene_score"])
	Glimmer["normalised_gene_score"] = Glimmer["Gene_score"]/max(Glimmer["Gene_score"])
	MGA["normalised_gene_score"] = MGA["Gene_score"]/max(MGA["Gene_score"])
	Prodigal["normalised_gene_score"] = Prodigal["Gene_score"]/max(Prodigal["Gene_score"])

	# merge all the DFs together
	final = GeneMarkS2.merge(GeneMark, how="outer").merge(GeneMarkS, how='outer').merge(Glimmer, how='outer').merge(MGA, how='outer').merge(Prodigal, how='outer').merge(df, on="id")
	final["Contig"] = 0

	print(f"{bcolours.OKGREEN}Calculating Gene Length{bcolours.ENDC}")
	# Calculate gene length because Glimmer uses a different system (3bp short)
	final["length"] = final[["Sstart", "Send"]].max(axis=1) - final[["Sstart", "Send"]].min(axis=1)

	# Create length Penality
	print(f"{bcolours.OKGREEN}Creating a Penality value for Short Genes{bcolours.ENDC}")
	score_bins = [0, 74, 89, 119, 149, 199, 99999]
	penality = [-999, -4, -3, -2, -1, 0]
	final["length_penality"] = pd.cut(final["length"], score_bins, labels=penality)

	# Create scoring system
	print(f"{bcolours.OKGREEN}Creating an initial score of each gene to select the best gene within a cluster{bcolours.ENDC}")
	final["score"] = final[["normalised_gene_score", "clstr_size", "length_penality"]].sum(axis = 1) # + final.e_value # create score
	idx = final.groupby(by=final["clstr"])["score"].transform(max)==final["score"] # grp by cluster keep clstr representative by highest score
	result_df = final[idx].drop_duplicates(subset=['clstr'], keep='first') # remove duplicate clstr representative if same score

	return result_df

def GeneOverlap(result_df, outfile, score):
	'''
	Calculate an overlap penality for each gene. However, only penalise genes with a lower gene cluster number.
	For example, a Gene that was found with all 6 genecallers should not be penalised if it overlaps with a gene that was only found with 1 genecaller.
	This is done by calculating the range between the start and end of each gene and comparing it against every other gene.
	If an overlap is recorded, not the size of the overlap and refering to the scoring system below to penalise that gene.
	However, if the overlap is either 1,4 or 8 bps in size, this is an operon. Some genes work as operons and overlap in this way.
	If an overlap of 1,4 or 8 is recored, instead award +1 point, as operons and ubiqutuis within viral and bacterial genomes.
	Scoring System
	--------------
	If overlap is 1, 4 or 8 bps, +1 point
	If overlap is not 1, 4 or 8 bps, then:
	If overlap is between 10 and 40 bps, -1 point
	If overlap is between 40 and 70 bps, -2 points
	If overlap is between 70 and 100 bps, -3 points
	If overlap is larger then 100 bps, -inf points. (This is because genes that overlap by more that 100 bps have never been observed in nature).
	This does not take into account if genes starting in the same location on the forward strand must have a 50bp gap between a gene starting in the same location on the reverse strand.
	This is to account for the space needed for transcription factors to be encoded.
	'''
	# Create dfs by cluster size aka consenus
	print(f"{bcolours.OKGREEN}Creating a penality system for genes that overlap{bcolours.ENDC}")
	agree_6 = result_df[result_df["clstr_size"] == 6].reset_index(drop=True)
	agree_5 = result_df[result_df["clstr_size"] == 5].reset_index(drop=True)
	agree_4 = result_df[result_df["clstr_size"] == 4].reset_index(drop=True)
	agree_3 = result_df[result_df["clstr_size"] == 3].reset_index(drop=True)
	agree_2 = result_df[result_df["clstr_size"] == 2].reset_index(drop=True)
	agree_1 = result_df[result_df["clstr_size"] == 1].reset_index(drop=True)

	print(f"{bcolours.OKGREEN}Checking to see if genes with the highest consenus overlap{bcolours.ENDC}")
	list_of_datasets = [agree_6, agree_5, agree_4, agree_3, agree_2, agree_1] # create a list of dataframes to loop through

	# Create
	for index, dataset in enumerate(list_of_datasets): # for each dataset compare against other datasets
		overlap = [] # set overlap counter to 0
		operon = [] # set operon counter to 0
		dataset["operon"] = 0 # create a column in each DF that records operon bonus overlap, resetting it if it already exists
		dataset["overlap"] = 0 # create a column in each DF that records overlap penalitye, resetting it if it already exists
		for i, row in agree_6.iterrows(): # for every row in the consensus 6 Df
			for j, row in dataset.iterrows(): # compare it to each row in the other 6 dataframes
				idx1 = pd.Index(range(agree_6[["Sstart", "Send"]].min(axis=1)[i], agree_6[["Sstart", "Send"]].max(axis=1)[i])) # calculate the first gene range in the consenus 6 DF
				idx2 = pd.Index(range(dataset[["Sstart", "Send"]].min(axis=1)[j], dataset[["Sstart", "Send"]].max(axis=1)[j])) # compare that to the first gene of the next dataframe
				if len(idx1.intersection(idx2)) == 1 or len(idx1.intersection(idx2)) == 4 or len(idx1.intersection(idx2)) == 8: # if the overlap of the genes are 1,4 or 8, this is an operon
					operon.append(1) # so count it as operon
					dataset["operon"] += operon # and append it to the DF
				else:
					overlap.append(len(idx1.intersection(idx2))) # otherwise if it is any other score, calcualte the length of the overlap
			dataset["overlap"] += overlap # and add it to the DF column as a score
			overlap = [] # reset the overlap counter and repeat this for the next DF (in this case the next DF is consenus 5)
	agree_6["overlap"] -= agree_6["length"] # because you are comparing this with the same gene you need deduct the overlap for measuring the same genes overlap

	print(f"{bcolours.OKGREEN}Checking to see if genes with the 2nd highest consenus overlap{bcolours.ENDC}")
	list_of_datasets = [agree_5, agree_4, agree_3, agree_2, agree_1]# repeat this process comparing consenus 5's genes with all other genes except consenus 6 which has a higest consenus score

	for index, dataset in enumerate(list_of_datasets):
		overlap = []
		operon = []
		dataset["operon"] += 0 # create a column in each DF that records operon penality, but if it already exist keep the score
		dataset["overlap"] += 0 # create a column in each DF that records overlap penality, but if it already exist keep the score
		for i, row in agree_5.iterrows():
			for j, row in dataset.iterrows():
				idx1 = pd.Index(range(agree_5[["Sstart", "Send"]].min(axis=1)[i], agree_5[["Sstart", "Send"]].max(axis=1)[i]))
				idx2 = pd.Index(range(dataset[["Sstart", "Send"]].min(axis=1)[j], dataset[["Sstart", "Send"]].max(axis=1)[j]))
				if len(idx1.intersection(idx2)) == 1 or len(idx1.intersection(idx2)) == 4 or len(idx1.intersection(idx2)) == 8:
					operon.append(1)
					dataset["operon"] += operon
				else:
					overlap.append(len(idx1.intersection(idx2)))
			dataset["overlap"] += overlap
			overlap = []
	agree_5["overlap"] -= agree_5["length"]

	print(f"{bcolours.OKGREEN}Checking to see if genes with the 3rd highest consenus overlap{bcolours.ENDC}")
	list_of_datasets = [agree_4, agree_3, agree_2, agree_1]

	for index, dataset in enumerate(list_of_datasets):
		overlap = []
		operon = []
		dataset["operon"] += 0
		dataset["overlap"] += 0
		for i, row in agree_4.iterrows():
			for j, row in dataset.iterrows():
				idx1 = pd.Index(range(agree_4[["Sstart", "Send"]].min(axis=1)[i], agree_4[["Sstart", "Send"]].max(axis=1)[i]))
				idx2 = pd.Index(range(dataset[["Sstart", "Send"]].min(axis=1)[j], dataset[["Sstart", "Send"]].max(axis=1)[j]))
				if len(idx1.intersection(idx2)) == 1 or len(idx1.intersection(idx2)) == 4 or len(idx1.intersection(idx2)) == 8:
					operon.append(1)
					dataset["operon"] += operon
				else:
					overlap.append(len(idx1.intersection(idx2)))
			dataset["overlap"] += overlap
			overlap = []
	agree_4["overlap"] -= agree_4["length"]

	print(f"{bcolours.OKGREEN}Checking to see if genes with the 4th highest consenus overlap{bcolours.ENDC}")
	list_of_datasets = [agree_3, agree_2, agree_1]

	for index, dataset in enumerate(list_of_datasets):
		overlap = []
		operon = []
		dataset["operon"] += 0
		dataset["overlap"] += 0
		for i, row in agree_3.iterrows():
			for j, row in dataset.iterrows():
				idx1 = pd.Index(range(agree_3[["Sstart", "Send"]].min(axis=1)[i], agree_3[["Sstart", "Send"]].max(axis=1)[i]))
				idx2 = pd.Index(range(dataset[["Sstart", "Send"]].min(axis=1)[j], dataset[["Sstart", "Send"]].max(axis=1)[j]))
				if len(idx1.intersection(idx2)) == 1 or len(idx1.intersection(idx2)) == 4 or len(idx1.intersection(idx2)) == 8:
					operon.append(1)
					dataset["operon"] += operon
				else:
					overlap.append(len(idx1.intersection(idx2)))
			dataset["overlap"] += overlap
			overlap = []
	agree_3["overlap"] -= agree_3["length"]

	print(f"{bcolours.OKGREEN}Checking to see if genes with the 5th highest consenus overlap{bcolours.ENDC}")
	list_of_datasets = [agree_2, agree_1]

	for index, dataset in enumerate(list_of_datasets):
		overlap = []
		operon = []
		dataset["operon"] += 0
		dataset["overlap"] += 0
		for i, row in agree_2.iterrows():
			for j, row in dataset.iterrows():
				idx1 = pd.Index(range(agree_2[["Sstart", "Send"]].min(axis=1)[i], agree_2[["Sstart", "Send"]].max(axis=1)[i]))
				idx2 = pd.Index(range(dataset[["Sstart", "Send"]].min(axis=1)[j], dataset[["Sstart", "Send"]].max(axis=1)[j]))
				if len(idx1.intersection(idx2)) == 1 or len(idx1.intersection(idx2)) == 4 or len(idx1.intersection(idx2)) == 8:
					operon.append(1)
					dataset["operon"] += operon
				else:
					overlap.append(len(idx1.intersection(idx2)))
			dataset["overlap"] += overlap
			overlap = []
	agree_2["overlap"] -= agree_2["length"]

	print(f"{bcolours.OKGREEN}Checking to see if genes with the lowest consenus overlap with themselves{bcolours.ENDC}")
	list_of_datasets = [agree_1]

	for index, dataset in enumerate(list_of_datasets):
		overlap = []
		operon = []
		dataset["operon"] += 0
		dataset["overlap"] += 0
		for i, row in agree_1.iterrows():
			for j, row in dataset.iterrows():
				idx1 = pd.Index(range(agree_1[["Sstart", "Send"]].min(axis=1)[i], agree_1[["Sstart", "Send"]].max(axis=1)[i]))
				idx2 = pd.Index(range(dataset[["Sstart", "Send"]].min(axis=1)[j], dataset[["Sstart", "Send"]].max(axis=1)[j]))
				if len(idx1.intersection(idx2)) == 1 or len(idx1.intersection(idx2)) == 4 or len(idx1.intersection(idx2)) == 8:
					operon.append(1)
					dataset["operon"] += operon
				else:
					overlap.append(len(idx1.intersection(idx2)))
			dataset["overlap"] += overlap
			overlap = []
	agree_1["overlap"] -= agree_1["length"]

	# Merge the dataframes back together
	results_df = agree_6.merge(agree_5, how='outer').merge(agree_4, how='outer').merge(agree_3, how='outer').merge(agree_2, how='outer').merge(agree_1, how='outer')

	# Create a an overlap scoring bin
	print(f"{bcolours.OKGREEN}Penalising genes that overlap{bcolours.ENDC}")
	score_bins = [-99999, 10, 40, 70, 100, 99999]
	penality = [0, -1, -2, -3, -99999]
	results_df["overlap_pen"] = pd.cut(results_df["overlap"], score_bins, labels=penality)

	# Create a total scoring system based on the previous scoring system plus the operon and overlap penality
	print(f"{bcolours.OKGREEN}Creating the final scoring system for representative genes of each cluster{bcolours.ENDC}")
	results_df["total_score"] = results_df[["overlap_pen", "operon", "score"]].sum(axis = 1)

	# Only include the gene if the score is greater than 3
	print(f"{bcolours.OKGREEN}Only keeping genes with a score higher then {score}{bcolours.ENDC}")
	final_df = results_df[results_df["total_score"] >= score]

	# Sort the dataframe by gene start
	sorted_df = final_df.sort_values(by=["Sstart"]).reset_index(drop=True)

	# Save the sorted_df to file
	print(f"{bcolours.OKGREEN}Gene statistics can be found here: {outfile}{bcolours.ENDC}")
	sorted_df.to_csv(f"{outfile}", sep = "\t", index = None, header = True)

	return sorted_df

def execute(bash):
	process = subprocess.Popen(bash, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	(stdout, stderr) = process.communicate()

#	print(stdout.decode("utf-8"), stderr.decode("utf-8")) # this should go to a log file instead
	return stdout, stderr


def main(args):
	runGeneCallers(args.infile)
	renameGeneCalls()
	getGeneStats()
	runCD_hit()
	GeneOverlap(combineGeneClusters(), args.outfile, args.score)

if __name__=="__main__":
	args = my_args()
	main(args)
