#!/usr/bin/env python3

import subprocess
import os
import glob
from Bio import SeqIO
import pandas as pd
import argparse
from functools import reduce
from itertools import groupby
import python_colours as col

def my_args():
    # Create argument parser
    parser = argparse.ArgumentParser(description = "Consenus Genecalling of a FASTA file")

    # Positional mandatory arguments
    parser.add_argument("--infile", "-i", dest = "infile", required = True, help = "Define the input FASTA file in nuclotide format")
    parser.add_argument("--prefix", "-p", dest = "prefix", required = True, help = "Define the output prefix for the consenus gene predictions")

    # Optional arguments
    parser.add_argument('--score', "-s", dest = "score", default = 3, type = int, help = "The minimum score a gene needs to be included in the final DF. Default = 3")

    # Parse arguments
    args = parser.parse_args()

    return args

def runGeneCallers(infile):
    print(f"{col.colours.BBlue}Running GeneCallers {col.colours.Colour_Off}")
    def SplitByContig(infile):
        '''
        File is split by contig to ensure overlap function doesn't mistake different contigs for the same genes.
        This also improves accuracy of cluster formation as repeat genes tend to form the same cluster.
        Also multithreading can be implimented in the future.
        '''
        print(f"{col.colours.BGreen}Splitting FASTA by contig {col.colours.Colour_Off}")

        contigsplit = f"""awk '{{if (substr($0, 1, 1)==">") {{filename=(substr($0,2) "_testing.fa")}} print $0 >filename}}' {infile}"""
        stdout, stderr = execute(contigsplit)

    def runProdigal():
        '''
        Run the prodigal gene caller on the input FASTA file and output the gene calls in nuclotide format
        as well as the gene statistics. Use meta as prodigal requires over 100K length sequences to be accurate.
        '''
        print(f"{col.colours.BGreen}Running Prodigal{col.colours.Colour_Off}")

        Prodigal = f"""for i in *_testing.fa; do prodigal -i $i -d $(basename $i .fa)_Prodigal.fnn -o $(basename $i .fa)_Prodigal.lst -f gff -p meta; done"""
        stdout, stderr = execute(Prodigal)

    def runGlimmer():
        '''
        Run the genecaller glimmer on the input FASTA file and output the gene calls in nuclotide format
        as well as gene statistics.
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
        print(f"{col.colours.BGreen}Running Glimmer{col.colours.Colour_Off}")

        Glimmer = f"""for i in *_testing.fa; do \
            long-orfs -n -t 1.15 $i run3.longorfs;
            extract -t $i run3.longorfs > run3.train;
            build-icm -r run3.icm < run3.train;
            glimmer3 -o50 -g110 -t30 $i run3.icm run3.run1;
            tail -n +2 run3.run1.predict > run3.coords;
            upstream-coords.awk 25 0 run3.coords | extract $i - > run3.upstream;
            elph run3.upstream LEN=6 | get-motif-counts.awk > run3.motif;
            startuse="$(start-codon-distrib -3 $i run3.coords)";
            glimmer3 -o50 -g110 -t30 -b run3.motif -P $startuse $i run3.icm FINAL;
            extract -t $i FINAL.predict > $(basename $i .fa)_Glimmer.fnn;
            mv FINAL.predict $(basename $i .fa)_Glimmer.lst;
            done"""
        stdout, stderr = execute(Glimmer)

    def runMGA():
        '''
        Run the gene caller MGA on the input FASTA file and output the genecalls in nuclotide format
        as well as the gene statistics. Use the `extract` function from glimmer to produce the nuclotide
        gene calls from the output file.
        '''
        print(f"{col.colours.BGreen}Running MGA{col.colours.Colour_Off}")

        MGA = f'for i in *_testing.fa; do mga $i -s > $(basename $i .fa)_MGA.lst; extract $i $(basename $i .fa)_MGA.lst > $(basename $i .fa)_MGA.fnn; done'
        stdout, stderr = execute(MGA)


    def runGeneMarkS():
        '''
        Run the gene caller GeneMarkS on the input FASTA file to create the .mat file for GeneMark and
        output the gene calls in nuclotide format as well as the gene caller statistics.
        '''

        print(f"{col.colours.BGreen}Running GeneMarkS{col.colours.Colour_Off}")

        GeneMarkS = f"""for i in *_testing.fa; do gmsn.pl --gm --clean --fnn $i --format GFF;
            mv $i.fnn $(basename $i .fa)_GeneMarkS.fnn;
            mv $i.gff $(basename $i .fa)_GeneMarkS.lst;
            done"""
        stdout, stderr = execute(GeneMarkS)

    #def runGeneMark(infile):
    #    '''
    #    Run the gene caller GeneMark on the input FASTA file and the .mat file from GeneMarkS and
    #    output the gene calls in nuclotide format as well as the gene caller statistics. Also combine the renamed
    #    GeneMark genes to the gene statistics for later analysis
    #    '''
    #
    #    print(f"{col.colours.BGreen}GeneMark Complete! {col.colours.Colour_Off}")
    #
    #    GeneMark = f'gm -onq -loq -m *.mat -D {infile}; mv {infile}.orf GeneMark.fnn'
    #    stdout, stderr = execute(GeneMark)

    def runGeneMarkS2():
        '''
        Run the gene caller GeneMarkS2 on the input FASTA file and
        output the gene calls in nuclotide format as well as the gene caller statistics. Also combine the renamed
        '''

        print(f"{col.colours.BGreen}GeneMarkS2 Complete! {col.colours.Colour_Off}")

        GeneMarkS2 = f"""for i in *_testing.fa; do \
            gms2.pl --seq $i --genome-type auto --output $(basename $i .fa)_GeneMarkS2.lst --fnn $(basename $i .fa)_GeneMarkS2.fnn --format gff;
            done"""
        stdout, stderr = execute(GeneMarkS2)

    def runPHANOTATE():
        '''
        '''

        print(f"{col.colours.BGreen}PHANOTATE Complete! {col.colours.Colour_Off}")

        PHANOTATE = f"""for i in *_testing.fa; do \
            phanotate.py $i -o $(basename $i .fa)_PHANOTATE.out;
            awk '!/#/ {{print "gene_"NR-2"\t", $0}}' $(basename $i .fa)_PHANOTATE.out > $(basename $i .fa)_PHANOTATE.lst;
            extract $i $(basename $i .fa)_PHANOTATE.lst > $(basename $i .fa)_PHANOTATE.fnn;
            done"""
        stdout, stderr = execute(PHANOTATE)

    SplitByContig(infile)
    runProdigal()
    runGlimmer()
    runMGA()
    runGeneMarkS()
    #runGeneMark(infile)
    runGeneMarkS2()
    runPHANOTATE()

def renameGeneCalls():
    '''
    Rename the outputted nuclotide genecalls into something sensible so they can be analysed in the future and combined with different genecalls
    that use the same name for the first protein. For example GeneMarkS and GeneMarkS2 have the same name for their first genecalled protein.
    '''
    print(f"{col.colours.BPurple}Renaming the called genes into something more sensible {col.colours.Colour_Off}")

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

def getGeneStats(infile):
    '''
    Combine the renamed genes with their gene statistics for comparison between gene callers
    '''
    print(f"{col.colours.BPurple}Collating the Gene Statistics into a list that is parseable  {col.colours.Colour_Off}")
    Stats = f"""
        # Prodigal;
        for i in *_testing_Prodigal.lst; do \
        grep ^[^#] $i | cut -f1,4,5,6,7 > $(basename $i .lst).trm.lst;
        grep ">" $(basename $i .lst)_reformatted.fa | sed "s/>//g" | paste -d "\t" - $(basename $i .lst).trm.lst > $(basename $i .lst)_genecalls.tsv;
        rm $(basename $i .lst).trm.lst $(basename $i .lst).lst $(basename $i .lst).fnn $(basename $i .lst)_mapping.tsv;
        done;
        # Glimmer;
        for i in *_testing_Glimmer.lst; do \
        awk '/>/{{filename=NR"_glimmer.txt"}}; {{print >filename}}' $i;
        for j in *_glimmer.txt; do sed "s/orf/$(grep ">" $j)\torf/g" $j | tail -n +2 | sed 's/>//g' >> $(basename $i .lst).lst2; done
        sed "s/-[1-3]/-/g; s/+[1-3]/+/g; s/ \+ /\t/g" $(basename $i .lst).lst2 > $(basename $i .lst).trm.lst # change strand number to + or - for standardisation vs other gene callers;
        grep ">" $(basename $i .lst)_reformatted.fa | sed "s/>//g" | paste -d "\t" - $(basename $i .lst).trm.lst > $(basename $i .lst)_genecalls.tsv;
        rm FINAL.detail run3.* $(basename $i .lst).lst $(basename $i .lst).trm.lst *_glimmer.txt $(basename $i .lst).lst2 $(basename $i .lst)_mapping.tsv $(basename $i .lst).fnn;
        done;
        # MGA;
        for i in *_testing_MGA.lst; do \
        egrep -v "# gc|# self" $i | awk '/# /{{filename=NR"_mga.txt"}}; {{print >filename}}';
        for j in *_mga.txt; do sed "s/gene_/$(grep "# " $j)\tgene_/g" $j | tail -n +2 | sed 's/# //g' >> $(basename $i .lst).lst2; done;
        cut -f1,2,3,4,5,8,12 $(basename $i .lst).lst2 > $(basename $i .lst).trm.lst;
        grep ">" $(basename $i .lst)_reformatted.fa | sed "s/>//g" | paste -d "\t" - $(basename $i .lst).trm.lst > $(basename $i .lst)_genecalls.tsv;
        rm $i $(basename $i .lst).trm.lst *_mga.txt $(basename $i .lst).lst2 $(basename $i .lst).fnn $(basename $i .lst)_mapping.tsv;
        done;
        # GeneMarkS
        for i in *_testing_GeneMarkS.lst; do \
        grep "^[^#]" $i | cut -f1,4,5,6,7 > $(basename $i .lst).trm.lst;
        grep ">" $(basename $i .lst)_reformatted.fa | sed "s/>//g" | paste -d "\t" - $(basename $i .lst).trm.lst > $(basename $i .lst)_genecalls.tsv;
        rm GeneMark_hmm_heuristic.mod GeneMark_hmm.mod $(basename $i .lst).trm.lst gms.log $i $(basename $i .lst).fnn $(basename $i .lst)_mapping.tsv;
        done;
        # GeneMark
        # grep ">" GeneMark.fnn | grep -o ',.*$' | sed 's/, //g; s/ - /\t/g; s/)//g' > GeneMark_sstart_send.tsv;
        # grep ">" GeneMark_reformatted.fa | sed "s/>//g" | paste -d "\t" - GeneMark_sstart_send.tsv > GeneMark_genecalls.tsv;
        rm GeneMark_heuristic.mat GeneMark.mat;
        # rm GeneMark_sstart_send.tsv {infile}.gdata {infile}.ldata {infile}.lst;
        # GeneMarkS2
        for i in *_testing_GeneMarkS2.lst; do \
        grep "^[^#]" $i | cut -f1,4,5,6,7 > $(basename $i .lst).trm.lst;
        grep ">" $(basename $i .lst)_reformatted.fa | sed "s/>//g" | paste -d "\t" - $(basename $i .lst).trm.lst > $(basename $i .lst)_genecalls.tsv;
        rm GMS2.mod $(basename $i .lst).trm.lst $(basename $i .lst).lst log $(basename $i .lst).fnn $(basename $i .lst)_mapping.tsv;
        done;
        # PHANOTATE;
        for i in *_testing_PHANOTATE.lst; do \
        grep ">" $(basename $i .lst)_reformatted.fa | sed "s/>//g" | paste -d "\t" - $(basename $i .lst).lst | cut -f1,3,4,5,6,7 > $(basename $i .lst)_genecalls.tsv;
        rm $(basename $i .lst).lst $(basename $i .lst).out $(basename $i .lst).fnn $(basename $i .lst)_mapping.tsv;
        done"""

    stdout, stderr = execute(Stats)

def runCD_hit(infile):
    '''
    Run cd-hit, a protein clustering software to determine if genecalled proteins from different protein softwares are the infact the same gene.
    Clusters are defined as 0.95 percentage identity over 90% of the gene with word search size of 10 and take description of gene to first defline(a space)
    '''
    print(f"""{col.colours.BPurple}Clustering simular genes into groups to create consenus gene calling.
    Cluster formation required 97% identity accross 97% of the gene.
    Higher numbers are better to prevent simular gene domains that are actually differnt genes from clustering{col.colours.Colour_Off}""")

    Cluster = f"""for i in $(grep ">" {infile} | sed "s/>//g"); do cat ${{i}}_testing*reformatted.fa > ${{i}}_Genecallers_combined.fa; done;
        for i in *_Genecallers_combined.fa; do \
        cd-hit-est -i $i -o $(basename $i .fa).out -c 0.97 -n 7 -s 0.97 -g 1 -d 0;
        clstr2txt.pl $(basename $i .fa).out.clstr > $(basename $i .fa)_clstr.txt;
        rm $(basename $i .fa).out $(basename $i .fa).out.clstr $(basename $i _Genecallers_combined.fa)_testing*reformatted.fa;
        done"""

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
    print(f"{col.colours.BPurple}Collating all genecalled genes from different algorithims together{col.colours.Colour_Off}")
    #import DFs
    file_list = glob.glob('*testing*')
    file_list.sort()

    clusters = glob.glob('*_Genecallers_combined_clstr.txt')
    clusters.sort()

    contig_list = [list(i) for j, i in groupby(file_list, lambda a: a.split('_')[0])]

    for i in range(len(contig_list)):
        GeneMarkS2 = pd.read_csv(contig_list[i][1], sep = "\t", names = ["id", "contig", "Sstart", "Send", "Gene_score", "Direction"], index_col=False)
        #GeneMark = pd.read_csv('GeneMark_genecalls.tsv', sep = "\t", names = ["id", "Sstart", "Send"], index_col=False)
        GeneMarkS = pd.read_csv(contig_list[i][2], sep = "\t", names = ["id", "contig","Sstart", "Send", "Gene_score", "Direction"], index_col=False)
        Glimmer = pd.read_csv(contig_list[i][3], sep = "\t", names = ["id", "contig", "old_name", "Sstart", "Send", "Direction", "Gene_score"], index_col=False)
        MGA = pd.read_csv(contig_list[i][4], sep = "\t", names = ["id", "contig", "old_name", "Sstart", "Send", "Direction", "Gene_score", "rbs_score"], index_col=False)
        PHANOTATE = pd.read_csv(contig_list[i][5], sep = "\t", names = ["id", "Sstart", "Send", "Direction", "contig", "Gene_score"], index_col=False)
        Prodigal = pd.read_csv(contig_list[i][6], sep = "\t", names = ["id", "contig", "Sstart", "Send", "Gene_score","Direction"], index_col=False)
        df = pd.read_csv(clusters[i], sep = "\t")

        # normalise the genescores
        print(f"{col.colours.BPurple}Normalising Gene Scores based on highest value{col.colours.Colour_Off}")
        GeneMarkS2["normalised_gene_score"] = GeneMarkS2["Gene_score"]/max(GeneMarkS2["Gene_score"])
        GeneMarkS["normalised_gene_score"] = GeneMarkS["Gene_score"]/max(GeneMarkS["Gene_score"])
        Glimmer["normalised_gene_score"] = Glimmer["Gene_score"]/max(Glimmer["Gene_score"])
        MGA["normalised_gene_score"] = MGA["Gene_score"]/max(MGA["Gene_score"])
        Prodigal["normalised_gene_score"] = Prodigal["Gene_score"]/max(Prodigal["Gene_score"])
        PHANOTATE["normalised_gene_score"] = PHANOTATE["Gene_score"]/min(PHANOTATE["Gene_score"])

        # merge all the DFs together
        dataframe=[GeneMarkS2, GeneMarkS, Glimmer, MGA, Prodigal, PHANOTATE]
        test=reduce(lambda x, y: pd.merge(x, y, how='outer'), dataframe)#.merge(GeneMark, how = "outer")
        final = pd.merge(df, test, on ='id')
        final["rbs_score"] = 0 # 0 because those without an rbs_score will mess everything up

        print(f"{col.colours.BGreen}Calculating Gene Length{col.colours.Colour_Off}")
        # Calculate gene length because Glimmer uses a different system (3bp short)
        final["length"] = final[["Sstart", "Send"]].max(axis=1) - final[["Sstart", "Send"]].min(axis=1)

        # Create length Penality
        print(f"{col.colours.BGreen}Creating a Penality value for Short Genes{col.colours.Colour_Off}")
        score_bins = [0, 74, 89, 119, 149, 199, 99999]
        penality = [-9, -4, -3, -2, -1, 0]
        final["length_penality"] = pd.cut(final["length"], score_bins, labels=penality)

        # Create scoring system
        print(f"{col.colours.BGreen}Creating an initial score of each gene to select the best gene within a cluster{col.colours.Colour_Off}")
        final["score"] = final[["normalised_gene_score", "clstr_size", "length_penality"]].sum(axis = 1) # + final.e_value # create score
        idx = final.groupby(by=final["clstr"])["score"].transform(max)==final["score"] # grp by cluster keep clstr representative by highest score
        result_df = final[idx].drop_duplicates(subset=['clstr'], keep='first') # remove duplicate clstr representative if same score
        result_df.to_csv(clusters[i].split('_Genecallers')[0]+"_result_df.tsv", sep = "\t")

def GeneOverlap(prefix, score):
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
    print(f"{col.colours.BGreen}Creating a penality system for genes that overlap{col.colours.Colour_Off}")

    contig_list = glob.glob("*_result_df.tsv")
    df_sorted = pd.DataFrame()
    df_final = pd.DataFrame()

    for contig in range(len(contig_list)):
        print(f"{col.colours.BGreen}Analysing {contig_list[contig].split('_result')[0]}{col.colours.Colour_Off}")
    #    print(f'{col.colours.BBlue}Checking if genes with a cluster size of {len(list_of_datasets) - i} overlap with genes that have a cluster size of {dataset["clstr_size"].iloc[0]}{col.colours.Colour_Off}')
        result_df = pd.read_csv(contig_list[contig], sep = "\t")

        result_df["operon"] = 0
        result_df["overlap"] = result_df["length"]*-1

        list_of_datasets = [result_df[result_df["clstr_size"] ==  i ].reset_index(drop=True) for i in reversed([n+1 for n in range(result_df["clstr_size"].max())])]
        for i in range(len(list_of_datasets)):
            print(f"""{col.colours.BGreen}Analysing consenus {list_of_datasets[i]["clstr_size"][0]}{col.colours.Colour_Off}""")
            while list_of_datasets[i].empty:
                break
            else:
                for INDEX, ROW in list_of_datasets[i].iterrows():
                    for dataset in list_of_datasets[i:]:
                        for index, row in dataset.iterrows():
                            range1 = pd.Index(range(list_of_datasets[i][["Sstart", "Send"]].min(axis=1)[INDEX], list_of_datasets[i][["Sstart", "Send"]].max(axis=1)[INDEX]))
                            range2 = pd.Index(range(dataset[["Sstart", "Send"]].min(axis=1)[index], dataset[["Sstart", "Send"]].max(axis=1)[index]))
                            if len(range1.intersection(range2)) == [4,8]:
                                dataset["operon"][index] += 1
                            elif list_of_datasets[i][["Sstart", "Send"]].min(axis=1)[INDEX] == dataset[["Sstart", "Send"]].max(axis=1)[index]:
                                dataset["operon"][index] += 1
                            elif list_of_datasets[i][["Sstart", "Send"]].max(axis=1)[INDEX] == dataset[["Sstart", "Send"]].min(axis=1)[index]:
                                dataset["operon"][index] += 1
                            else:
                                dataset["overlap"][index] += len(range1.intersection(range2))
        results_df = reduce(lambda x,y: pd.merge(x,y, how='outer'), list_of_datasets)
        print(f"{col.colours.BGreen}Analysing {contig_list[contig].split('_result')[0]} done!{col.colours.Colour_Off}")

        # Create a an overlap scoring bin
        print(f"{col.colours.BGreen}Penalising genes that overlap{col.colours.Colour_Off}")
        score_bins = [-999999, 10, 40, 70, 100, 999999]
        penality = [0, -1, -2, -3, -9]
        results_df["overlap_pen"] = pd.cut(results_df["overlap"], score_bins, labels=penality)

        # Create a total scoring system based on the previous scoring system plus the operon and overlap penality
        print(f"{col.colours.BGreen}Creating the final scoring system for representative genes of each cluster{col.colours.Colour_Off}")
        results_df["total_score"] = results_df[["overlap_pen", "operon", "score"]].sum(axis = 1)

        # Sort the dataframe by gene start
        print(f"{col.colours.BGreen}Mark 1{col.colours.Colour_Off}")
        sorted_df = results_df.sort_values(by=["Sstart"]).reset_index(drop=True)

        # Only include the gene if the score is greater than 3
        print(f"{col.colours.BGreen}Mark 2{col.colours.Colour_Off}")
        final_df = sorted_df[sorted_df["total_score"] >= score]

        # Save the sorted_df to file
        final_df.to_csv(f"""{contig_list[contig].split('_result')[0]}_{score}_gene_stats.tsv""", sep = "\t", index = None, header = True)
        print(f"{col.colours.BGreen}Mark 3{col.colours.Colour_Off}")
        Genes = f"""rm {contig_list[contig].split('_result')[0]}_contig_seaphages.fna;
            for i in $(cut -f2 {contig_list[contig].split('_result')[0]}_{score}_gene_stats.tsv); do \
            echo $i | seqtk subseq {contig_list[contig].split('_result')[0]}*_Genecallers_combined.fa - >> {contig_list[contig].split('_result')[0]}_contig_seaphages.fna;
            done;
            rm {contig_list[contig].split('_result')[0]}_{score}_gene_stats.tsv \
            {contig_list[contig].split('_result')[0]}_Genecallers_combined.fa \
            {contig_list[contig].split('_result')[0]}_testing.fa \
            {contig_list[contig].split('_result')[0]}_result_df.tsv \
            {contig_list[contig].split('_result')[0]}*_genecalls.tsv \
            {contig_list[contig].split('_result')[0]}_Genecallers_combined_clstr.txt"""

        stdout, stderr = execute(Genes)
        print(f"{col.colours.BGreen}Mark 4{col.colours.Colour_Off}")

        df_sorted = sorted_df.append(sorted_df)
        df_final = final_df.append(final_df)
        print(f"{col.colours.BGreen}Mark 5{col.colours.Colour_Off}")

    print(f"{col.colours.BGreen}Only keeping genes with a score higher then {score}{col.colours.Colour_Off}")
    df_sorted.to_csv(f"""{prefix}_all_scores_gene_stats.tsv""", sep = "\t", index = None, header = True)
    df_final.to_csv(f"""{prefix}_{score}_gene_stats.tsv""", sep = "\t", index = None, header = True)

    print(f"{col.colours.BGreen}Gene statistics can be found here: {prefix}_gene_stats.tsv{col.colours.Colour_Off}")

    Combine = f"""cat *_contig_seaphages.fna > {prefix}_seaphages.fna; rm *_contig_seaphages.fna"""
    stdout, stderr = execute(Combine)
    print(f"{col.colours.BGreen}Gene calls can be found here: {prefix}_seaphages.fna{col.colours.Colour_Off}")


def execute(bash):
    process = subprocess.Popen(bash, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    (stdout, stderr) = process.communicate()

    print(stdout.decode("utf-8"), stderr.decode("utf-8")) # this should go to a log file instead
    return stdout, stderr


def main(args):
    runGeneCallers(args.infile)
    renameGeneCalls()
    getGeneStats(args.infile)
    runCD_hit(args.infile)
    combineGeneClusters()
    GeneOverlap(args.prefix, args.score)

if __name__=="__main__":
    args = my_args()
    main(args)
