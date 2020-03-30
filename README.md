# SEAPHAGES

## Quickstart
SEAPHAGES.py -i <input_FASTA> -p <output_prefix>

Returns 2 files:
  <prefix>\_seaphages.fna which is the consensus GeneCallers
  <prefix>\_gene_stats.tsv which an a tab separated file containing the scoring of each gene from each cluster for manual curation

## Description
SEAPHAGES is a consenus genecaller, using multiple genecalling softwares to determine gene calls. It uses Prodigal, Glimmer, MGA, GeneMark and PHANOTATE to call genes and cd-hit to determine genecalls that are the same. Then using multiple criteria derived from <URL_TO_SEAPHAGES> it scores genes and selects resulting genes that are highly probably to be genes.

## Dependancies
Prodigal v2.6.3  
https://github.com/hyattpd/Prodigal  
`conda install -c bioconda prodigal`  

Glimmer v3.02b  
http://ccb.jhu.edu/software/glimmer/index.shtml  
`wget http://ccb.jhu.edu/software/glimmer/glimmer302b.tar.gz  
cd glimmer3.02/src
make`
Include /glimmer3.02/scripts and /glimmer3.02/bin in PATH
change awk shebang in /glimmer3.02/scripts/upstream-coords.awk and /glimmer3.02/scripts/get-motif-counts.awk to #!/usr/bin/awk

MGA (MetaGeneAnnotator) 2008/8/19  
http://metagene.nig.ac.jp/  
`conda install -c bioconda metagene_annotator`  

GeneMarkS v4.6b  
GeneMarkS2 v1.10  
http://exon.gatech.edu/GeneMark/license_download.cgi  
`cp gm_key_64 ~/.gm_key `

cd-hit v4.8.1  
http://weizhongli-lab.org/cd-hit/  
`conda install -c bioconda cd-hit`  

ELPH (Estimated Locations of Pattern Hits) v1.0.1  
https://cbcb.umd.edu/software/ELPH/  
`wget ftp://ftp.cbcb.umd.edu/pub/software/elph/ELPH-1.0.1.tar.gz
 tar -xzf ELPH.tar.gz
 cd ELPH/sources
 make`  

Seqtk v1.3-r106  
https://github.com/lh3/seqtk  
`git clone https://github.com/lh3/seqtk.git
cd seqtk
make`  

PHANOTATE v1.2.2  
https://github.com/deprekate/PHANOTATE  
`git clone --recursive https://github.com/deprekate/PHANOTATE.git; cd PHANOTATE; make`
conda install -c bioconda trnascan-se


### Python Modules  
argparse v1.4.0  
BioPython v1.72  
glob  
os  
pandas v0.24.1  
subprocess  
functools  
itertools

## Work to be done  
Include ARAGORN/ teRNAse / a tRNA caller  
Include gggenes for gene annotation  
Loop break when comparing different contigs in gene overlap  
