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


Glimmer v3.02b
http://ccb.jhu.edu/software/glimmer/index.shtml


MGA (MetaGeneAnnotator) 2008/8/19
http://metagene.nig.ac.jp/


GeneMarkS v4.6b
GeneMarkS2 v1.10
http://exon.gatech.edu/GeneMark/


cd-hit v4.8.1
http://weizhongli-lab.org/cd-hit/


ELPH (Estimated Locations of Pattern Hits) v1.0.1
https://cbcb.umd.edu/software/ELPH/


Seqtk v1.3-r106
https://github.com/lh3/seqtk

PHANOTATE v1.2.2
https://github.com/deprekate/PHANOTATE

### Python Modules
argparse v1.4.0

BioPython v1.72

glob

os

pandas v0.24.1

subprocess

functools

## Work to be done
Include ARAGORN/ teRNAse / a tRNA caller

Include gggenes for gene annotation

Loop break when comparing different contigs in gene overlap
