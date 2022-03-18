# Simulated data set
The data and information in this directory is from the [14th QTL-MAS](http://jay.up.poznan.pl/qtlmas2010/).
See also [https://bmcproc.biomedcentral.com/articles/10.1186/1753-6561-5-S3-S3] (https://bmcproc.biomedcentral.com/articles/10.1186/1753-6561-5-S3-S3) 
It represents simulated geneomic mouse data for public use and which has been edited and
formatted for bayesR3 usage:


## Pedigree [49.9 KB] (pedigree.txt)
The pedigree consists of 3226 individuals in 5 generations (F0-F4). There are 20 founders: 5 males and 15 females. 
A female mates once and gives birth to about 30 progeny. A pedigree file is in the form "ID MaleParent FemaleParent Sex". M=male, F=female, 0=missing parent. Note, the pedigree file was corrected March 18, 2010. In the previous version of the file, two individuals (89 and 374) appeared with misspecified gender. Correct sex for the individual 89 is F (female) and correct sex for the individual  374 is M (male). We thank Javad Nadaf for finding the errors.
## Phenotypes [30.7 KB] (simPhen.txt)
Two traits are observed. Trait Q is a quantitative trait, whereas trait B is a binary trait. Young individuals  (generation F4: individuals 2327 to 3226) have no phenotypic records.
##   Genome [5 MB] (simGen.txt.gz)
It is assumed that genome is about 500 mln bp long. It consists of 5 chromosomes, each of about 100 mln bp.  The genome was sequenced, hence a position of each observed SNP is known. The position of each SNP is included in the 'marker-info' file. Column 1 is for SNP marker number, column 2 is chromosome number where the SNP is located, column 3 is the chromosomal position in bp.
## Marker data [11 MB] (genotypes.txt.gz)
Each individual is genotyped for 10031 biallelic SNPs. There are 3226 rows in the genotype file, each row for a single individual. First column is ID for an individual. Next columns contain alleles.

### Problems
Geneticists may want to recover the genetic architecture of traits Q and B, and genetic links between the traits. 
Breeders may be interested in breeding value of each individual.

True simulated breeding values, phased SNP data and QTL positions can be found in Download.
If you have any questions concerning the dataset please contact: qtlmas@jay.up.poznan.pl 
