# The QTL-MAS 2010 simulated livestock-style pedigree dataset
The data and information in this directory is from the [14th QTL-MAS](http://jay.up.poznan.pl/qtlmas2010/).
See also (https://bmcproc.biomedcentral.com/articles/10.1186/1753-6561-5-S3-S3) 
It represents simulated geneomic data for public use and which has been edited and
formatted for bayesR3 usage and extended to include a VCF file representation and an associated map file.


## Pedigree [49.9 KB] (pedigree.txt)
The pedigree consists of 3226 individuals in 5 generations (F0-F4). There are 20 founders: 5 males and 15 females. 
A female mates once and gives birth to about 30 progeny. A pedigree file is in the form "ID MaleParent FemaleParent Sex". M=male, F=female, 0=missing parent. Note, the pedigree file was corrected March 18, 2010. In the previous version of the file, two individuals (89 and 374) appeared with misspecified gender. Correct sex for the individual 89 is F (female) and correct sex for the individual  374 is M (male). We thank Javad Nadaf for finding the errors.
## Phenotypes [30.7 KB] (simPhen.txt)
Two traits are observed. Trait Q is a quantitative trait, whereas trait B is a binary trait. Young individuals  (generation F4: individuals 2327 to 3226) have no phenotypic records.
##   Genome [5 MB] (simGen.txt.gz)
It is assumed that genome is about 500 mln bp long. It consists of 5 chromosomes, each of about 100 mln bp.  The genome was sequenced, hence a position of each observed SNP is known. The position of each SNP is included in the 'marker-info' file. Column 1 is for SNP marker number, column 2 is chromosome number where the SNP is located, column 3 is the chromosomal position in bp.
## Marker data [11 MB] (genotypes.txt.gz)
Each individual is genotyped for 10031 biallelic SNPs. There are 3226 rows in the genotype file, each row for a single individual. First column is ID for an individual. Next columns contain alleles.

## simVCF.vcf.gz
`simVCF.vcf` 
is a phased VCF representation of the full QTL-MAS 2010 `simData` SNP genotype panel. 
It is generated from `genotypes-true-phase.txt.gz` and `marker-info.txt`, preserving all individuals and all biallelic SNP markers in the phased form `paternal|maternal`. 
The genotype file provides the paternal and maternal allele for each marker for each individual, while `marker-info.txt` provides marker identifier, chromosome, and physical position. 
A companion `simVCF_cm.map` file can be generated using the simulation assumption `1 Mb = 1 cM`, with cM positions calculated within chromosome. 
This VCF is intended as a whole-dataset phased input for PBWT/TPBWT testing.
## simVCF_cm.map.gz

`simVCF_cm.map` is a companion genetic-map file for `simVCF.vcf`. It is generated from `marker-info.txt` using the QTL-MAS simulation assumption `1 Mb = 1 cM` between adjacent markers. The file stores chromosome-local cumulative cM positions, with the first marker on each chromosome assigned `0 cM`.

This map is intended for use with the whole-dataset phased `simVCF.vcf` in PBWT/TPBWT testing and related tools that expect marker-level cM coordinates.

The file has four tab-separated columns and no header:

- column 1: chromosome name, with `chr` prefix added if needed
- column 2: marker ID from `marker-info.txt`
- column 3: chromosome-local cumulative cM position
- column 4: physical base-pair position from `marker-info.txt`

If recombination fractions are required, they should be calculated from adjacent interval distances, not from the cumulative cM value itself. For adjacent markers `i−1` and `i`:

	d_cM = cM_i - cM_(i-1)

Convert to Morgans:

	d_M = d_cM / 100

Then, for example, Haldane’s mapping function gives:

	r = 0.5 * (1 - exp(-2d_M))

and Kosambi’s mapping function gives:

	r = 0.5 * tanh(2d_M)

### Problems
Geneticists may want to recover the genetic architecture of traits Q and B, and genetic links between the traits. 
Breeders may be interested in breeding value of each individual.

True simulated breeding values, phased SNP data and QTL positions can be found in Download.
If you have any questions concerning the dataset please contact: qtlmas@jay.up.poznan.pl 
