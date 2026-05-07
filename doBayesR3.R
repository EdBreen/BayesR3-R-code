# BayesR3 by Edmond J. Breen 01/04/2022
# Department of Jobs Precincts and Regions Victoria
# Agriculture Victoria

#########################################################################################
# This code is for illustrative purpose and not really suitable for large data sets 
# Its main purpose is to demonstrate the block Gibbs sampling version of Bayes R known as bayesR3
# It to accompany the manuscript:

    # BayesR3: fast MCMC blocked processing for large scale 
    # multi-trait genomic prediction and QTN mapping analysis
    # 
    # E.J. Breen.1, I.M. MacLeod.1, P.N. Ho1, M. Haile-Mariam.1, 
    # J.E. Pryce.1,2, C.D. Thomas.1, H.D. Daetwyler.1,2, and M.E. Goddard.1,3
    # 
    # 1 Agriculture Victoria, AgriBio, Centre for AgriBioscience, Bundoora, Victoria 3083, Australia 
    # 2 School of Applied Systems Biology, La Trobe University, Bundoora, Victoria 3083, Australia 
    # 3 Faculty of Veterinary and Agricultural Sciences, University of Melbourne, 
    #  Parkville, VIC, 3052 Australia

#########################################################################################
# BayesR3 uses Gibbs sampling to solve for the following mixed effects genetic model
#
#	          y=Xu+Vg+Za+e	
#
# where y is column vector of phenotype values, of length n_R, which is the number of records. 
#
# X is a (n_R X n_F) fixed effects design matrix and u is a column vector of length n_F of fixed effects, 
# which are assumed to be normally distributed u_f~N(0,σ_f^2). 
#
# V is a coded genotype (n_R X n_M) matrix, with codes 0,1 and 2 and n_M is the number of SNPs. 
# V represents the observed genotypes for each individual/record across n_M markers.  
# Therefore, g is a column vector of length n_M containing the SNP effects.
#
# Z is an identity matrix (n_R X n_R) and a is a column vector of random genetic effects not explained by the SNPs 
# with polgenic variance  σ_a^2; such that a~N(0,Aσ_a^2), 
# and A is the additive relationship matrix.  Here A will be generated from a user provided pedigree.
# Also it is possible to provide a GRM file, such as the G-matrix presented by 
# VanRaden PM. Efficient methods to compute genomic predictions. J Dairy Sci. 2008;91(11):4414-23.
# Also, the sampling of the polygenic genetic effects is recommended when the SNPs do not explain all 
# the genetic variance or when it is desired to fit all the SNPs with small effects as well as a small 
# number of SNPs with larger effects. In the former case one can use the A matrix derived 
# from the pedigree and in the latter case a GRM based on all the SNPs 
#
# Note also, that e ~ N(0,W^(-1) σ_e^2), 
#
# W is a diagonal weight matrix and σ_e^2 is the residual variance.
#
# Note both the fixed effects component, Xu, and the polygenic component Za are optional.
#
# BayesR3 breaks up the Vg component into n_B non-overlapping blocks such that:
#    
#  Vg=V_1g_1 + V_2g_2 + ... V_n_Bg_n_B
# 
# The number of blocks, n_B, is determined from by the block size, bs, and therefore,  n_B, 
# is the least integer greater or equal to n_M/bs. Where n_M is the number of SNPs
# All blocks are the same size, except the last block which can be smaller.
#
# The variances Va, Ve and Vg are sampled from inverse-chisquare distributions.

#########################################################################################
# INPUT FILE FORMATS
# All input files are expected to be single-space delimited plain ASCII text files. 
# Missing values in all files types can be represented by either a comma , an asterisk * or NA

#-------------------
# The phenotype file
# A phenotype file is a single-space delimited text file (txt), with no header line. 
# All the phenotype data for a record should be presented on one line. 
# The file should contain a minimum of at least two columns.
# The first column (column 1) should contain the record IDs. 
# Subsequent columns should contain the phenotypes (traits) for each record, 
# and any corresponding weights. Often phenotype columns are interlaced with weight columns

#-------------------
# The reference file
# This is an optional text file, and when specified it is expected to contain only the ID’s 
# for the records to include in the analysis. The file is expected to contain a single column of IDs; 
# that is, one ID per line. However, if the file contains multiple columns of text, 
# then they are expected to be separated by spaces but only the first column of text
# is used as the identifying reference set. 
# The reference file is used to select a subset of animals to analyze rather than analyzing all
# the animals in the phenotype type. It is useful for doing training and validation.

#-------------------
# The genotype file
# This file is also a single-space delimited text file with no header line. 
# All the genotype data for a record should be presented on one line. 
# The line (row) order of the genotype file does not need to match the order given in the phenotype file. 
# The first column in the genotype file is expected to hold the record IDs. Subsequent columns hold the genotypes.
# Only records corresponding to records specified in the phenotype file or a reference file (see below) 
# will be extracted and used in the analysis. If the genotype file has more records than the specified in
# the phenotype file or reference file, only those records matching the phenotype/reference records are used.  

#-------------------
# Fixed effects file
# This file is optional, it is also a single-space delimited text file with no header. 
# All the fixed effects data for a record should be presented on one line.
# The first column contains the record IDs. Subsequent columns specify the fixed 
# effects associated with each ID. It is assumed that the fixed effects is actually
# a design matrix as made by R.

#-------------------
# The Pedigree file
# This file is optional, it is also a single-space delimited text file with no header line, 
# where the first three columns are expected to represent record ID, sire ID and dam ID and in that order. 
# It is also assumed that youngest individuals are towards the bottom of the file, while the older individuals 
# are towards the top. See http://rstudio-pubs-static.s3.amazonaws.com/382572_35603135de154a79a72b9f04064be2f2.html

#-------------------
# The GRM file
# A GRM file or an A-matrix file is expected to be a single-space delimited file without any header line. 
# It is also expected to contain a n×n matrix of record relationships as determined from some set of SNPs, 
# or from some pedigree.  This GRM/A-matrix file is also expected to contain n+1 column, where the first
# column contains the record IDs. 
# 
# The GRM matrix should be as specified by:
# VanRaden PM. Efficient methods to compute genomic predictions. J Dairy Sci. 2008;91(11):4414-23.

#########################################################################################
# EXAMPLE DATA SET:
# A small example data set is provided using simulated data that was supplied to 
# the 14th QTL-MAS workshop (http://jay.up.poznan.pl/qtlmas2010/index.html). 
# Szydlowski, M., Paczyńska, P. QTLMAS 2010: simulated dataset. BMC Proc 5, S3 (2011). 
# (https://doi.org/10.1186/1753-6561-5-S3-S3)
# 
# In brief: an out breed population was simulated using the LDSO software (Ytournel, 2008), 
# with 1000 generations of 1000 individuals, followed by 30 generations of 150 individuals.
# Nine thousand nine hundred ninety SNP markers were distributed on 5 chromosomes. Each chromosome
# had a size of 1 Morgan and carried 1998 SNP equally distributed (1 SNP every 0.05 cM).
# Data correspond to the last generations of this pedigree, i.e. 20 sire families 
# with 10 mates by sire and 15 progeny per dam.

#----------------------------------------
# The pedigree file [ pedigree.txt, 47K]
# It consists of 3226 individuals in 5 generations (F0-F4). There are 20 founders: 5 males and 15 females. 
# A female mates once and gives birth to about 30 progeny. A pedigree file is in the form 
# "ID MaleParent FemaleParent Sex". M=male, F=female, 0=missing parent.

#----------------------------------------
# The phenotype file: [simPhen.txt, 31K]
# Contains 2326 records, each with 2 traits: trait Q is a quantitative trait, and 
# trait B is a binary trait. Note BayesR3 doesn't inherently work with binary traits.
# Note. Young individuals (generation F4: individuals 2327 to 3226) have no phenotypic records. 

#----------------------------------------
# The genotype file [simGen.txt, 45M]
# It contains the genotypes associated with the phenotype records. 
# It contains 10032 columns, where the first column is record ID. 


#----------------------------------------
# The map file: [marker-info.txt, 165]
# The position of each SNP is included in the 'marker-info' file. 
# Column 1 is for SNP marker number, column 2 is chromosome number where the SNP is located, 
# column 3 is the chromosomal position in bp. 

##########################################################################################
# WORK FLOW:
#
#           getData --> imputeGT --> centreScale --> makeBlocks --> bayesR3
#
# The first thing is to load the data for analysis using getData
# The arguments that are passed to getData are as explained:
#
# getData: input parameters
#     genofile          genotype file name and path   
#     phenofile         phenotype file name and path 
#     tc                trait column in phenotype file 
#     wtc = NULL        optional weight column in phenotype file
#     reffile = NULL    optional reference file and path
#     fixedfile = NULL  optional fixed effect file name and path
#     pedfile = NULL    optional pedigree file name and path
#     grmfile = NULL    optional GRM or A matrix file  name and path
#     
# Example call: getData(genofile = "./simGen.txt", phenofile = "./simPhen.txt", tc = 2) 
# 
# getData exports the following variables to the hosting environment
#     geno:   matrix of genotypes used by makeBlocks and bayesR3
#     pheno:  vector of phenotypes used by bayesR3
#     wts: a  vector of weights used by makeBlocks and bayesR3
#     fixed:  fixed effects matrix, used by bayesR3  
#     iA22:   an optional inverse relationship matrix produced from a pedigree (ID, sire, dam),
#              or generated from a GRM or A-matrix 

# Once the data is loaded call imputeGT to impute any missing genotypes. 
# Missing genotypes are imputed as the mean of that SNP
#
# Example call: imputeGT() 

# Once the genotypes are complete then centre and scale using centreScale.
#
# centreScale
#     minFreq=0.002  # MAF minimal allele frequency. SNPs with a minor allele frequency 
#                      less than minFreq are set to 0
#  
# Example call: centreScale(minFreq=0.01)
#
# centreScale exports the following variables to the hosting environment
#   mono:   A vector of allele frequencies of length n_M
#   geno:   It modifies the geno matrix of the hosting environment
 
# Once the genotypes are centred and scaled call makeBlocks
# makeBlocks will create the blocks needed by bayesR3,  each of which will a bs X bs matrix.
# 
# makeBlocks input parameter
#         bs = NULL: optional value, as it allows users to define their own block sizes.
#              However, if left undefined makeBlocks will use the sqrt of the # of records
#  
# Example call: makeBlocks() 
# 
# makeBlocks exports the following variables to the hosting environment
#        bs: block size, 
#        blocks: the bayesR3

# Once the blocks have been made then bayesR3 is ready to be called
#
# Next call bayesR3, passing it the number of iterations, plus the heritability of the trait if known.
#  
# Example call  br3 = bayesR3(it = 10000) # this will take ~10 minutes
#
# timestamp(); a = bayesR3(10000); timestamp() minutes
##------ Wed Apr 06 05:28:51 2022 ------##
##------ Wed Apr 06 05:39:31 2022 ------##                

# bayesR3 
#         it=1000,              # Number of iterations
#         burnIn = 0.5 * it,    # Number of burn-in iterations
#         h2 = 0.5,             # heritability estimate
#         mix.pi = c(0.94,0.049,0.01,0.001),     # Proportion in each SNP effect distribution
#         mix.alloc = c(0.0,0.0001,0.001,0.01),  # Allocation of additive genetic variance
#         mix.alpha = c(1,1,1,1),                # Dirichlet prior
#         mix.nd = length(mix.pi)
#----------------------------------------
# The return value from bayes3 is an invisible list containing
# 
# lst = list(s = c(r2 = R-square from lm between the yhats and phenotypes (y-values)
#                  mu = mean 
#                  Vp = weighted phenotype variance
#                  Va = genetic variance explained by the SNPs 
#                  Ve = residual variance
#                  Vg = polygenic variance. Variance not explained by SNPs)
#            nMono = scalar representing the number of monomorphic SNPs
#            Nk = vector of SNPs count falling into mixture components 1,2,3 and 4,
#            Vk = vector of mixture component variances
#            yhat = model yhat
#            PiP =  posterior inclusion probabilities for the mixture components for each SNP
#            effects = model SNP effects
#                   
# If fixed effects were included for analysis
#
# lst = append(lst, 
#                 list(fixed = fixed effects, 
#                      Vf = variance of the fixed effects))
#
# If genomic relationships were included for analysis
#            
# lst = append(lst, list(u = polygenic effects))

#----------------------------------------
# The code

if(!require(DirichletReg)) {
    install.packages("DirichletReg")
    require(DirichletReg)
}
if(!exists("createA")) # from http://rstudio-pubs-static.s3.amazonaws.com/382572_35603135de154a79a72b9f04064be2f2.html
    source(file = "./createA.R")

exit <- function(msg) { print(msg); invokeRestart("abort") } 


getData = function(genofile ,    
                   phenofile ,
                   tc,
                   wtc = NULL,
                   reffile = NULL,
                   fixedfile = NULL,
                   pedfile = NULL,
                   grmfile = NULL)
    
{
    # This function reads both file.txt.gz and normal file.txt files
    
    # Automatically aligns genotypes, phenotypes, fiexed effects rows by ID
    # Assumes genotypes are 0,1,2 encoded
    # Missing genotypes are imputed using average SNP genotype
    # Missing phenotypes are used by setting their weight to 0
    # Missing values are one of: c(".", "*", "NA")
    
    doRef = !is.null(reffile)
    
    if(doRef) 
        ref = read.table(gzfile(reffile), sep=' ', header=F)[,1]
    
    pheno <<- read.table(gzfile(phenofile),sep = ' ',header = F,
                   na.strings = c(".", "*", "NA"))
    
    if(doRef) {
        id = match(y[,1], ref) # order according to phenotype file
        if(length(y[!is.na(id),1]) != length(ref))
            exit("Error: Reference record(s) not found in phenotype file")
        
        pheno <<- pheno[id,]
        
    } else {
        ref = pheno[,1]
        doRef = T
    }
    
    
    # read and create weight vector and adjust phenotypes appropriately
    if(!is.null(wtc)) 
        wts <<- pheno[,wtc]   
    else
        wts <<-  rep(1, length(pheno[,1]))
    
    id <<- is.na(pheno[,tc])
    pheno[id,tc] <<- 0
    wts[id] <<- 0
    
    pheno <<- as.numeric(pheno[, tc])
    
    
    geno <<-  read.table(gzfile(genofile), sep = ' ', header = F,
                    na.strings = c(".", "*", "NA"))
    
    #freq <<- colSums(m[,-1]) / (2 * nrow(m))
    
    id = match(ref,geno[,1])
    if(any(is.na(id))) 
        exit("Error: reference phenotype record(s) not found in  genotypefile")
    
    
    geno <<- geno[id,]
    
    colnames(geno) <<- NULL
    geno <<- data.matrix(geno[,-1]) * 1.0 # convert to float point

    
    
    
    if (!(is.null(fixedfile))) {
        fx = read.table(gzfile(fixedfile), sep = ' ', header = F) # no missing values allowed
        id = match(ref, fx[, 1])
        if (any(is.na(id)))
            exit("Error: reference phenotype record(s) not found in  fixed effects file")
        fx = fx[id,]
        fx = data.matrix(fx[,-1])
        fixed <<- fx
    } else
        fixed <<- NULL
    
    
    g = NULL
    
    if(!(is.null(pedfile))) {
        ped <<- read.table(gzfile(pedfile), sep = ' ', header = F)
        
        id = match(ref,ped[,1])
        if (any(is.na(id)))
            exit("Error: reference phenotype record(s) not found in  pedigree file")
        
        A = createA(ped[,1:3]) # only include first 3 cols
        iA22 <<- solve(data.matrix(A[id, id]))
        rm(A)
        
    } else if (!(is.null(grmfile))) {
        g = read.table(gzfile(grmfile), sep = ' ', header = F)
        id = match(ref, g[, 1])
        if (any(is.na(id)))
            exit("Error: reference phenotype record(s) not found in  grm file")
        
        iA22 <<- solve(data.matrix(g[id, id + 1])) #skip ID column
        rm(g)
    } else 
        iA22 <<- NULL
}
imputeGT = function()
{
    # impute missing genotypes using average SNP genotypes
    
    id = apply(geno,2, function(x) any(is.na(x)))
    if(length(which(id)) > 0) {
        #impute missing genotypes
        for(i in which(id)) {
            j = is.na(geno[,i])
            geno[j,i] <<- mean(geno[!j,i])
        }
    }   
}
centreScale = function(minFreq=0.002)
{
    nc = ncol(geno)
    nr = nrow(geno)
    maxFreq = 1 - minFreq
    
    p = colSums(geno) / (2 * nr)
    mono <<-  p <= minFreq | p >= maxFreq
    geno[, mono] <<-  0
    
    p2 = 2 * p
    psqrt = vector('numeric', length = nc)
    psqrt = sqrt(p2[!mono] * (1.0 - p[!mono]))
    p2 = p2[!mono]
    
    #geno[,!mono] <<-  sweep(geno[,!mono],1,p2,FUN='-')
    
    id = !mono
    for(i in 1:nrow(geno))
        geno[i,id] <<- (geno[i,id] - p2) / psqrt
    
    
}
makeBlocks=function(bs = NULL)
{
    # geno is a matrix and bs is block size
    # wts is a vector of weights
    
    if(is.null(bs))
        bs = as.integer(min(sqrt(length(pheno)), sqrt(ncol(geno))))
    
    bs <<- bs
    nc = ncol(geno)
    nr = nrow(geno)
    
    s1 = seq(1,nc, by = bs)
    s2 = s1 + bs -1
    ln = length(s1)
    s2[ln] = nc
    
    blocks <<- list()
    for(i in 1:ln) {
        m = as.matrix(geno[,s1[i]:s2[i]])
        blocks[[i]] <<- t(m) %*% sweep(m, 1, wts, FUN="*") 
    }
}

bayesR3 = function(it=100,               # Number of iterations
                   burnIn = 0.5 * it,    # Number of burnin iterations
                   h2 = 0.5,             # heritability
                   mix.pi = c(0.94,0.049,0.01,0.001),     # Proportion in each SNP effect distribution
                   mix.alloc = c(0.0,0.0001,0.001,0.01),  # Allocation of additive genetic variance
                   mix.alpha = c(1,1,1,1),                # Dirichlet prior
                   mix.nd = length(mix.pi))
{
    pheno = as.vector(pheno)
    snpDist = matrix(0, ncol=mix.nd, nrow=nrow(geno))
    
    set.seed(1234567)
    
    nxc = ncol(geno)
    nxr = nrow(geno)
    
    
    bs = ncol(blocks[[1]]) # block size
    ln = length(blocks) # number of blocks
    
    s1 = seq(1,nxc, by = bs)
    s2 = s1 + bs -1
    s2[ln] = nxc
    
    
    id = wts != 0
    wsum = sum(wts)
    nw = length(which(id)) # number of none zero weights
    
    mu = sum(pheno * wts)/wsum
    vara = sum( (pheno[id] - mu)^2 / (h2 + (1-h2)/wts[id]) ) 
    phenoVar = vara /(nw-1)
    vare = phenoVar * (1 - h2)
    vara = phenoVar * h2
    
    
    varDistStoreAcc = vector(mode='numeric',length = mix.nd)
    
    #fixed effects
    nfc = 0
    dgf = NULL
    fx = NULL
    fxAcc = NULL
    doFixed = F
    doMean = T
    
    if(!is.null(fixed)) {
        
        nfc = ncol(fixed)
        dgf = vector(mode="double", nfc)
        for(i in 1:nfc)
            dgf[i] = t(fixed[,i]) %*% (fixed[,i] * wts)
        
        doMean = F
        doFixed = T
        
        fx = rep(0, nfc)
        fxAcc =  rep(0, nfc)
        Vf = rep(0, nfc)
    }
    
    vgAcc  =0
    if(!is.null(iA22)) {
        
        doiA22 = T
        
        t = vara
        vara = t * 0.9
        varg = t * 0.1
        
        nu = nrow(iA22)
        u = vector("numeric", nu)
        ucollect = u
    } else {
        varg = 0
        doiA22 = F
        ucollect = F
    }
    
    
    a = vector('numeric', nxc) # SNP effects 
    acollect = a
    
    
    e = pheno # initialize error vector
    mu = 0    # reset to zero
    
    # create some accumulators
    muAcc =  0
    vaAcc =  0 
    veAcc =  0 #rep(0, nyc) 
    
    
    mix.distFreq = rep(0, mix.nd)
    mix.vk = mix.alloc * vara     # initialize mixture variance vk
    
    
    snpDist = matrix(0, nrow = nxc, ncol=mix.nd) # counters
    varDistStore = rep(0, mix.nd)
    
    outer = as.integer(it/bs)
    inner = bs
    bjburn = burnIn / inner;
    
    inc1 = max(1.0,ln/inner) # increment for updating mean,fixed effects etc
    update = 0
    bj = 0
    cnt = 0
    cnt2 = 0
    
    for(l in 1:outer) { # the main outer cycle iterations
        snpVarSum = 0
        beta = rep(0, mix.nd)
        if(l >= bjburn) varDistStore = beta
        for (j in 1:ln) {  # iterate over the blocks
            if(bj >= update) { # time to update fixed effects?
                #sample mean
                if (doMean) {
                    mold = mu
                    bl = (sum(e*wts) + wsum * mold)/ wsum
                    mu = rnorm(1, bl, sqrt(vare / wsum))
                    e = e - (mu - mold)
                }
                #sample fixed effects
                if (doFixed) {
                    for (k in 1:nfc) {
                        fold = fx[k]
                        bl = ((fixed[, k] %*% (e*wts)) + dgf[k] * fold) / dgf[k]
                        fx[k] = rnorm(1, bl, sqrt(vare / dgf[k]))
                        adiff = fx[k] - fold
                        if(abs(adiff) > 1e-8)
                            e = e -  fixed[, k] * adiff
                    }
                    mu = fx[1]
                }
            }
            
            dg = diag(blocks[[j]])  # get diagonal
            inc = length(dg)        # number of SNPs in block
            jns = s1[j] - 1         # starting SNP offset of block
            aseg = a[s1[j]:s2[j]]   # make copy of block effects
            xe =  t(geno[,s1[j]:s2[j]]) %*% (e * wts) # create block residuals
            
            for(i1 in 1:inner) {    # inner cycle
                for (ij in 1:inc) { # sample block SNP effects
                    ai = ij + jns   # get current SNP position
                    if(mono[ai]) next # skip monomorphic SNPs
                    
                    aold = a[ai]
                    zz = dg[ij]
                    ve = vare / zz
                    rhs = xe[ij] + zz * aold
                    gls = rhs / zz
                    gls2 = gls ^ 2
                    v = mix.vk + ve
                    
                    #pjk = exp(log(mix.pi) - (0.5 * gls2/v + log (sqrt(v))))
                    pjk = mix.pi/ (exp(0.5 * gls2/v) * sqrt(v))
                    pjk = pjk/sum(pjk)
                    
                    tk = mix.nd
                    accum = 0
                    uRand = runif(1, min = 0, max = 1)
                    
                    for(k in 1:mix.nd) { # categorical type function
                        accum = accum + pjk[k]
                        if(accum>=uRand) {
                            tk = k
                            break # break out of for loop
                        }
                    }
                    
                    beta[tk] = beta[tk] + 1
                    if(l >= bjburn)
                        snpDist[ai, tk] = snpDist[ai, tk] + 1;
                    
                    if(tk == 1) {
                        a[ai] = 0
                    } else {
                        v1 =  zz + vare/mix.vk[tk]
                        a[ai] = rnorm(1, rhs/v1, sqrt(vare/v1))
                        if(l >= bjburn) {
                            acollect[ai] = acollect[ai] + a[ai]
                            varDistStore[tk] = varDistStore[tk] + a[ai] * a[ai];
                        }
                    }
                    # update residuals
                    adiff = aold - a[ai]
                    if(abs(adiff) > 10^-8) # test for almost zero
                        xe = xe + blocks[[j]][,ij] * adiff # update block residuals
                } #end ij
                
            }# end inner loop
            
            e = e - geno[,s1[j]:s2[j]] %*% (a[s1[j]:s2[j]] - aseg)
            if(bj >= update) {
                
                update = update  + inc1
                if(doiA22) {
                    lambda = vare/varg
                    k = 1
                    for(j in 1:nu) {
                        uold = u[j]
                        lhs = wts[k] + lambda * iA22[j,j]
                        m = (e[k] * wts[k] - lambda * sum(iA22[,j]*u))/lhs + uold
                        u[j] = rnorm(1, m, sqrt(vare/lhs))
                        e[k] = e[k] - (u[j] - uold)
                        k = k + 1
                    }
                    varg =   sum(u*(iA22 %*% u))/rchisq(1, nu-2)
                }
                vare = (sum(e * e * wts)) / rchisq(1, wsum - 2)
                cnt2 = cnt2+1
            }
            
            bj = bj+1 # processed block counter
        } # end of blocks (ln)
        
        beta = beta/inner
        
        mix.pi = rdirichlet(1, beta + mix.alpha)
        
        { # sample vara
            s = sum(beta[-1])  # number of SNPs in the model
            if (s > 2)
                vara = (s * sum(a ^ 2)) / rchisq(1, s - 2)
            mix.vk = vara * mix.alloc
        }
        
        if ((l >= bjburn) | (l == outer)) {
            vaAcc = vaAcc + vara
            veAcc = veAcc + vare
            vgAcc = vgAcc + varg
            
            muAcc = muAcc + mu
            varDistStoreAcc = varDistStoreAcc + varDistStore/inner
            cnt = cnt + 1
            
            if(doiA22)
                ucollect = ucollect + u
            
            if (doFixed) {
                diff = fx - fxAcc;
                fxAcc = fxAcc + diff/cnt
                Vf = Vf + (fx - fxAcc) * diff
            } 
            
        }
        
    } # end of outer loop
    
    # collect results
    
    mix.distFreq = round(colSums(snpDist)/(inner *cnt),2)
    snpDist = snpDist / rowSums(snpDist)
    snpDist[is.nan(snpDist[,1]),] = 0
    
    # varDist <<- varDistStoreAcc/cnt
    # vare <<- veAcc/cnt
    # vara <<- vaAcc/cnt
    # mu <<- muAcc/cnt
    # 
    
    u = ucollect/cnt
    
    a = acollect/(inner * cnt)
    
    yhat = (geno %*% a)[,1]
    
    if(doiA22)
        yhat  = yhat + u
    
    if(doFixed)
        yhat = yhat +  fixed %*% fxAcc 
    
    r2 = summary(lm(pheno[id]~yhat[id]))$adj.r.squared
    lst = list(s = c(r2 = r2, mu = muAcc/cnt, Vp = phenoVar, Va = vaAcc/cnt, Ve = veAcc/cnt, Vg=vgAcc/cnt),
               nMono = length(which(mono)),
               Nk = mix.distFreq, 
               Vk = varDistStoreAcc/cnt,
               
               yhat = yhat,
               PiP = snpDist,
               effects = a)
    
    if(doFixed)
        lst = append(lst, list(fixed = fxAcc, Vf))
    if(doiA22)
        lst = append(lst, list(u = u))
    
    invisible(lst)
    
}













