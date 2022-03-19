# This code was written by EdBreen
# Department of Jobs Precincts and Regions Victoria
# Agriculture Victoria
# This code is under development 

if(!require(DirichletReg))
    install.packages("DirichletReg")
if(!exists("createA"))
  source(file = "./createA.R")

# mixture model
#results
#a = 0      # vector of SNP effects
#fx =0      # vector of fixed effects
#vara =0    # accumulative SNP effects variances
#vare =0    # error variance
#mix.vk  =0 # distribution variances
#mono = 0   # vector of monomorphic SNPs mono[x] == T if monomorphic
#bs = 0     # block size
# sequence order for functions calls
# getData(); imputeGT(); centreScale() ;makeBlocks(); bayesR3()


exit <- function(msg) { print(msg); invokeRestart("abort") } 


exampleDriver = function(it = 500)
{
  
  rm(list=ls())
  
  print("Getting Data")
  getData(genofile = "./simData/simGen.txt.gz", 
          phenofile = './simData/simPhen.txt', 
          tc = 2,
          pedfile = NULL) #"./simData/pedigree.txt")
  
  
  print("Doing  imputeGT")
  imputeGT(); # impute any missing genotypes
  print("Doing centreScale")
  centreScale(); # centre and Scale genotypes
  
  print(paste("Number of monomorphics:", length(which(mono))))
  
  print("Making blocks")
  
  makeBlocks()
  
  print(paste("Running BayesR3 using", it , "iterations"))
  bayesR3(it=it)  #returns an invisible list
}


createMap = function()
{
  map = read.table("./simData/marker-info.txt", header=F, sep=' ')
  map$label = paste(map$V3, map$V4, sep=':')
  colnames(map)[c(3,4)] = c("chr", "pos")
  map[,c(5,3,4)]
}

crtGRM = function(genofile)
{
    # create a GRM matrix file from the genofile
    # see: VanRaden PM. Efficient methods to compute genomic predictions. 
    #                   J Dairy Sci. 2008 Nov;91(11):4414-23. 
    m  =  read.table(genofile, sep = ' ', header = F)
    ids = m[,1]
    colnames(m) = NULL
    m  = data.matrix(m[,-1]) # convert to float point
    nr = nrow(m)
    p = colSums(m) / (2 * nr)
    denom = 2*sum(p*(1-p))
    a = round((m %*% t(m))/ denom, 3)
    df = data.frame(cbind(ids, a))
    #write output GRM file
    require(tools)
    grmfile = paste(file_path_sans_ext(genofile), "GRM.",  file_ext(genofile), sep = '')
    # write output GRM file
    write.table(df, file = grmfile,col.names = F, row.names = F,sep=' ')
    df
}


# minFreq = 0.002
# h2 = 0.5  # heritability
# DAV = T 


getData = function(genofile ,
                   phenofile ,
                   tc,
                   wtc = NULL,
                   reffile = NULL,
                   fixedfile = NULL,
                   pedfile = NULL,
                   grmfile = NULL,
                   bs = NULL)
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
    
    y = read.table(gzfile(phenofile),sep = ' ',header = F,
                   na.strings = c(".", "*", "NA"))
    
    
    if(doRef) {
      id = match(y[,1], ref) # order according to phenotype file
      if(length(y[!is.na(id),1]) != length(ref))
        exit("Error: Reference record(s) not found in phenofile")
      
      y = y[!is.na(id),]
      ref = y[,1]
      
    } else {
      ref = y[,1]
      doRef = T
    }
    

    # read and create weight vector and adjust phenotypes appropriately
    if(!is.null(wtc)) 
      wts <<- y[,wtc]   
    else
      wts <<-  rep(1, length(y[,1]))
    
    id = is.na(y[,tc])
    y[id,tc] = 0
    wts[id] = 0
    
    pheno <<- as.numeric(y[, tc]); 
    rm(y)
    
    
    m =  read.table(gzfile(genofile), sep = ' ', header = F,
                         na.strings = c(".", "*", "NA"))
    
    id = match(ref,m[,1])
    if(any(is.na(id))) 
      exit("Error: reference phenotype record(s) not found in  genofile")
    m = m[id, ]
    colnames(m) = NULL
    geno <<- data.matrix(m[,-1]) * 1.0 # convert to float point
    rm(m)
   
    if(is.null(bs))
      bs <<- as.integer(min(sqrt(length(pheno)), sqrt(ncol(geno))))
    
    
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
    psqrt[!mono] = sqrt(p2[!mono] * (1.0 - p[!mono]))
    
    
    #geno[,!mono] <<-  sweep(geno[,!mono],1,p2,FUN='-')
    
    geno[, !mono] <<-  mapply('-', data.frame(geno[,!mono]), p2[!mono])
    geno[, !mono] <<-  mapply('/', data.frame(geno[,!mono]), psqrt[!mono])
    
   

}
makeBlocks=function()
{
    # geno is a matrix and bs is block size
    # wts is a vector of weights
    
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

bayesR3 = function(it=100, burnIn = 0.5 * it,
                   h2 = 0.5,
                   DAV = T,
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
                    
                    pjk = mix.pi/ (exp(0.5 * gls2/v) * sqrt(v))
                    pjk = pjk/sum(pjk)
                    
                    tk = mix.nd
                    accum = 0
                    uRand = runif(1, min = 0, max = 1)
                    
                    for(k in 1:mix.nd) { # multinomial distribution
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
        
        if (DAV) { # sample vara
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





