## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  autodep = TRUE
)

## ----setup--------------------------------------------------------------------
library(gJLS2)

## ----gJLS2data----------------------------------------------------------------
data("chrXdat")
head(chrXdat)

## ----plot, fig.width=8, fig.height=6------------------------------------------
library(ggplot2)
ggplot(data = chrXdat, aes(x=PHENOTYPE, fill=as.factor(SEX))) + 
  geom_histogram(aes(y = ..density..),alpha=0.2, position="identity", 
                 binwidth = 0.2) + geom_density(alpha=0.2, size=0.2) +
  ggtitle("Historagm of Quantitative Trait") + xlab("Phenotype")

## ----summary------------------------------------------------------------------
print(summary(chrXdat$PHENOTYPE))
mean(chrXdat$PHENOTYPE); sd(chrXdat$PHENOTYPE)
library(moments)
skewness(chrXdat$PHENOTYPE)
kurtosis(chrXdat$PHENOTYPE)

## ---- echo=FALSE--------------------------------------------------------------
dt <- data.frame("CHR" = 23, "SNP" = c("rs5983012", 
                                      "rs986810",
                                      "rs180495",
                                      "rs5911042",
                                      "rs4119090"),
                 "A1" = c("A", "C", "G", "T", "G"),
                 "MAF" = c(0.1016, 0.2033, 0.3008, 0.4024, 0.4451))
print(dt)

## ----gJLS2_example_loc--------------------------------------------------------
locReg(GENO=chrXdat[,7:11], SEX=chrXdat$SEX, Y=chrXdat$PHENOTYPE, Xchr=TRUE);

## ----gJLS2_example_scale------------------------------------------------------
scaleReg(GENO=chrXdat[,7:11], SEX=chrXdat$SEX, Y=chrXdat$PHENOTYPE, Xchr=TRUE)
scaleReg(GENO=chrXdat[,7:11], SEX=chrXdat$SEX, Y=chrXdat$PHENOTYPE, Xchr=TRUE, loc_alg="OLS")

## ----gJLS2_example_gJLS-------------------------------------------------------
gJLS2(GENO=chrXdat[,7:11], Y=chrXdat$PHENOTYPE, SEX=chrXdat$SEX, Xchr=TRUE)

## ----dosage-------------------------------------------------------------------
library("MCMCpack")
N <- 300 
geno <- rbinom(N, 2, 0.3)
a <- 0.3 ## uncertainty
genPP <- rbind(rdirichlet(sum(geno==0),c(a,(1-a)/2,(1-a)/2)), 
        rdirichlet(sum(geno==1),c((1-a)/2,a,(1-a)/2)),
        rdirichlet(sum(geno==2),c((1-a)/2,(1-a)/2,a)))
head(genPP);
summary(rowSums(genPP))

## ----dosage_exm1--------------------------------------------------------------
sex <- rbinom(N, 1, 0.5)+1 ## using PLINK coding
y <- rnorm(N)
covar <- matrix(rnorm(N*10), ncol=10)
gJLS2(GENO=list(genPP), SEX=sex, Y=y, COVAR=covar, genotypic = TRUE) ## geno probabilities
gJLS2(GENO=list(genPP), SEX=sex, Y=y, COVAR=covar) ## geno dosage
try(gJLS2(GENO=list(genPP), SEX=sex, Y=y, COVAR=covar, origLev = TRUE)) ## cannot perform Levene's test

## ----related_samples----------------------------------------------------------
gJLS2(GENO=geno, SEX=sex, Y=y, COVAR=covar, related=TRUE, clust = rep(1:3, c(N/2, N/4, N/4)))

## ----xchr---------------------------------------------------------------------
genoX <- NA
genoX[sex==2] <- rbinom(sum(sex==2), 2, 0.3)
genoX[sex==1] <- rbinom(sum(sex==1), 1, 0.3)
table(genoX, sex)

## ----xchr_loc-----------------------------------------------------------------
locReg(GENO=genoX, SEX=sex, Y=y, COVAR=covar, Xchr=TRUE)

## ----xchr_scale---------------------------------------------------------------
gJLS2(GENO=genoX, SEX=sex, Y=y, COVAR=covar, Xchr=TRUE)

## ----chr_scale----------------------------------------------------------------
gJLS2(GENO=geno, SEX=sex, Y=y, COVAR=covar, origLev=TRUE)

## ---- eval=FALSE--------------------------------------------------------------
#  #!/usr/bin/env Rscript
#  
#  if("optparse" %in% rownames(installed.packages()) == FALSE) {
#  print("optparse not installed, trying to intall now ...")
#  install.packages("optparse", repos='http://cran.us.r-project.org')
#  }
#  
#  require("optparse")
#  
#  option_list = list(
#    make_option(c("-b", "--bfile"), type="character", default=NULL,
#                help="genotype dataset file name", metavar="character"),
#    make_option(c("-p", "--pfile"), type="character", default=NULL,
#                help="pheno and covariate dataset file name", metavar="character"),
#    make_option(c("-m", "--pheno"), type="character", default=NULL,
#                help="phenotype name", metavar="vector"),
#    make_option(c("-c", "--covar"), type="character", default=NULL,
#                help="covariate name", metavar="vector"),
#    make_option(c("-q", "--center"), type="character", default="median",
#                help="center option in gJLS2", metavar="character"),
#    make_option(c("-g", "--genotypic"), type="logical", default="TRUE",
#                help="genotypic option in gJLS2", metavar="character"),
#    make_option(c("-t", "--transform"), type="logical", default="FALSE",
#                help="transform option in gJLS2", metavar="character"),
#    make_option(c("-x", "--Xchr"), type="logical", default="FALSE",
#                help="Xchr option in gJLS2", metavar="character"),
#    make_option(c("-w", "--write"), type="integer", default= 50,
#                help="writing in chunk size", metavar="integer"),
#    make_option(c("-o", "--out"), type="character", default="out.txt",
#                help="output file name [default= %default]", metavar="character")
#  );
#  
#  
#  opt_parser <-OptionParser(option_list=option_list)
#  arguments <- parse_args (opt_parser, positional_arguments=TRUE)
#  opt <- arguments$options
#  args <- arguments$args
#  
#  bimf <- opt$bfile
#  phenof <-opt$pfile
#  
#  phenoNames <- opt$pheno
#  covarNames <- strsplit(opt$covar, ",")[[1]]
#  chunk_size <- opt$write
#  cat(paste("Writing in chunk size of", chunk_size, "\n"))
#  
#  ## additional options:
#  
#  centre <- opt$center;
#  cat(paste("Using center option", centre, "\n"))
#  genotypic <- opt$genotypic
#  cat(paste("Using genotypic option", genotypic, "\n"))
#  transform <- opt$transform
#  cat(paste("Using transform option", transform, "\n"))
#  xchr <- opt$Xchr
#  cat(paste("Using Xchr option", xchr, "\n"))
#  out <- opt$out
#  
#  
#  if("gJLS2" %in% rownames(installed.packages()) == FALSE) {
#  cat("gJLS2 not installed, trying to intall now ...")
#  install.packages("gJLS2", repos='http://cran.us.r-project.org')
#  }
#  
#  require("gJLS2")
#  
#  
#  ## checking pheno file
#  
#  
#  if("BEDMatrix" %in% rownames(installed.packages()) == FALSE) {
#  print("BEDMatrix not installed, trying to intall now ...")
#  install.packages("BEDMatrix", repos='http://cran.us.r-project.org')
#  }
#  
#  if("BGData" %in% rownames(installed.packages()) == FALSE) {
#  print("BGData not installed, trying to intall now ...")
#  install.packages("BGData", repos='http://cran.us.r-project.org', dependencies=T)
#  }
#  
#  ## checking inputs to be bed, fam, bim files
#  
#  require("BGData")
#  require("BEDMatrix")
#  bedFiles <- BEDMatrix(bimf)
#  
#  
#  cat(paste("linking phenotype file", phenof, "\n"))
#  
#  bg <- as.BGData(bedFiles, alternatePhenotypeFile = paste0(phenof))
#  	
#  ## CHECKING ALL INPUT FILES AGAIN:
#  
#  pheno_dat <- pheno(bg)
#  geno_dat <- geno(bg)
#  
#  
#  if (!is.null(covarNames)){
#  
#  if (sum(grepl("sex|SEX|Sex", covarNames)) > 0){
#  
#  	SEX_cov <- pheno_dat[,names(pheno_dat) %in% covarNames][grepl("sex|SEX|Sex", covarNames)][,1]
#  	covarNames_use <- covarNames[!grepl("sex|SEX|Sex", covarNames)]
#  	SEX_cov_PLINK <- ifelse(SEX_cov==0, 2, SEX_cov)
#  
#  cat(paste("Covariates include", covarNames, " from covariate/pheno file \n"))
#  
#  } else {
#  
#  	SEX_cov <- pheno_dat$SEX
#  	SEX_cov_PLINK <- ifelse(SEX_cov==0, 2, SEX_cov)
#  	covarNames_use <- covarNames
#  
#  cat(paste("Covariates did not include SEX, taking SEX from .fam file\n"))
#  
#  }
#  
#  cat(paste("Writing results to output", out, "\n"))
#  
#  ## writing results by chunks of 50 SNPs to avoid loss in interruption
#  
#  iter <- round(dim(geno_dat)[2]/chunk_size)
#  
#  if (xchr) {
#  write.table(gJLS2(GENO = geno_dat[,(1):(chunk_size)], Y = pheno_dat[,names(pheno_dat) %in% phenoNames[1]], SEX = SEX_cov_PLINK, COVAR = pheno_dat[,names(pheno_dat) %in% covarNames_use], Xchr=xchr), file = out, col.names=T, row.names=F, quote=F, sep="\t")	
#  } else {
#  write.table(gJLS2(GENO = geno_dat[,(1):(chunk_size)], Y = pheno_dat[,names(pheno_dat) %in% phenoNames[1]], COVAR = pheno_dat[,names(pheno_dat) %in% covarNames], Xchr=xchr), file = out, col.names=T, row.names=F, quote=F, sep="\t")		
#  }
#  
#  
#  for (j in 2:iter){
#  	
#  	cat(paste("Running the", j, "th chunk", "\n"))
#  
#  if (j == iter){
#  
#  if (xchr) {
#  write.table(gJLS2(GENO = geno_dat[,(1 + chunk_size*(iter-1)):dim(geno_dat)[2]], Y = pheno_dat[,names(pheno_dat) %in% phenoNames[1]], SEX = SEX_cov_PLINK, COVAR = pheno_dat[,names(pheno_dat) %in% covarNames_use], Xchr=xchr), file = out, col.names=F, row.names=F, quote=F, append=TRUE, sep="\t")
#  } else {
#  write.table(gJLS2(GENO = geno_dat[,(1 + chunk_size*(iter-1)):dim(geno_dat)[2]], Y = pheno_dat[,names(pheno_dat) %in% phenoNames[1]], COVAR = pheno_dat[,names(pheno_dat) %in% covarNames], Xchr=xchr), file = out, col.names=F, row.names=F, quote=F, append=TRUE, sep="\t")	
#  }
#  
#  } else {	
#  if (xchr){
#  write.table(gJLS2(GENO = geno_dat[,(1 + chunk_size*(j-1)):(chunk_size*j)], Y = pheno_dat[,names(pheno_dat) %in% phenoNames[1]], SEX = SEX_cov_PLINK, COVAR = pheno_dat[,names(pheno_dat) %in% covarNames_use], Xchr=xchr), file = out, col.names=F, row.names=F, quote=F, append=TRUE, sep="\t")
#  } else {
#  write.table(gJLS2(GENO = geno_dat[,(1 + chunk_size*(j-1)):(chunk_size*j)], Y = pheno_dat[,names(pheno_dat) %in% phenoNames[1]], COVAR = pheno_dat[,names(pheno_dat) %in% covarNames], Xchr=xchr), file = out, col.names=F, row.names=F, quote=F, append=TRUE, sep="\t")	
#  }
#  }
#  }
#  
#  # locReg(GENO = geno_dat[,1:3], Y = pheno_dat[,names(pheno_dat) %in% phenoNames[1]], SEX = SEX_cov_PLINK, COVAR = pheno_dat[,names(pheno_dat) %in% covarNames_use], Xchr=xchr)
#  # scaleReg(GENO = geno_dat[,1:3], Y = pheno_dat[,names(pheno_dat) %in% phenoNames[1]], SEX = SEX_cov_PLINK, COVAR = pheno_dat[,names(pheno_dat) %in% covarNames_use], Xchr=xchr)
#  #gJLS2(GENO = geno_dat[,1:dim(geno_dat)[2]], Y = pheno_dat[,names(pheno_dat) %in% phenoNames[1]], SEX = SEX_cov_PLINK,#COVAR = pheno_dat[,names(pheno_dat) %in% covarNames_use], Xchr=xchr)
#  
#  } else {
#  
#  cat(paste("Writing results to output", out, "\n"))
#  	
#  
#  iter <- round(dim(geno_dat)[2]/chunk_size)
#  
#  if (xchr){
#  write.table(gJLS2(GENO = geno_dat[,(1):(chunk_size)], Y = pheno_dat[,names(pheno_dat) %in% phenoNames[1]], SEX = SEX_cov_PLINK,  Xchr=xchr), file = out, col.names=T, row.names=F, quote=F, sep="\t")
#  } else {
#  write.table(gJLS2(GENO = geno_dat[,(1):(chunk_size)], Y = pheno_dat[,names(pheno_dat) %in% phenoNames[1]],  Xchr=xchr), file = out, col.names=T, row.names=F, quote=F, sep="\t")	
#  }
#  
#  for (j in 2:iter){
#  
#  	cat(paste("Running the", j, "th chunk", "\n"))
#  
#  if (j == iter){
#  
#  if (xchr){	
#  write.table(gJLS2(GENO = geno_dat[,(1 + chunk_size*(iter-1)):dim(geno_dat)[2]], Y = pheno_dat[,names(pheno_dat) %in% phenoNames[1]], SEX = SEX_cov_PLINK,  Xchr=xchr, head=FALSE), file = out, col.names=F, row.names=F, quote=F, append=TRUE, sep="\t")
#  } else {
#  write.table(gJLS2(GENO = geno_dat[,(1 + chunk_size*(iter-1)):dim(geno_dat)[2]], Y = pheno_dat[,names(pheno_dat) %in% phenoNames[1]], Xchr=xchr, head=FALSE), file = out, col.names=F, row.names=F, quote=F, append=TRUE, sep="\t")	
#  }
#  
#  
#  } else {	
#  
#  if (xchr){
#  write.table(gJLS2(GENO = geno_dat[,(1 + chunk_size*(j-1)):(chunk_size*j)], Y = pheno_dat[,names(pheno_dat) %in% phenoNames[1]], SEX = SEX_cov_PLINK,  Xchr=xchr, head= FALSE), file = out, col.names=F, row.names=F, quote=F, append=TRUE, sep="\t")
#  } else {
#  write.table(gJLS2(GENO = geno_dat[,(1 + chunk_size*(j-1)):(chunk_size*j)], Y = pheno_dat[,names(pheno_dat) %in% phenoNames[1]], Xchr=xchr, head= FALSE), file = out, col.names=F, row.names=F, quote=F, append=TRUE, sep="\t")
#  }
#  }
#  }
#  }

## -----------------------------------------------------------------------------
Rplink <- function(PHENO,GENO,CLUSTER,COVAR){
 require(gJLS2)
 
  f1 <- function(s) 
       {    
      r <-  gJLS2(GENO=s, SEX=ifelse(COVAR[,1]==0, 2, COVAR[,1]), Y=PHENO, COVAR=COVAR[,-1], Xchr=TRUE)
      rr <- as.numeric(r[3:5])
      c( length(rr) , rr )
       }
      apply( GENO , 2 , f1 )
}

