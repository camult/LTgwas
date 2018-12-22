#' @title Association Teste using Linear Transformation of Genomic Evaluations
#' 
#' @description Standardized test of association where genomic breeding values from a mixed model including fixed effects and a genomic information are linearly transformed to estimate marker effects and pvalues.
#' 
#' @param file name of all plink files (.bed, bim and .fam)
#' @param varG genetic variance 
#' @param varE residual variance
#' @param gebv genomic estimated breeding value
#' 
#' @return Returns SNP effects and SNP variance
#' 
#' @author Fernando Brito Lopes
#' 
#' @references VanRaden, Paul M. Efficient methods to compute genomic predictions. Journal of dairy science 91.11 (2008): 4414-4423.
#' 
#' @references Gualdron Duarte JL, Cantet RJ, Bates RO, Ernst CW, Raney NE, Steibel JP. Rapid screening for phenotype-genotype associations by linear transformations of genomic evaluations. BMC Bioinformatics (2014). 19;15:246. doi: 10.1186/1471-2105-15-246.
#' 
#' @examples
#'
#' # GWAS <- ltgwas(file = "filename",
#' #                varG=10,
#' #                varE=20,
#' #                gebv=gebv)
#' 
#' @export ltgwas
#' @importFrom stats pchisq pnorm pt
#' @useDynLib LTgwas, .registration = TRUE
ltgwas <- function(file=NULL,varG=NULL,varE=NULL,gebv=NULL) {
     bim <- read.table(paste0(file,".bim"), header=FALSE)
     fam <- read.table(paste0(file,".fam"), header=FALSE)
     fnRAW <- paste0(file,".raw")
     if(file.exists(fnRAW)) file.remove(fnRAW)
     if(file.exists("fort.6")) file.remove("fort.6")
     n <- nrow(fam)
     m <- nrow(bim)
     cls <- 1:m
     rws <- 1:n
     keep <- rep(TRUE,m)
     nr <- length(rws)
     nc <- length(cls)
     nbytes <- ceiling(n/4)
     append <- 0
     fnBED <- paste0(file,".bed")  
     if(.Platform$OS.type=="windows") fnRAW <- tolower(gsub("/","\\",fnRAW,fixed=T))    
     bedList <- .Fortran("bed2raw",
                         m = as.integer(m),
                         cls = as.integer(keep),
                         nbytes = as.integer(nbytes),
                         append = as.integer(append),
                         fnBED = as.character(fnBED),
                         fnRAW = as.character(fnRAW),
                         PACKAGE = 'LTgwas')
     ltgList <- .Fortran("ltgebed",
                         n = as.integer(n),
                         m = as.integer(m),
                         rws = as.integer(rws),
                         cls = as.integer(cls),
                         nbytes = as.integer(nbytes),
                         fnRAW = as.character(fnRAW),
                         varG = as.double(varG),
                         varE = as.double(varE),
                         gebv = as.double(gebv),
                         SNPeff = as.double(rep(0,m)),
                         SNPvar = as.double(rep(0,m)),
                         nI = as.matrix(diag(n) - 1/(n-1)),
                         tol = as.double(1.8e-03),
                         PACKAGE = 'LTgwas')
     SNPeff <- ltgList$SNPeff
     SNPvar <- ltgList$SNPvar
     std_SNP <- SNPeff/sqrt(SNPvar)
     out <- list(Name=bim[,2],
                 Chr=bim[,1],
                 Pos=bim[,4],
                 SNPeff=SNPeff,
                 SNPvar=SNPvar,
                 p.Zscore=2*pnorm(q=abs(std_SNP), mean=0, sd=1, lower.tail=F, log.p=F),
                 p.ttest=2*pt(q=abs(std_SNP), df=n-1, lower.tail=F, log.p=F),
                 p.Wald=pchisq(q=(std_SNP^2), df=1, lower.tail=F, log.p=F))
     if(file.exists(fnRAW)) file.remove(fnRAW)
     if(file.exists(fnRAW)) file.remove("fort.6")
     return(out)
}
