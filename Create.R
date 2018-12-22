library("Rcpp")
library("devtools")
library("roxygen2")

setwd("C:\\Users\\Camult\\Documents\\GitHub\\LTgwas")
compileAttributes(verbose=TRUE)
roxygenise()



library("Rcpp")
library("devtools")
library("roxygen2")
setwd("C:\\Users\\Camult\\Documents\\GitHub")
system("R CMD check LTgwas --no-manual --no-examples")
system("R CMD build LTgwas")
system("R CMD INSTALL LTgwas_1.0.0.tar.gz --no-multiarch")

library("Rd2md")
setwd("C:\\Users\\Camult\\Documents\\GitHub")
Rd2markdown(rdfile="C:\\Users\\Camult\\Documents\\GitHub\\LTgwas\\man\\ltgwas.Rd",
            outfile="C:\\Users\\Camult\\Documents\\GitHub\\LTgwas\\README.md")

