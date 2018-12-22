# `ltgwas`: Association Teste using Linear Transformation of Genomic Evaluations

## Description


 Standardized test of association where genomic breeding values from a mixed model including fixed effects and a genomic information are linearly transformed to estimate marker effects and pvalues.


## Usage

```r
ltgwas(file = NULL, varG = NULL, varE = NULL, gebv = NULL)
```


## Arguments

Argument      |Description
------------- |----------------
```file```     |     name of all plink files (.bed, bim and .fam)
```varG```     |     genetic variance
```varE```     |     residual variance
```gebv```     |     genomic estimated breeding value

## Value


 Returns SNP effects and SNP variance


## Author


 Fernando Brito Lopes


## Examples

```r 
 
 # GWAS <- ltgwas(file = "filename,
 #                varG=10,
 #                varE=20,
 #                gebv=gebv)
 
 ``` 

