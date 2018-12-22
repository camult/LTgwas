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


 Returns SNP effects, SNP variance and p-values using t test, Z score and Wald test.


## Author


 Fernando Brito Lopes


## References


 VanRaden, Paul M. Efficient methods to compute genomic predictions. Journal of dairy science 91.11 (2008): 4414-4423.
 
 Gualdron Duarte JL, Cantet RJ, Bates RO, Ernst CW, Raney NE, Steibel JP. Rapid screening for phenotype-genotype associations by linear transformations of genomic evaluations. BMC Bioinformatics (2014). 19;15:246. doi: 10.1186/1471-2105-15-246.


## Examples

```r 
 
 # GWAS <- ltgwas(file = "filename",
 #                varG=10,
 #                varE=20,
 #                gebv=gebv)
 
 ``` 

