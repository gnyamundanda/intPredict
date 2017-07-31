
# intPredict

An R package for gene selection using a combination of different gene selection and class prediction methods.

# INSTALLING intPredict.
To install indtPredict you require devtools package. Here are two lines of codes to install intPredict from our syspremed account on github.

require(devtools)

install_github("syspremed/intPredict") 

# Running the example in the paper. 
The following lines of R code demonstrate how to run intPredict using the example in the paper.

require(intPredict) 

data(edata)                                     # CRCassigner training data with 387 CRC samples

data(labs)                                      # Corresponding CRCassigner classes for 387 samples

seqp=c(seq(2,44,2),47)                          # A sequence for top number of genes to be selected 

fac=2/3                                         # A fraction of samples to be used as training samples

intpred(edata, labs, seqp, fac, seed=166)
