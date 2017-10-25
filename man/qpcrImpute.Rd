\name{qpcrImpute}
\alias{qpcrImpute}
\title{Impute Non-detects in qPCR data}
\description{
  This function models the missing data mechanism and uses an EM
  algorithm to impute the non-detect values in qPCR data. 
}
\usage{
qpcrImpute(object, dj=NULL, pyfit=NULL, groupVars=NULL, 
          batch=NULL, tol=1, iterMax=100, 
          outform=c("Single","Param","Multy"), 
          vary_fit=TRUE, vary_model=TRUE, add_noise=TRUE, 
          formula=NULL, numsam=5, 
          linkglm = c("logit", "probit", "cloglog"))
}
\arguments{
  \item{object}{a qPCRset}
  \item{dj}{normalization values. If NULL, features with "control"
    in featureType(object) are used to normalize the data. If no control
    features are found, the data are not normalized.}
  \item{pyfit}{initial estimate of the relationship between the
    probability of a non-detect and average expression. If NULL, this
    relationship is estimated from the data.} 
  \item{groupVars}{which columns in pData(object) should be used to
    determine replicate samples. If NULL, all columns are used.} 
  \item{batch}{amatrix with control samples for each batch, if NULL, batch effect is not taken into account.}
  \item{tol}{likelihood convergence criterion of the EM algorithm.}
  
  \item{iterMax}{maximimum number of iterations of the EM algorithm.}
  \item{outform}{the form of the output requested.If "Single" performes a single imputation of missing values. If "Param" returnes estimated model parameters: mean and variance. If "Multy" performes a multiple imputation of missing values, and creats multiple data sets with imputed values.}
  
  \item{vary_fit}{if outform="Multy", includes the model uncertainty due to the logit of the probability of being missing. The default value is "TRUE".}
  \item{vary_model}{if outform="Multy", includes the model uncertainty due to the estimating mean of the data. The default value is "TRUE".}
  \item{add_noise}{if outform="Multy", introduses the variance component due to the random noise. The default value is "TRUE".}
  
  \item{formula}{specifies the model.}
  \item{numsam}{number of the datasets to be created if outform="Multy". The default value is 5.}
  \item{linkglm}{a link used for estimation of the missing data mechanism.}

}
\value{
  The function returns a qPCRset object with non-detects replaced by their imputed values. 
}
\author{Valeriia Sherina}
\examples{
  data(sagmb2011)
  tst <- qpcrImpute(sagmb2011, groupVars="sampleType",
                    outform=c("Single"), batch=NULL, linkglm = c("logit"))
}
\keyword{manip}

