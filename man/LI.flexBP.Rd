\name{LI.flexBP}
\alias{LI.flexBP}
\title{
  Compute local influence measures for flexBP model
}
\description{
  This function can be used to compute the local influence measures under certain perturbations for flexible quasi-beta prime regression models based in FREITAS et. al. (2024).
}
\usage{
  LI.flexBP(model, pert = c("case-weight", "response"))
}
\arguments{
  \item{model}{
    an R object returned by flexBP.
  }
  \item{pert}{
    a character string specifying the perturbation scheme. The following are permitted: "case-weight" and "response".
  }
}
\details{
}
\value{
  A vector is returned with the local influence values.
}
\references{
  FREITAS, João; NOBRE, Juvêncio; AZEVEDO, Caio. Flexible quasi-beta prime regression models for dependent continuous positive data. Statistics and Its Interface, v. 17, n. 4, p. 715-731, 2024.
}
\author{
  João de Freitas, jvbfreitas@ime.unicamp.br.
}
\examples{

}
