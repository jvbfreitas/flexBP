\name{flexBP}
\alias{flexBP}
\title{
  Fit flexible quasi-beta prime regression models
}
\description{
  The flexBP function fits flexible quasi-beta prime regression models (FREITAS et al., 2024).
}
\usage{
  flexBP(formula, data, id, tol = 0.001,
        maxiter = 25, compm, linkmu = "log", tuning = 1)
}
\arguments{
  \item{formula}{
    an object of class "formula" (or one that can be coerced to that class) related to the mean predictor: a expression of the form "response ~ predictors". See the documentation of \code{\link{glm}} and \code{\link{formula}} for details.
  }
  \item{data}{
    data frame containing the variables in the formula.
  }
  \item{tol}{
    the tolerance used in the fitting algorithm.
  }
  \item{maxiter}{
    the maximum number of iterations.
  }
  \item{compm}{
    a list of structured matrices to be used on the estimation of the covariance. For details see \code{\link{mcglm}}.
  }
  \item{link}{
    The link function for mean model. The following are permitted: "log" and "identity".
  }
  \item{tuning}{
    The value of the tuning parameter, for more detail see FREITAS et. al. (2024)
  }
}
\details{

}
\value{
  An object of class "flexBP" representing the fit.
}
\references{
  FREITAS, João; NOBRE, Juvêncio; AZEVEDO, Caio. Flexible quasi-beta prime regression models for dependent continuous positive data. Statistics and Its Interface, v. 17, n. 4, p. 715-731, 2024.
}
\author{
  João de Freitas, jvbfreitas@ime.unicamp.br.
}

\examples{

}
