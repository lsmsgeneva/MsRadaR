% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Processing_Functions.R
\name{Annotate_Fragments}
\alias{Annotate_Fragments}
\title{Annotate Fragments with Elemental Formulas}
\usage{
Annotate_Fragments(
  Data,
  mDA_Error = 20,
  Elements = Elements_Allowed,
  minForm = minFormula,
  maxForm = Precursor_DF$Elemental_Formula,
  maxDBE = Precursor_DF$DBE + 2,
  minDBE = -0.5
)
}
\arguments{
\item{Data}{A data frame with at least a column for mass/charge values (\code{Mass.Charge})
to be annotated with elemental formulas.}

\item{mDA_Error}{Numeric. Maximum allowed error for mass decomposition in mDa (default: 20).}

\item{Elements}{A vector of allowed elements for formula generation.}

\item{minForm}{Character. Minimum formula constraint for fragment annotation.}

\item{maxForm}{Character. Maximum formula constraint based on precursor formula (default: Precursor_DF).}

\item{maxDBE}{Numeric. Maximum allowed double bond equivalent for annotation.}

\item{minDBE}{Numeric. Minimum allowed double bond equivalent for annotation (default: -0.5).}
}
\value{
A data frame of annotated fragments with mass/charge, elemental formula,
error, DBE, and additional calculated values based on fragmentation models.
}
\description{
Annotates MS/MS spectra with possible elemental formulas for each fragment
mass by calculating and filtering based on user-defined parameters.
}
\examples{
# Example of annotating data fragments with default parameters
annotated_fragments <- Annotate_Fragments(Data = spectra_data)

}
