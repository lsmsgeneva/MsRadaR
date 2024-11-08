% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Load_Examples.R
\name{Load_Example_Info}
\alias{Load_Example_Info}
\title{Load Example Spectra Information}
\usage{
Load_Example_Info()
}
\value{
This function does not return any value. It sets the global variable
\code{Example_Info} with the path to the folder of the example spectra and sample
information excel file.
}
\description{
This function loads all required information for batch processing from the
directory of the MsRadaR package. This includes the folder location and file
location of the "Spectra_Informations.xlsx" information file.
}
\examples{
# Load example spectra information
Load_Example_Info()

}
