% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Create_Summary.R
\name{Create_Summary}
\alias{Create_Summary}
\title{Create a Summary Plot for Double Bonds from Batch Results Data}
\usage{
Create_Summary(Data, Expected_DB = NA, Plot.Order = "Order")
}
\arguments{
\item{Data}{A data frame containing the batch results with necessary columns
and structures, including annotated loss data.}

\item{Expected_DB}{An optional parameter indicating the number of double
bonds column within the info.file. If NA, it calculates expected double
bonds based on the data (default is NA).}

\item{Plot.Order}{A character string specifying the column within the
info.file that indicates the desired plot order (default is "Order").}
}
\value{
A ggplot object representing the summary plot of double bonds.
}
\description{
This function generates a summary plot based on the input batch results
data, focusing on double bonds. It allows for optional parameters to
specify expected double bonds and the desired plot order.
}
\examples{
# Example usage of Create_Summary function
summary_plot <- Create_Summary(Batch_Results, Expected_DB = "Max_DB_Col")

}
