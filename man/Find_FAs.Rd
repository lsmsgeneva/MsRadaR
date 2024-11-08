% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Processing_Functions.R
\name{Find_FAs}
\alias{Find_FAs}
\title{Identify Potential Fatty Acids (FAs) from Monoglyceride Data}
\usage{
Find_FAs(Data, Threeshold = 10, Oxygen = 3, median_Factor = 4)
}
\arguments{
\item{Data}{A data frame containing the monoglyceride data, including columns for
\code{Height}, \code{C} (number of carbon atoms), and \code{O} (number of oxygen atoms).}

\item{Threeshold}{Numeric. The minimum height value for a peak to be considered
(default: 10).}

\item{Oxygen}{Numeric vector. The specific oxygen counts to filter the data on
(default: 3).}

\item{median_Factor}{Numeric. A factor used to scale the median when determining
threshold conditions (default: 4).}
}
\value{
A list where each element corresponds to the specified oxygen count and
contains carbon counts that exceed the calculated threshold. Returns "No Peak"
if no data meets the criteria for a given oxygen count.
}
\description{
This function identifies potential fatty acids from monoglyceride (MG) data based
on specified thresholds and oxygen counts. It filters data based on height and
oxygen content, summarizes the height by carbon count, and determines which carbon
counts exceed a calculated threshold.
}
\examples{
# Example usage of Find_FAs function
identified_fas <- Find_FAs(Data = mg_data, Oxygen = c(1, 2, 3))

}
