# MsRadaR, A Tool for De-Novo Structural Elucidation of Lipids from MS/MS Spectra (EDP-CID, EAD)
# Copyright (C) 2024 P. Mueller, G. Hopfgartner
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

#===============================================================================
#Data Loading for Single File Processing (Example)
#===============================================================================

#' Load Example Spectra Path
#'
#' This function generates the path to a specified example spectrum file
#' from the "Example Files" directory of the MsRadaR package based on
#' the provided sample name and file type.
#'
#' @param Sample A character string specifying the sample name. This is used
#'   to construct the filename for the example spectrum.
#' @param file.type A character string specifying the file type. Default is ".txt".
#'
#' @return A character string representing the path to the specified example
#'   spectrum file.
#'
#' @examples
#' # Load the path for the  sample "LLL.txt"
#' Example_Spectra_Path <- Load_Example_Path(Sample = "LLL")
#'
#' @export
Load_Example_Path <- function(Sample, file.type = ".txt") {

  Spectra <- list.files(system.file("Example Files", package = "MsRadaR"), pattern = ".txt")

  Example_File <- paste(Sample, file.type, sep = "")

  system.file("Example Files", Example_File, package = "MsRadaR")
}

#===============================================================================
##Data Loading for Batch Processing (Example)
#===============================================================================

#' Load Example Spectra Information
#'
#' This function loads all required information for batch processing from the
#' directory of the MsRadaR package. This includes the folder location and file
#' location of the "Spectra_Informations.xlsx" information file.
#'
#' @return This function does not return any value. It sets the global variable
#'   `Example_Info` with the path to the folder of the example spectra and sample
#' information excel file.
#' @examples
#' # Load example spectra information
#' Load_Example_Info()
#'
#' @export
Load_Example_Info <- function() {

  Example_Folder <<- system.file("Example Files", package = "MsRadaR")

  Example_Info <<- system.file("Example Files", "Spectra_Informations.xlsx", package = "MsRadaR")
}

