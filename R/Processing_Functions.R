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


#' Load Default Parameters for Elemental Formula Calculation
#'
#' Loads default parameters for elemental formula calculation and data loading.
#' These settings can be adjusted according to user-defined preferences.
#'
#' @param Threshold Numeric. Intensity threshold for data loading (default: 2).
#' @param mDaError Numeric. Allowed mass error in mDa (default: 20).
#' @param Kendricks_Defect Numeric. Kendrick mass defect (default: 14 / 14.01565).
#' @param minFormula_Lipid Character. Minimum elemental formula for lipids (default: "C3H8O2").
#' @param min_formula Character. Minimum elemental formula (default: "C0H1").
#' @param max_Formula Character. Maximum elemental formula (default: "C200H800O6").
#' @param max_Formula_FA Character. Maximum elemental formula for fatty acids (default: "C20H40O4").
#' @param Allowed_Elements Character. String of allowed elements (default: "CHO").
#'
#' @return Sets global variables for threshold, allowed mass error, Kendrick defect,
#' and elemental formula constraints based on specified inputs.
#'
#' @examples
#' Load_Parameters()
#' Load_Parameters(Threshold = 5, mDaError = 15)
#'
#' @export
Load_Parameters <- function(Threshold = 2,
                            mDaError = 20,
                            Kendricks_Defect = 14/14.01565,
                            minFormula_Lipid = "C3H8O2",
                            min_formula = "C0H1",
                            max_Formula = "C200H800O6",
                            max_Formula_FA = "C20H40O4",
                            Allowed_Elements = "CHO") {

  Thres <<- Threshold
  mDA_Error <<- mDaError
  Kendrick_CH2 <<- Kendricks_Defect
  minFormula_TG <<- minFormula_Lipid
  minFormula <<- min_formula
  maxFormula <<- max_Formula
  maxFormula_FA <<- max_Formula_FA
  Elements <- unlist(strsplit(Allowed_Elements, "(?=[A-Z])", perl = TRUE))
  Elements_Allowed <<- initializeElements(Elements)
}

#===============================================================================
# Load_Data
#===============================================================================
#' Load Data and Perform Elemental Formula Assignment
#'
#' Loads data from a specified file and performs elemental formula assignment.
#' This function allows for either batch processing or single-spectrum processing.
#'
#' @param Dat.Name Character. Name of the file containing the data.
#' @param Precursor.Mass Numeric or logical. If a numeric value, it sets the precursor mass;
#'   if \code{FALSE}, uses the mass with the highest intensity (default: \code{FALSE}).
#' @param mz Character. Name of the column representing mass/charge values (default: "Mass.Charge").
#' @param Intensity Character. Name of the column representing intensity values (default: "Height").
#'
#' @return A data frame with the processed data, including relative intensity, precursor formula, and mass loss.
#'   If no valid formula match is found, an error message is displayed.
#'
#' @examples
#' # Assuming data is available in a file named "sample_data.txt"
#' Load_Data("sample_data.txt")
#'
#' @export
Load_Data <- function(Dat.Name = "", Precursor.Mass = FALSE, mz = "Mass.Charge", Intensity = "Height") {
  suppressMessages(reset())# required loading for CHNOSZ. CHNOSZ should be replaced.
  Data.Name <- paste(Dat.Name, sep="")
  Data <- read.delim(Data.Name)
  Data <- Data[, colnames(Data) %in% c(mz, Intensity)]

  # Set precursor mass or use the highest intensity peak
  if (is.numeric(Precursor.Mass)) {
    Prec_Mass <- Precursor.Mass
  } else {
    Prec_Mass <- Data$Mass.Charge[which.max(Data$Height)]
  }

  # Set intensity threshold if not globally defined
  if (Thres == "") {
    Thres <- Find_Threeshold(data = Data)
  }

  # Filter data and prepare for further processing
  Data <- Data[Data$Height > Thres,]
  Data$Mass.Charge <- round(Data$Mass.Charge, 4)
  Data$Relative_Intensity <- rescale_max(Data$Height) * 100
  Data <- Data[Data$Mass.Charge < (Prec_Mass + 0.5),]
  Data$loss <- round(Data$Mass.Charge - Prec_Mass, 4)

  # Elemental formula decomposition of the precursor mass
  # Decompose precursor mass into possible elemental formulas
  Prec_Formula <- decomposeMass(abs(Prec_Mass),
                                mzabs = mDA_Error / 1000,
                                maxisotopes = 1,
                                elements = Elements_Allowed,
                                minElements = minFormula_TG,
                                maxElements = maxFormula)

  # Filter out formulas with invalid double bond equivalent (DBE) values
  Precursor <- which(Prec_Formula$DBE >= 0)
  Precursor <- Precursor[1]  # Select first valid formula match

  # Process if a valid formula match is found
  if (!is.na(Precursor)) {
    Precursor_DF <- c(round(Prec_Formula$exactmass[Precursor], 4),
                      Prec_Formula$formula[Precursor],
                      round((Prec_Mass - Prec_Formula$exactmass[Precursor]), 4),
                      Prec_Formula$DBE[Precursor])

    # Assign column names to precursor formula data
    names(Precursor_DF) <- c("exactmass", "Elemental_Formula", "Error", "DBE")

    # Convert precursor formula data to a data frame and expand to show elemental composition
    Precursor_DF2 <- as.data.frame(t(Precursor_DF))
    Precursor_DF <- cbind(Precursor_DF2, as.data.frame(t(makeup(Precursor_DF2$Elemental_Formula))))

    # Ensure DBE is numeric for compatibility
    Precursor_DF$DBE <- as.numeric(Precursor_DF$DBE)

    # Assign precursor data frame globally
    Precursor_DF <<- Precursor_DF

    # Return processed data
    Data
  } else {
    # Print error message if no valid formula match is found
    print(paste("An error occurred while processing the following data file:", Dat.Name))
    print("Precursor mass does not match elemental formula assignment criteria")
  }
}

#===============================================================================
# Kendrick Transformation
#===============================================================================
#' Kendrick Mass Defect Transformation
#'
#' Performs the Kendrick mass defect transformation on the input spectra data
#' and calculates the relative intensities of fragment ions.
#'
#' @param Kendricks_Data A data frame containing spectra data with at least
#'   columns for mass/charge values and intensity.
#' @param Kendrick_Defect Numeric. Kendrick mass defect, usually CH2 group (default: 14 / 14.01565).
#' @param Mass.Charge Character. Column name for mass/charge values in `Kendricks_Data`.
#'
#' @return A data frame with the transformed Kendrick mass defect values,
#'   relative intensities, and other calculated columns for further analysis.
#'
#' @examples
#' # Example of transforming data with default Kendrick defect for CH2 group
#' transformed_data <- Kendricks_Transformation(Kendricks_Data = your_data)
#'
#' @export
Kendricks_Transformation <- function(Kendricks_Data, Kendrick_Defect = 14 / 14.01565, Mass.Charge = Mass.Charge) {

  # Round m/z values to nearest integer
  Kendricks_Data$round <- as.integer(Kendricks_Data$Mass.Charge)

  # Calculate mass defect
  Kendricks_Data$defect <- Kendricks_Data$Mass.Charge - Kendricks_Data$round

  # Calculate Kendrick mass
  Kendricks_Data$Ke_CH2 <- Kendricks_Data$Mass.Charge * Kendrick_Defect

  # Calculate Kendrick mass defect
  Kendricks_Data$KDM_CH2 <- Kendricks_Data$Ke_CH2 - as.integer(Kendricks_Data$Ke_CH2)

  # Round Kendrick mass defect to 3 decimal places
  Kendricks_Data$KDM_CH2 <- round(Kendricks_Data$KDM_CH2, 3)

  # Scale intensity values (0-100%)
  Kendricks_Data$Rel.Int <- rescale_max(Kendricks_Data$Height) * 100

  # Add an alpha threshold column for further processing
  Kendricks_Data$alpha_thres <- NA

  # Return transformed data
  Kendricks_Data
}

#===============================================================================
# Annotate Fragments
#===============================================================================
#' Annotate Fragments with Elemental Formulas
#'
#' Annotates MS/MS spectra with possible elemental formulas for each fragment
#' mass by calculating and filtering based on user-defined parameters.
#'
#' @param Data A data frame with at least a column for mass/charge values (`Mass.Charge`)
#'   to be annotated with elemental formulas.
#' @param mDA_Error Numeric. Maximum allowed error for mass decomposition in mDa (default: 20).
#' @param Elements A vector of allowed elements for formula generation.
#' @param minForm Character. Minimum formula constraint for fragment annotation.
#' @param maxForm Character. Maximum formula constraint based on precursor formula (default: Precursor_DF).
#' @param maxDBE Numeric. Maximum allowed double bond equivalent for annotation.
#' @param minDBE Numeric. Minimum allowed double bond equivalent for annotation (default: -0.5).
#'
#' @return A data frame of annotated fragments with mass/charge, elemental formula,
#'   error, DBE, and additional calculated values based on fragmentation models.
#'
#' @examples
#' # Example of annotating data fragments with default parameters
#' annotated_fragments <- Annotate_Fragments(Data = spectra_data)
#'
#' @export
Annotate_Fragments <- function(Data,
                               mDA_Error = 20,
                               Elements = Elements_Allowed,
                               minForm = minFormula,
                               maxForm = Precursor_DF$Elemental_Formula,
                               maxDBE = Precursor_DF$DBE + 2,
                               minDBE = -0.5) {

  suppressMessages(reset())# required loading for CHNOSZ. CHNOSZ should be replaced.
  # Set maximum formula constraints if undefined
  if (is.na(maxForm)) {
    maxForm <- maxFormula
    maxDBE <- 50
  }

  # Calculate elemental formulas for each fragment based on exact mass
  Formula_List <- lapply(Data$Mass.Charge, function(z) {

    # Decompose fragment mass to obtain possible elemental formulas
    Decomposed <- decomposeMass(abs(z),
                                mzabs = mDA_Error / 1000,
                                maxisotopes = 1,
                                elements = Elements,
                                minElements = minForm,
                                maxElements = maxForm)

    # Remove the last element in the list (usually metadata)
    Decomposed <- Decomposed[-length(Decomposed)]
    Decomposed <- as.data.frame(Decomposed)

    # Remove elemental formulas with double bond equivalent values above the maximum threshold
    if (min(Decomposed$DBE) <= maxDBE) {
      Decomposed <- Decomposed[Decomposed$DBE <= maxDBE,]
    }

    # Remove elemental formulas with double bond equivalent values below the minimum threshold
    if (max(Decomposed$DBE) >= minDBE) {
      Decomposed <- Decomposed[Decomposed$DBE >= -1,]
    }

    # Annotate fragment with elemental formulas if a valid elemental formula has been calculated
    if (nrow(Decomposed) > 0) {
      Fragment_Error <- Decomposed$exactmass[1] - abs(z)
      data.frame(Mass.Charge = z, Elemental_Formula = Decomposed$formula[1],
                 Error = Fragment_Error, DBE = Decomposed$DBE[1])
    } else {
      # Assign NA if no valid elemental formula has been calculated
      data.frame(Mass.Charge = z, Elemental_Formula = NA, Error = NA, DBE = NA)
    }
  })

  # Combine all annotated fragments into a single data frame
  Formula_List2 <- do.call(rbind, Formula_List)

  # Remove entries which were not assigned any elemental formula
  Formula_List3 <- Formula_List2[!is.na(Formula_List2$Elemental_Formula),]

  # Merge elemental formulas to original input data
  Formula_List3 <- merge(Data, Formula_List2, by = "Mass.Charge")

  # Remove empty rows
  Formula_List4 <- Formula_List3[!is.na(Formula_List3$Elemental_Formula),]

  # Calculate C, H, O content for each elemental formula
  Formula_List5 <- cbind(Formula_List4, do.call(rbind, makeup(Formula_List4$Elemental_Formula, count.zero = TRUE)))

  #=============================================================================
  # Calculate CHx difference from expected values for common fragmentation types
  #=============================================================================

  # Loss from the CH3 side end of a fatty acid
  Formula_List5$L_R <- (Formula_List5$C * 2 + 1) - Formula_List5$H

  # Losses from the carboxy (CO2) side of a fatty acid
  Formula_List5$R_L <- ((Formula_List5$C - 1) * 2) - Formula_List5$H

  # Loss from both sides (possible mixed fragmentation)
  Formula_List5$RnL <- NA
  Formula_List5$RnL[Formula_List5$C > 1] <- ((Formula_List5$C[Formula_List5$C > 1] - 1) * 2 + 1) -
    Formula_List5$H[Formula_List5$C > 1]

  #=============================================================================

  # Remove duplicate elemental formula assignments and keep only the maximal intensity fragment ion
  Unify_Formulas(Formula_List5, Parameter = "Elemental_Formula")
}

#===============================================================================
# Unify Formulas
#===============================================================================
#' Unify Duplicate Elemental Formulas by Maximum Intensity
#'
#' This function identifies and unifies rows with duplicate elemental formulas
#' based on a specified parameter (default is "Elemental_Formula") by keeping
#' only the entry with the highest intensity for each duplicate.
#'
#' @param Data A data frame containing at least an intensity column (`Height`)
#'   and the column specified in `Parameter`.
#' @param Parameter Character. The name of the column used to identify duplicates
#'   (default: "Elemental_Formula").
#'
#' @return A data frame where each unique value in the specified column has
#'   only the row with the maximum intensity retained.
#'
#' @examples
#' # Example of unifying duplicate elemental formulas
#' unified_data <- Unify_Formulas(Data = annotated_data)
#'
#' @export
Unify_Formulas <- function(Data, Parameter = "Elemental_Formula") {

  # Identify duplicates based on the specified parameter (e.g., "Elemental_Formula")
  Check <- Data[duplicated(Data[Parameter]), Parameter]

  # For each duplicate elemental formula, keep only the entry with the highest intensity
  filtered <- lapply(Check, function(x) {

    # Subset data for rows with the current duplicate elemental formula
    Output <- Data[Data[, Parameter] == x, ]

    # Select the row with the maximum intensity ("Height") among duplicates
    Output[which.max(Output$Height), ]
  })

  # Combine the filtered rows (one per duplicate elemental formula) into a data frame
  filtered2 <- do.call(rbind, filtered)

  # Keep rows without duplicates in the specified parameter
  unique_rows <- Data[!Data[, Parameter] %in% Check, ]

  # Combine unique rows and filtered rows with maximum intensity for duplicates
  rbind(unique_rows, filtered2)
}


#===============================================================================
# Annotate Losses
#===============================================================================
#' Annotate Elemental Formula Losses in MS/MS Spectra
#'
#' This function annotates elemental formula losses of a given MS/MS spectrum
#' using user-defined limits for elemental formula annotations and filtering.
#'
#' @param Data A data frame containing MS/MS spectra data, including a column
#'   for mass loss values.
#' @param mDA_Error Numeric. The allowed mass error in milliDaltons (default: 20).
#' @param Elements Character vector specifying the elements to consider
#'   (default: Elements_Allowed).
#' @param minForm Character. Minimum elemental formula for decomposition
#'   (default: minFormula).
#' @param maxForm Character. Maximum elemental formula constraint
#'   (default: maxFormula_FA).
#' @param DBE_Filter Numeric. The minimum Double Bond Equivalent (DBE)
#'   threshold (default: -0.5).
#' @param Keep_C0 Logical. Whether to keep entries without carbon in the
#'   annotation (default: FALSE).
#'
#' @return A data frame containing annotated elemental formula losses, including
#'   derived information on Carbon position, series representation, intensity
#'   normalization, and calculated slopes.
#'
#' @examples
#' # Example usage of Annotate_Losses function
#' annotated_losses <- Annotate_Losses(Data = ms_data, mDA_Error = 20, Keep_C0 = TRUE)
#'
#' @export
Annotate_Losses <- function(Data,
                            mDA_Error = 20,
                            Elements = Elements_Allowed,
                            minForm = minFormula,
                            maxForm = maxFormula_FA,
                            DBE_Filter = -0.5,
                            Keep_C0 = F) {
  suppressMessages(reset())# required loading for CHNOSZ. CHNOSZ should be replaced.
  # If maxForm is NA, set it to maxFormula
  if (is.na(maxForm)) {
    maxForm <- maxFormula
  }

  # Decompose each loss into elemental formulas by looping through the data
  Loss_Information <- lapply(Data$loss, function(z) {

    # Decompose the mass loss into possible elemental formulas
    Decomposed <- decomposeMass(abs(z),
                                mzabs = mDA_Error / 1000,
                                maxisotopes = 1,
                                elements = Elements_Allowed,
                                minElements = minForm,
                                maxElements = maxForm)

    # Remove the last element of the decomposition result
    Decomposed <- Decomposed[-length(Decomposed)]

    # Convert the decomposed results to a data frame
    Decomposed <- as.data.frame(Decomposed)

    # Filter based on the Double Bond Equivalent (DBE) criteria
    if (max(Decomposed$DBE) >= DBE_Filter) {
      Decomposed <- Decomposed[Decomposed$DBE >= DBE_Filter - 0.5, ]
    }

    # If there are valid elemental formula losses, create a data frame
    if (nrow(Decomposed) > 0) {
      Fragment_Error <- Decomposed$exactmass[1] - abs(z)  # Calculate error
      data.frame(Loss_Mass = z, Loss_Formula = Decomposed$formula[1],
                 Loss_Error = Fragment_Error, DBE = Decomposed$DBE[1])
    } else {
      # Return NA values if no valid elemental formula was assigned
      data.frame(Loss_Mass = z, Loss_Formula = NA, Loss_Error = NA, DBE = NA)
    }
  })

  # Bind the list of elemental loss information to a data frame
  Loss_Information2 <- do.call(rbind, Loss_Information)

  # Combine elemental formula losses with the original input data
  Loss_Data <- merge(Data, Loss_Information2, by.x = "loss", by.y = "Loss_Mass")

  # Remove entries with unidentified elemental formulas
  Loss_Data2 <- Loss_Data[!is.na(Loss_Data$Loss_Formula), ]

  # Calculate the number of Carbon (C), Hydrogen (H), and Oxygen (O) atoms of each elemental formula loss
  Loss_Data3 <- cbind(Loss_Data2, do.call(rbind, (makeup(Loss_Data2$Loss_Formula, count.zero = T))))



  # Set Oxygen count to 0 if it does not exist in the data
  if (max(colnames(Loss_Data3) %in% "O") == 0) {
    Loss_Data3$O <- 0
  }

  # Define carbon position using the "n-x" nomenclature for lipids
  Loss_Data3$Carbon_Position <- paste("n-", Loss_Data3$C + 1, sep = "")

  # Convert Carbon_Position to a factor for ordered plotting
  Loss_Data3$Carbon_Position <- factor(Loss_Data3$Carbon_Position,
                                       levels = unique(Loss_Data3$Carbon_Position)[order(unique(Loss_Data3$C), decreasing = TRUE)])

  # Create a simplified series representation of the elemental formula
  Loss_Data3$Carbon_Series <- sub("O[0-9]+", "", Loss_Data3$Loss_Formula)
  Loss_Data3$Carbon_Series <- sub("O", "", Loss_Data3$Carbon_Series)

  # Remove losses without carbon unless Keep_C0 is TRUE (C0 = precursor)
  if (Keep_C0 == F) {
    Loss_Data4 <- Loss_Data3[!Loss_Data3$C == 0, ]
  } else {
    Loss_Data4 <- Loss_Data3
  }

  # Initialize a scaled intensity column
  Loss_Data4$scaled_intensity <- NA

  # Calculate the hydrogen difference based on the assumption of losing CH3(CH2)n
  Loss_Data4$delta_H <- (Loss_Data4$C * 2 + 1) - Loss_Data4$H

  # Create labels based on the hydrogen difference
  Loss_Data4$Label <- paste(Loss_Data4$delta_H, "H", sep = "")
  Loss_Data4$Label[Loss_Data4$delta_H == 0] <- "CH2"

  # Generate a unique identifier for looping over series
  Loss_Data4$Identifier <- paste(Loss_Data4$O, Loss_Data4$delta_H)

  # Normalize intensities and calculate the slope for each unique identifier
  Loss_Data5 <- lapply(unique(Loss_Data4$Identifier), function(x) {

    Data_Proces <- Loss_Data4[Loss_Data4$Identifier == x, ]
    Data_Proces <- Data_Proces[order(Data_Proces$Carbon_Position, decreasing = TRUE), ]
    Data_Proces$Slope <- NA
    Data_Proces$Slope_Scaled <- NA

    # Calculate the slope if there are multiple entries
    if (nrow(Data_Proces) > 1) {
      Data_Proces$Slope[2:nrow(Data_Proces)] <- diff(Data_Proces$Height)
      Data_Proces$Slope_Scaled <- rescale_max(Data_Proces$Slope)
    }

    # Scale intensities of a specified elemental formula loss series
    Data_Proces$scaled_intensity <- rescale_max(Data_Proces$Height)
    Data_Proces
  })

  # Combine the processed data into a final data frame
  Loss_Data6 <- do.call(rbind, Loss_Data5)

  # Return the data.frame
  Loss_Data6
}

#===============================================================================
# Peakwindow Adapt
#===============================================================================
#' Identify Maxima within Series with Low Data Points
#'
#' This function is adapted from the `peakwindow` function in the `cardiates`
#' package to identify peaks within series that contain low data points. It
#' allows for the detection of maximum points with flexibility in peak limits.
#'
#' @param x Numeric vector. The x-coordinates of the data points (e.g., time or
#'   m/z values).
#' @param y Numeric vector. The y-coordinates of the data points (e.g., intensity
#'   values).
#' @param xstart Numeric. The starting x-coordinate to consider for peak
#'   detection (default: 0).
#' @param xmax Numeric. The maximum x-coordinate for peak detection (default: max(x)).
#' @param minpeak Numeric. The minimum threshold for peak detection, specified as
#'   a fraction of the maximum y value (default: 0.1).
#' @param mincut Numeric. Threshold for cutting insignificant peaks; values below
#'   this fraction of the minimum value in each peak range are ignored (default: 0.382).
#'
#' @return A list of class `cardiPeakwindow`, containing the following components:
#' \item{peaks}{A data frame of detected peaks with indices, left and right bounds,
#'   x and y coordinates, and intensities.}
#' \item{data}{A data frame containing the original x and y data points.}
#' \item{smd.max.index}{The index of the peak with the maximum smoothed y value.}
#' \item{smd.max.x}{The x-coordinate of the peak with the maximum smoothed y value.}
#' \item{smd.indices}{Indices of the smoothed region.}
#' \item{smd.x}{x-coordinates of the smoothed region.}
#' \item{smd.y}{y-coordinates of the smoothed region.}
#' \item{peakid}{Peak identifiers for each point in the y data.}
#'
#' @examples
#' # Example usage of peakwindow_adapt function
#' peaks <- peakwindow_adapt(x = seq(1, 100), y = rnorm(100))
#'
#' @export
peakwindow_adapt <- function (x, y = NULL, xstart = 0, xmax = max(x), minpeak = 0.1,
                              mincut = 0.382) {

  # Helper function to count peaks within the data
  numpeaks <- function(ft, y) {
    npeaks <- length(ft)
    ndata <- length(y)
    if (npeaks > 0) {
      fp <- x1 <- x2 <- numeric(npeaks)
      peakid <- numeric(length(y))
      for (i in 1:npeaks) {
        fp[i] <- ft[i]
        mp <- match(fp[i], ft)
        if (mp == 1) x1[i] <- 1 else x1[i] <- ft[mp - 1] + which.min(y[ft[mp - 1]:fp[i]]) - 1
        if (mp == npeaks) x2[i] <- ndata else x2[i] <- fp[i] + which.min(y[fp[i]:ft[mp + 1]]) - 1
        peakid[x1[i]:x2[i]] <- i
      }
    } else {
      fp <- ft
      x1 <- x2 <- NULL
      peakid <- rep(0, ndata)
    }
    list(fp = fp, x1 = x1, x2 = x2, id = peakid)
  }

  xy <- xy.coords(x, y)
  x <- xy$x
  y <- xy$y
  iend <- max(c(1, which(x <= xmax)))
  tp <- turnpoints(y)
  ft <- tp$tppos

  # Modification: Ensure tp$firstispeak is not NA
  if (all(is.na(ft)) | is.na(tp$firstispeak)) {
    x1 <- 1
    x2 <- length(y)
    fp <- NULL
    ft <- NULL
  } else {
    ft <- ft[seq(2 - 1 * tp$firstispeak, length(ft), by = 2)]
    if (y[1] > y[2]) ft <- c(1, ft)
    if (y[iend] > y[iend - 1]) ft <- c(ft, iend)
    ft <- ft[ft <= iend]
    ftsmall <- which(y < minpeak * max(y))
    ft <- ft[!ft %in% ftsmall]
    dosearch <- ifelse(length(ft) > 1, TRUE, FALSE)
    while (dosearch) {
      fl <- length(ft)
      eli <- NULL
      for (n in 1:(fl - 1)) {
        km <- which.min(y[ft[n:(n + 1)]])
        if (min(y[ft[n:(n + 1)]]) * mincut < min(y[ft[n]:ft[n + 1]]))
          eli <- c(eli, n + km - 1)
      }
      if (length(eli) == 0) dosearch <- FALSE else ft <- ft[-eli]
      if (length(ft) < 2) dosearch <- FALSE
    }
    if (length(ft) > 1) {
      if (any(x[ft] > xstart)) {
        fp <- min(ft[x[ft] > xstart])
      } else {
        fp <- max(ft[x[ft] <= xstart])
      }
      mp <- match(fp, ft)
      if (mp == 1) x1 <- 1 else x1 <- ft[mp - 1] + which.min(y[ft[mp - 1]:fp]) - 1
      if (mp == length(ft)) x2 <- iend else x2 <- fp + which.min(y[fp:ft[mp + 1]]) - 1
    } else {
      fp <- ft
      x1 <- 1
      x2 <- iend
    }
  }

  allpeaks <- numpeaks(ft, y)
  ret <- list(peaks = data.frame(index = ft, xleft = allpeaks$x1,
                                 x = x[ft], xright = allpeaks$x2, y = y[ft]),
              data = data.frame(x = x, y = y),
              smd.max.index = fp, smd.max.x = x[fp],
              smd.indices = x1:x2, smd.x = x[x1:x2], smd.y = y[x1:x2],
              peakid = allpeaks$id)

  class(ret) <- c("list", "cardiPeakwindow")
  ret
}

#===============================================================================
# Find Fatty Acid Losses from DG Fragments
#===============================================================================
#' Identify Carbon Losses in MS/MS Spectra
#'
#' This function identifies intensity maxima of carbon losses with specified
#' oxygen content, allowing for the determination of the number of carbon atoms
#' in fatty acids from a given MS/MS spectrum.
#'
#' @param Data A data frame containing the MS/MS spectrum data, including columns
#'   for `Height` (intensity), `C` (number of carbon atoms), and `O` (number of
#'   oxygen atoms).
#' @param Threeshold Numeric. The minimum intensity threshold for identifying peaks
#'   (default: 10).
#' @param Oxygen Numeric vector. A vector of allowable oxygen content values to
#'   filter the data (default: c(1, 2)).
#' @param median_Factor Numeric. A factor used to normalize the maximum summed
#'   height against the median (default: 4).
#'
#' @return A vector of carbon fragment counts (`C`) that exceed the specified
#'   criteria. If no peaks meet the criteria, returns the message "No Peak".
#'
#' @examples
#' # Example usage of Find_FA_Losses function
#' carbon_losses <- Find_FA_Losses(Data = spectrum_data)
#'
#' @export
Find_FA_Losses <- function(Data, Threeshold = 10, Oxygen = c(1, 2), median_Factor = 4) {

  # Filter the Data to only include rows where Height is greater than the threshold
  # and the Oxygen content is in the specified range
  Data <- Data[Data$Height > Threeshold & Data$O %in% Oxygen, ]

  # Check if there are any rows left after filtering
  if (nrow(Data) > 0) {

    # Summarize the total height for each carbon fragment (grouped by Carbon count 'C')
    Sum_Sum <- Data %>%
      group_by(C) %>%
      summarise(sum = sum(Height, na.rm = TRUE))  # Calculate the sum of heights, ignoring NA values

    # Convert the summarized data into a data frame
    Sum_Sum <- as.data.frame(Sum_Sum)

    # Calculate the median of the summed heights
    med <- median(Sum_Sum$sum, na.rm = TRUE)

    # Find the maximum summed height
    max <- max(Sum_Sum$sum, na.rm = TRUE)

    # Identify carbon fragments that exceed a threshold based on the median and maximum values
    Carbons <- Sum_Sum$C[Sum_Sum$sum > max / med / median_Factor * med]

    # Return the carbon fragments that meet the criteria
    Carbons
  } else {
    # Return a message indicating no peaks were found if no data passed the filtering
    "No Peak"
  }
}



#===============================================================================
# Extract Fatty Acid Candidates from DG Fragments
#===============================================================================
#' Extract Fatty Acid (FA) Candidates from Peak Data
#'
#' This function extracts potential fatty acid candidates from a given dataset based
#' on specified peaks and conditions. It generates combinations of peaks, filters
#' them based on precursor data, and computes related properties. Only combinations
#' that match the precursor mass will be considered.
#'
#' @param Data A data frame containing the MS/MS spectrum data, including columns
#'   for `C` (number of carbon atoms) and `O` (number of oxygen atoms).
#' @param FA_Peaks A numeric vector of fatty acid peaks to be extracted.
#' @param Fa_Number Numeric. The number of fatty acid components to consider,
#'   calculated as `Precursor_DF$O - 3`.
#' @param Precursor_Df A data frame containing precursor information, including
#'   columns for carbon and oxygen counts.
#' @param Backbone_C Numeric. The carbon count of the backbone, which is subtracted
#'   from the precursor carbon count (default: 3).
#' @param O_Count Numeric vector. A vector of allowable oxygen counts to filter the
#'   data (default: c(1, 2)).
#'
#' @return A data frame containing potential fatty acid candidates based on the
#'   specified conditions. Returns NA if no fatty acid peaks are provided.
#'
#' @examples
#' # Example usage of Extract_FA_Losses function
#' potential_fa <- Extract_FA_Losses(Data = spectrum_data, FA_Peaks = c(12, 14, 16),
#'                                    Precursor_Df = precursor_data)
#'
#' @export
Extract_FA_Losses <- function(Data, FA_Peaks, Fa_Number = (Precursor_DF$O - 3), Precursor_Df,
                              Backbone_C = 3, O_Count = c(1, 2)) {

  # Check if there are any fatty acid peaks provided
  if (length(FA_Peaks) > 0) {
    # Repeat the fatty acid peaks according to the specified Fa_Number
    Vec <- rep(FA_Peaks, Fa_Number)

    # Generate all combinations of the repeated peaks of size Fa_Number
    Comb_Mat <- combn(Vec, Fa_Number)

    # Sort each combination and convert to a data frame
    Comb_Mat_Sort <- as.data.frame(apply(Comb_Mat, 2, sort))

    # Remove duplicated columns to ensure unique combinations
    Comb_Mat_Uni <- as.matrix(Comb_Mat_Sort[, !duplicated(t(Comb_Mat_Sort))])

    # Calculate the sum of each unique combination
    sum_vec <- apply(Comb_Mat_Uni, 2, sum)

    # Filter unique combinations where the sum matches the precursor carbon count minus the glycerol carbon count
    Pot_Fa_C <- unique(Comb_Mat_Uni[, sum_vec %in% (Precursor_Df$C - Backbone_C)])

    # Filter the original data for rows that match the potential fatty acid carbon counts and specified oxygen counts
    FA_df <- Data[Data$C %in% Pot_Fa_C & Data$O %in% O_Count, ]

    # Calculate the delta H for each row in the filtered data frame
    FA_df$delta_H <- ((FA_df$C - 1) * 2 + 1) - FA_df$H

    # Determine if delta_H is even and create a corresponding logical column
    FA_df$even <- FA_df$delta_H %% 2 == 0

    # Return the processed data frame with potential fatty acids
    return(FA_df)
  } else {
    # Return NA if no fatty acid peaks are provided
    return(NA)
  }
}


#===============================================================================
# Find Potential Fatty Acids from Monoglyceride Data
#===============================================================================
#' Identify Potential Fatty Acids (FAs) from Monoglyceride Data
#'
#' This function identifies potential fatty acids from monoglyceride (MG) data based
#' on specified thresholds and oxygen counts. It filters data based on height and
#' oxygen content, summarizes the height by carbon count, and determines which carbon
#' counts exceed a calculated threshold.
#'
#' @param Data A data frame containing the monoglyceride data, including columns for
#'   `Height`, `C` (number of carbon atoms), and `O` (number of oxygen atoms).
#' @param Threeshold Numeric. The minimum height value for a peak to be considered
#'   (default: 10).
#' @param Oxygen Numeric vector. The specific oxygen counts to filter the data on
#'   (default: 3).
#' @param median_Factor Numeric. A factor used to scale the median when determining
#'   threshold conditions (default: 4).
#'
#' @return A list where each element corresponds to the specified oxygen count and
#'   contains carbon counts that exceed the calculated threshold. Returns "No Peak"
#'   if no data meets the criteria for a given oxygen count.
#'
#' @examples
#' # Example usage of Find_FAs function
#' identified_fas <- Find_FAs(Data = mg_data, Oxygen = c(1, 2, 3))
#'
#' @export
Find_FAs <- function(Data, Threeshold = 10, Oxygen = 3, median_Factor = 4) {

  Output <- lapply(Oxygen, function(x){

    # Filter data to include only rows with Height above the threshold and the specified oxygen count
    Data <- Data[Data$Height > Threeshold & Data$O == x, ]

    # Check if there are any rows left after filtering
    if (nrow(Data) > 0) {
      # Summarize the total height by carbon count (C)
      Sum_Sum <- Data %>%
        group_by(C) %>%
        summarise(sum = sum(Height, na.rm = TRUE))

      # Convert the summarized data to a data frame
      Sum_Sum <- as.data.frame(Sum_Sum)

      # Calculate the median of the summed heights
      med <- median(Sum_Sum$sum, na.rm = TRUE)

      # Calculate the maximum of the summed heights
      max <- max(Sum_Sum$sum, na.rm = TRUE)

      # Identify carbon counts where the summed height exceeds a calculated threshold
      Carbons <- Sum_Sum$C[Sum_Sum$sum > max / med / median_Factor * med]

      # Return the identified carbon counts
      return(Carbons)
    } else {
      # Return a message if no data meets the criteria
      return("No Peak")
    }
  })

  names(Output) <- Oxygen
  as.vector(Output)
}

#===============================================================================
# Extract Fatty Acids from Monoacylglycerol and Fatty Acid Fragments
#===============================================================================
#' Extract Carbon Numbers of Fatty Acid Fragments
#'
#' This function extracts the carbon number of fatty acid fragments based on
#' monoacylglycerol (MG) and fatty acid fragments. It matches combinations of carbon
#' atoms to the precursor carbon count after glycerol carbon deduction.
#'
#' @param Data A data frame containing the data to be filtered, including columns for
#'   `C` (number of carbon atoms), `O` (number of oxygen atoms), and other relevant
#'   measurements.
#' @param FA_Peaks A numeric vector of fatty acid peaks identified in the dataset.
#' @param Fa_Number Numeric. The number of fatty acids to consider in combinations
#'   (default: calculated as `Precursor_DF$O - 3`).
#' @param Precursor_Df A data frame containing precursor information, including carbon counts.
#' @param Backbone_C Numeric. The number of carbons in the glycerol backbone (default: 3).
#' @param O_Count A numeric vector of oxygen counts to filter the data (default: 1:3).
#'
#' @return A data frame containing entries corresponding to identified potential fatty
#'   acids, including columns for `FA` (indicating potential fatty acid) and `delta_H`.
#'   Returns NA if no peaks are found.
#'
#' @examples
#' # Example usage of Extract_FAs function
#' extracted_fas <- Extract_FAs(Data = mg_data, FA_Peaks = fa_peaks)
#'
#' @export
Extract_FAs <- function(Data, FA_Peaks, Fa_Number = (Precursor_DF$O - 3), Precursor_Df,
                        Backbone_C = 3, O_Count = c(1:3)) {

  # Remove "No Peak" entries from FA_Peaks
  FA_Peaks <- FA_Peaks[!FA_Peaks %in% "No Peak"]

  if (length(FA_Peaks) > 0) {
    Glycerol_Adj <- O_Count[O_Count > 1]

    if (length(Glycerol_Adj) > 0) {
      Ox <- unlist(FA_Peaks[!names(FA_Peaks) %in% Glycerol_Adj])
      Ox2 <- unlist(FA_Peaks[names(FA_Peaks) %in% Glycerol_Adj]) - 3

      FA_Peaks <- if (length(Ox) > 0) c(Ox, Ox2) else Ox2
    }

    # Generate combinations of unique carbon counts for the specified number of fatty acids
    Vec <- rep(unique(FA_Peaks), Fa_Number)
    Comb_Mat <- combn(Vec, Fa_Number)
    Comb_Mat_Sort <- as.data.frame(apply(Comb_Mat, 2, sort))

    # Remove duplicated columns
    Comb_Mat_Uni <- as.matrix(Comb_Mat_Sort[, !duplicated(t(Comb_Mat_Sort))])

    # Calculate the sum of each combination
    sum_vec <- apply(Comb_Mat_Uni, 2, sum)

    # Identify potential carbon counts for fatty acids based on precursor data
    Pot_Fa_C <- unique(Comb_Mat_Uni[, sum_vec %in% (Precursor_Df$C - Backbone_C)])

    # Filter data to find entries corresponding to identified potential fatty acids
    FA_df <- Data[Data$C %in% c(Pot_Fa_C, Pot_Fa_C + 3) & Data$O %in% O_Count, ]

    # Mark the entries that correspond to potential fatty acids
    FA_df$FA <- FA_df$C %in% Pot_Fa_C

    # Calculate the change in hydrogen count (delta_H)
    FA_df$delta_H <- ((FA_df$C - 1) * 2 + 1) - FA_df$H

    # Check if delta_H is even
    FA_df$even <- FA_df$delta_H %% 2 == 0

    # Return the resulting data frame
    return(FA_df)
  } else {
    # Return NA if no peaks are found
    return(NA)
  }
}

#===============================================================================
# Find Double Bonds in Data
#===============================================================================
#' Identify and Plot Data Related to Double Bonds
#'
#' This function identifies and plots data related to double bonds in a given
#' dataset. It filters the data based on specified conditions and creates slope
#' and height plots.
#'
#' @param Data A data frame containing columns for carbon counts (C), height,
#'   slope, and delta_H values. The data should include measurements related to
#'   fatty acids.
#' @param delta_H A numeric vector of delta_H values to filter the data.
#'   Default is c(0, 2, 4, 6).
#' @param Threeshold A numeric value representing the minimum height threshold
#'   for filtering data. Default is 10.
#' @param Oxygen A numeric value representing the oxygen count to filter data.
#'   Default is 0.
#' @param mincut A numeric value used as a threshold for peak detection.
#'   Default is 0.382.
#' @param xstart A numeric value representing the starting point for peak
#'   detection. Default is 0.
#' @param minpeak A numeric value for the minimum height of detected peaks.
#'   Default is 0.1.
#' @param text_size A numeric value for text size in plots. Default is 10.
#'
#' @return A list containing:
#'   - `Peaks_Height_Based`: Data related to height peaks.
#'   - `Peaks_Slope_Based`: Data related to slope peaks.
#'   - `DB_Figure`: A combined plot of height and slope for each delta_H value.
#'
#' @examples
#' # Example usage of Find_DB function
#' result <- Find_DB(Data = your_data_frame)
#' # This will filter your_data_frame and create plots for identified double bonds.
#'
#' @export
Find_DB <- function(Data, delta_H = c(0, 2, 4, 6), Threeshold = 10,
                    Oxygen = 0, mincut = 0.382, xstart = 0,
                    minpeak = 0.1, text_size = 10) {

  # Filter data to only include entries with positive carbon counts
  Data <- Data[Data$C > 0, ]

  # Further filter data based on the specified delta_H values
  Data <- Data[Data$delta_H %in% delta_H, ]

  # Replace NA values in Slope column with 0
  Data$Slope[is.na(Data$Slope)] <- 0

  # Calculate the slope based on the height differences
  Data$Slope[2:nrow(Data)] <- diff(Data$Height)

  # Generate plots for each unique delta_H value
  Output_List <- lapply(unique(Data$delta_H), function(x) {

    # Filter output for the current delta_H value and other conditions
    Output <- Data[Data$delta_H == x & Data$Height > Threeshold & Data$O == Oxygen, ]

    # Create parsed titles for the plots based on the delta_H value
    if (x != 0) {
      Parsed_Height <- parse(text = paste(expression('M' * '-' * 'CH'[3] * '(CH'[2] * ')'[n] * ' ' * Delta), x, "H", expression(" " * "(" * "Intensity" * " " * "Plot)"), sep = "*"))
      Parsed_Slope <- parse(text = paste(expression('M' * '-' * 'CH'[3] * '(CH'[2] * ')'[n] * ' ' * Delta), x, "H", expression(" " * "(" * "Slope" * " " * "Plot)"), sep = "*"))
    } else {
      Parsed_Slope <- parse(text = paste(expression('M' * '-' * 'CH'[3] * '(CH'[2] * ')'[n]), expression(" " * "(" * "Slope" * " " * "Plot)"), sep = "*"))
      Parsed_Height <- parse(text = paste(expression('M' * '-' * 'CH'[3] * '(CH'[2] * ')'[n]), expression(" " * "(" * "Intensity" * " " * "Plot)"), sep = "*"))
    }

    # Check if there are any valid data entries left
    if (nrow(Output) > 0) {

      # Identify peaks in height and slope using the custom peak detection function
      peaks_Height <- peakwindow_adapt(Output$C, Output$Height, mincut = mincut, xstart = xstart)
      peaks_Slope <- peakwindow_adapt(Output$C, Output$Slope, mincut = mincut, xstart = xstart)

      # Create Slope Plot
      A <- ggplot(Output, aes(x = as.numeric(C), y = Slope)) +
        geom_hline(yintercept = 0, linetype = "dashed", colour = "darkgrey", size = 1) +
        geom_point(size = 3) +
        theme(legend.position = "none",
              text = element_text(size = text_size),
              panel.background = element_rect(fill = "white", colour = "black"),
              panel.grid.minor.y = element_blank(),
              panel.border = element_blank(),
              axis.ticks.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.major.y = element_blank(),
              strip.background = element_rect(fill = "grey90", colour = "black")) +
        xlab("Carbon Loss") +
        ggtitle(Parsed_Slope) +
        scale_x_continuous(limits = c(1, 20), breaks = seq(1, max(Output$C), by = 1))

      # Add peak rectangles if peaks were detected
      if (nrow(peaks_Slope$peaks) > 0) {
        A <- A + geom_rect(data = peaks_Slope$peaks, inherit.aes = FALSE,
                           aes(xmin = x - 0.5, xmax = x + 0.5, ymin = -Inf, ymax = Inf,
                               group = index, fill = as.character(index)),
                           alpha = 0.4, colour = "black") +
          scale_fill_brewer(palette = "Set1")
      }

      # Create Height Plot
      B <- ggplot(Output, aes(x = as.numeric(C), y = Height)) +
        geom_point(size = 3) +
        theme(legend.position = "none",
              text = element_text(size = text_size),
              panel.background = element_rect(fill = "white", colour = "black"),
              panel.grid.minor.y = element_blank(),
              panel.border = element_blank(),
              axis.ticks.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.major.y = element_blank(),
              strip.background = element_rect(fill = "grey90", colour = "black")) +
        xlab("Carbon Loss") +
        ggtitle(Parsed_Height) +
        scale_x_continuous(limits = c(1, 20), breaks = seq(1, max(Output$C), by = 1))

      # Add peak rectangles for height if detected
      if (nrow(peaks_Height$peaks) > 0) {
        B <- B + geom_rect(data = peaks_Height$peaks, inherit.aes = FALSE,
                           aes(xmin = x - 0.5, xmax = x + 0.5, ymin = -Inf, ymax = Inf,
                               group = index, fill = as.character(index)),
                           alpha = 0.4, colour = "black") +
          scale_fill_brewer(palette = "Set1")
      }

      # Return the height peaks, slope peaks, and combined plots
      Output <- list(peaks_Height, peaks_Slope, B | A)
      names(Output) <- c("Peaks_Height_Based", "Peaks_Slope_Based", "DB_Figure")
      Output
    } else {
      # Return a message if no peaks were found
      "No Peak"
    }
  })

  # Name the output list based on unique delta_H values
  names(Output_List) <- unique(Data$delta_H)

  # Return the output list sorted by delta_H values
  Output_List[order(as.numeric(names(Output_List)))]
}

#===============================================================================
# Extract n-th Element from List
#===============================================================================
#' Extract the n-th Element from Each Item in a List
#'
#' This function extracts the n-th element from each item in a given list.
#'
#' @param lst A list from which elements will be extracted. Each item in the list
#'   should be a vector or a list itself.
#' @param n An integer representing the index of the element to extract from each
#'   item in the list. Indices start from 1.
#'
#' @return A vector containing the n-th elements extracted from each item in the
#'   list. If n exceeds the length of any item, NA will be returned for that item.
#'
#' @examples
#' # Example usage of Extract_List function
#' my_list <- list(c(1, 2, 3), c(4, 5, 6), c(7, 8, 9))
#' result <- Extract_List(my_list, 2)
#' # result will be c(2, 5, 8)
#'
#' @export
Extract_List <- function(lst, n) {
  sapply(lst, `[`, n)
}

#===============================================================================
# Add Secondary Double Bond Series
#===============================================================================
#' Add Secondary Double Bond Series Based on Specified Losses
#'
#' This function adds secondary double bond series to the provided data based
#' on specified losses. It calculates double bond equivalents and filters the
#' data accordingly.
#'
#' @param Data A data frame containing the original data. Default is a variable
#'   named `Data`.
#' @param Losses A string specifying the losses (chemical formulas or mass values)
#'   to be processed. Default is an empty string, which calculates losses based
#'   on FA peaks.
#' @param DBE_min Minimum double bond equivalent to consider. Default is -20.
#' @param C0 A boolean indicating whether to keep carbon count of zero. Default
#'   is TRUE.
#' @param FA_Peak_vec A vector of fatty acid peaks. Default is `FA_Peaks`.
#'
#' @return A list containing processed data frames and any invalid inputs. The
#'   list may contain:
#'   - `Sec_Data_DB_filtered`: Filtered data for calculated losses.
#'   - `Sec_Data_Loss_H3`: Filtered data with carbon count equal to zero.
#'   - Messages for invalid inputs if applicable.
#'
#' @examples
#' # Example usage of Find_DB2 function
#' result <- Find_DB2(Data = my_data, Losses = "C18H36O2")
#'
#' @export
Find_DB2 <- function(Data = Data, Losses = "", DBE_min = -20, C0 = TRUE, FA_Peak_vec = FA_Peaks) {

  # If no losses are provided, calculate losses based on FA peaks.
  if (Losses == "") {
    H_Calc <- (FA_Peak_vec - 1) * 2 + 1
    Losses <- paste("C", FA_Peak_vec, "H", H_Calc, "O2", sep = "", collapse = ",")
  }

  # Split the Losses string into individual loss components
  Sec_string <- unlist(strsplit(Losses, ","))

  # Separate losses into strings and numeric values
  Sec_Frags_string <- Sec_string[is.na(as.numeric(Sec_string))]
  Sec_Frags <- Sec_string[!is.na(as.numeric(Sec_string))]

  # Calculate mass for elemental formulas parsed from strings
  if (length(Sec_Frags_string) > 0) {
    Sec_Frags_string[grepl("[aA-zZ]+", Sec_Frags_string)] <-
      apply(as.matrix(Sec_Frags_string[grepl("[aA-zZ]+", Sec_Frags_string)]), 1, function(x) { getMolecule(x)$exactmass })

    Sec_Frags <- append(Sec_Frags, Sec_Frags_string)
  }

  # Round the mass values
  Sec_Frags <- round(as.numeric(Sec_Frags), 4)

  # Calculate the starting double bond mass based on precursor mass
  DB_Start <- as.numeric(Precursor_DF$exactmass) - Sec_Frags

  # Filter out invalid starting values (negative mass)
  Invalid <- DB_Start > 0
  DB_Start <- DB_Start[Invalid]

  # Process each valid starting double bond mass
  Secondary_Frag_List <- lapply(DB_Start, function(Z) {

    Sec_Data <- Data
    Sec_Data$loss <- Z - Sec_Data$Mass.Charge  # Calculate the loss based on the mass difference

    # Annotate and unify the data based on the calculated losses
    Sec_Data_DB <- Annotate_Losses(Sec_Data, DBE_Filter = -20, Keep_C0 = TRUE)
    Sec_Data_DB <- Unify_Formulas(Sec_Data_DB, Parameter = "Loss_Formula")

    # Filter the data for losses less than the starting mass
    Sec_Data_DB_filtered <- Sec_Data_DB[Sec_Data_DB$Mass.Charge < Z, ]

    # Extract entries with carbon count equal to zero
    Sec_Data_Loss_H <- Sec_Data_DB[Sec_Data_DB$C == 0, ]

    # Filter duplicates based on the loss formula
    Sec_Data_Loss_H2 <- lapply(unique(Sec_Data_Loss_H$Loss_Formula), function(x) {
      Output <- Sec_Data_Loss_H[Sec_Data_Loss_H$Loss_Formula %in% x, ]
      Output[which.max(Output$Height), ]
    })

    # Combine the filtered results into a single data frame
    Sec_Data_Loss_H3 <- do.call(rbind, Sec_Data_Loss_H2)

    # Return the filtered data sets
    list(Sec_Data_DB_filtered, Sec_Data_Loss_H3)
  })

  # Combine results from all processed double bond series into a single list
  Secondary_Frag_List2 <- unlist(Secondary_Frag_List, recursive = FALSE)
  names(Secondary_Frag_List2) <- paste(c("Sec_Data_DB_filtered:", "Sec_Data_Loss_H3"), rep(Sec_string[Invalid], each = 2), "loss")

  # Check for any invalid input (negative mass difference) and append to the results
  if (!length(Sec_string) == length(DB_Start)) {
    List <- as.list(paste("Mass Difference is negative for:", Sec_string[!Invalid]))
    names(List) <- paste("Invalid input:", Sec_string[!Invalid])

    Secondary_Frag_List2 <- append(Secondary_Frag_List2, List)
  }

  # Return the combined results
  Secondary_Frag_List2
}


#===============================================================================
# Translate Name to Elemental Formula
#===============================================================================
#' Convert Fatty Acid Name to Elemental Formula
#'
#' This function takes a string representing fatty acid names (e.g., "18:1") and
#' converts it into its corresponding elemental formula.
#'
#' @param Name A character string representing the fatty acid names in the format
#'   "number of carbons:number of double bonds" (e.g., "18:1").
#'
#' @return A character vector containing the elemental formula in the format
#'   "C[number of carbons]H[number of hydrogens]O2".
#'
#' @examples
#' # Example usage of Name_to_Formula function
#' formula <- Name_to_Formula("18:1")
#'
#' @export
Name_to_Formula <- function(Name) {

  # Extract unique fatty acid notations (e.g., "18:1")
  Unique_FAs <- unique(unlist(str_extract_all(Name, "\\d{1,2}:\\d+\\b")))

  # Split the extracted FA strings into two columns (number of carbons and double bonds)
  ADJ <- as.data.frame(do.call(rbind, (str_split(Unique_FAs, ":"))))

  # Convert the columns to numeric
  ADJ$V1 <- as.numeric(ADJ$V1)  # Number of carbons
  ADJ$V2 <- as.numeric(ADJ$V2)  # Number of double bonds

  # Calculate the elemental formula
  FA_Formula <- paste("C", ADJ$V1, "H", (ADJ$V1 - 1) * 2 + 1 - 2 * ADJ$V2, "O2", sep = "")

  FA_Formula  # Return the computed elemental formula
}
