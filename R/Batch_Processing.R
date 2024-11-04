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

#' Process Multiple Mass Spectrometry Data Files
#'
#' This function processes multiple mass spectrometry data files based on user-defined parameters.
#'
#' @param Spectra_Location A character string specifying the location of the spectra files. Defaults to the current working directory.
#' @param Sample_Information A character string or logical indicating whether sample information is provided (e.g., an Excel file). Defaults to FALSE.
#' @param Fragment_Annotation A logical indicating whether to perform fragment annotation. Defaults to TRUE.
#' @param Kendricks_Plot A logical indicating whether to generate Kendrick plots. Defaults to TRUE.
#' @param MS_Spectra A logical indicating whether to generate mass spectrum plots. Defaults to TRUE.
#' @param Find_FA_Attachment A logical indicating whether to find and extract fatty acid peaks. Defaults to FALSE.
#' @param Find_DB_from_FA A logical indicating whether to find double bonds from fatty acid data. Defaults to FALSE.
#' @param Loss_Annotation A logical indicating whether to annotate elemental formula losses. Defaults to TRUE.
#' @param Find_DBs A logical indicating whether to find double bonds based on fragmentation data. Defaults to TRUE.
#' @param Find_Additional_DB A logical indicating whether to find additional double bonds based on losses. Defaults to TRUE.
#' @param Create_Summary A logical indicating whether to create a summary of the results. Defaults to TRUE.
#' @param File.Column A character string specifying the column name in the sample information that contains file names. Defaults to "File.Name".
#' @param Analyte.Column A character string specifying the column name for analyte codes. Defaults to "Letter_Code".
#' @param Precursor.Column A character string specifying the column name for precursor masses. Defaults to "Precursor".
#' @param Loss.Column A character string specifying the column name for additional losses. Defaults to "Add_Losses".
#' @param Composition.Column A character string specifying the column name for analyte compositions. Defaults to "Analyte".
#'
#' @return A list containing processed data and plots for each spectrum file. The list includes:
#'   - Input_Data: The loaded data for the current spectrum.
#'   - Precursor_Data: The precursor data.
#'   - Annotated_Data: The annotated fragment data.
#'   - Kendricks_Data: The Kendrick transformed data.
#'   - MG_Fragment_Data: Data related to fatty acid fragments.
#'   - Annotated_Loss_Data: Data related to annotated losses.
#'   - DG_Fragment_Data: Data for double bond fragments.
#'   - DB_Data: Data for double bonds.
#'   - DB_After_Loss_Data: Data for additional double bonds after loss.
#'
#' @examples
#' # Basic usage
#' results <- Batch_Processing(Spectra_Location = "path/to/spectra", Sample_Information = "path/to/sample_info.xlsx")
#'
#' @export
Batch_Processing <- function(Spectra_Location = getwd(),
                             Sample_Information = F,
                             Fragment_Annotation = T,
                             Kendricks_Plot = T,
                             MS_Spectra = T,
                             Find_FA_Attachment = F,
                             Find_DB_from_FA = F,
                             Loss_Annotation = T,
                             Find_DBs = T,
                             Find_Additional_DB = T,
                             Create_Summary = T,
                             File.Column = "File.Name",
                             Analyte.Column = "Letter_Code",
                             Precursor.Column = "Precursor",
                             Loss.Column = "Add_Losses",
                             Composition.Column = "Analyte") {
  # Initialize a list to store batch conditions/settings
  Batch_Conditions <- list(Batch_Conditions = c(
    Spectra_Location = Spectra_Location,
    Sample_Information = Sample_Information,
    Fragment_Annotation = Fragment_Annotation,
    Kendricks_Plot = Kendricks_Plot,
    MS_Spectra = MS_Spectra,
    Find_FA_Attachment = Find_FA_Attachment,
    Find_DB_from_FA = Find_DB_from_FA,
    Loss_Annotation = Loss_Annotation,
    Find_DB = Find_DB,
    Find_Additional_DB = Find_Additional_DB,
    Create_Summary = Create_Summary,
    File_Column = File.Column,
    Analyte_Column = Analyte.Column,
    Precursor_Column = Precursor.Column,
    Loss_Column = Loss.Column,
    Composition.Column = Composition.Column
  ))

  if(!Spectra_Location == getwd()){

    setwd(Spectra_Location)

  }


  # Load spectra files based on whether sample information is provided
  if (Sample_Information == F) {
    # If no sample information is provided, list all .txt files in the specified location
    Spectra_Files <- list.files(path = Spectra_Location, pattern = ".txt")
  } else {
    # If sample information is provided, read it from the specified Excel file
    Info_File <- as.data.frame(read_xlsx(Sample_Information))
    Info_File[is.na(Info_File)] <- ""  # Replace NA values with empty strings
    Spectra_Files <- Info_File[, File.Column]  # Extract file names from specified column



    # Adjust any loss annotations in the Info_File
    Names_To_Adjust <- which(grepl(":", Info_File[, Loss.Column]))  # Find rows with loss annotations
    for (i in Names_To_Adjust) {
      Info_File[i, Loss.Column] <- paste(Name_to_Formula(Info_File[i, Loss.Column]), collapse = ",")  # Convert names to formulas
    }
    rm(Names_To_Adjust)  # Remove temporary variable
  }

  #Add folder location to spectra files if it is not within the working directory


  #=============================================================================
  # Loop Start: Process each spectrum file
  #=============================================================================
  Output <- lapply(Spectra_Files, function(Spectra) {
    Spec_Pos <- which(Spectra_Files == Spectra)  # Find the index of the current spectrum
    Perc <- Spec_Pos / length(Spectra_Files) * 100  # Calculate progress percentage

    # Print the current processing status
    print(paste(Spectra, " (Spectra ", Spec_Pos, " of ", length(Spectra_Files), ")", sep = ""))

    # Load data for the current spectrum file
    Data <- Load_Data(Dat.Name = Spectra, Precursor.Mass = ifelse(exists("Info_File"), Info_File[Spec_Pos, Precursor.Column], F))

    #==========================================================================
    # Without Fragment Annotation, most functions will not work
    #==========================================================================
    if (Fragment_Annotation) {
      Elemental_Formula_df <- Annotate_Fragments(Data)  # Annotate fragments for the loaded data

      if (Kendricks_Plot) {
        # Generate Kendrick plot if enabled
        Kendricks_Data <- Kendricks_Transformation(Kendricks_Data = Elemental_Formula_df)
        Kendrick_gg_colored <- Plot_Kendricks_Plot(Kendricks_Data, toggle_alpha = 0, Rel.Thres = 1 / 10, Abs.Thres = Thres, plot_legend = F, point_size = 1)
      }  # Kendrick Plot

      if (Find_FA_Attachment) {
        # Find and extract fatty acid peaks if enabled
        MG_Peaks <- Find_FAs(Elemental_Formula_df, Oxygen = 1:3)
        MG_DF <- Extract_FAs(FA_Peaks = MG_Peaks, Precursor_Df = Precursor_DF, Data = Elemental_Formula_df)
        MG_Plot <- Plot_FAs(Data = MG_DF)


        if (Find_DB_from_FA) {
          # Find double bonds from fatty acid data if enabled
          FA_DB_Series <- Plot_DB_from_FA(Data = Elemental_Formula_df, column = "R_L", Max_C = max(MG_DF$C))
        }

      }  # MG Peaks

      if (MS_Spectra) {
        # Generate mass spectrum plot if enabled
        MS_Spec <- Plot_MS2(Kendricks_Data, Axis_Split_Size = 1, Split_Text_Size = 8, Label_Size = 8, Line_Size = 1)
        MS_SPEC_KENDRICK <- MS_Spec / Kendrick_gg_colored
      }  # MS Spec


    }  # Fragment Annotation related

    #==========================================================================
    # Without loss Annotation, most will not be executed
    #==========================================================================
    if (Loss_Annotation) {
      # Annotate elemental formula losses
      if (exists(Composition.Column)) {
        FA_Formula <- Name_to_Formula(Info_File[Spec_Pos, Composition.Column])  # Convert composition names to formulas
      }  # Composition Based adjustments

      Fragmentation_Data <- Annotate_Losses(Data = Data, maxForm = ifelse(exists(Composition.Column), max(FA_Formula), maxFormula_FA))

      if (Find_FA_Attachment) {
        # From elemental formula losses leading to DG fragments
        FA_Peaks <- Find_FA_Losses(Fragmentation_Data)
        FA_DF <- Extract_FA_Losses(Data = Fragmentation_Data, FA_Peaks = FA_Peaks, Precursor_Df = Precursor_DF)
        M_FA_Plot <- Plot_FA_Losses(Data = FA_DF)
      }  # Fatty Acid Attachment

      if (Find_DBs) {
        # Find double bonds based on fragmentation data
        Double_Bond_Data <- Find_DB(Fragmentation_Data, Threeshold = Thres, delta_H = 0:9)
        Extracted_Plots <- Extract_DB_Plots(Double_Bond_Data)
      }  # Find DB
    }  # Annotate Losses

    if (Find_Additional_DB & Find_FA_Attachment | !(Info_File[Spec_Pos, Loss.Column] == "")) {
      # Find additional double bonds based on losses or if fatty acid attachment is found
      Additional_DB_Series <- Find_DB2(Data = Data, Losses = if (exists("Info_File")) Info_File[Spec_Pos, Loss.Column] else "", FA_Peak_vec = FA_Loss_Peaks)

      Alt_DB_Series <- Extract_DB2(Data = Additional_DB_Series)
    }  # Additional DB Series

    # Compile data and plots into lists for output
    Data_List <- list(
      Input_Data = if (exists("Data")) Data else NULL,
      Precursor_Data = if (exists("Precursor_DF")) Precursor_DF else NULL,
      Annotated_Data = if (exists("Elemental_Formula_df")) Elemental_Formula_df else NULL,
      Kendricks_Data = if (exists("Kendricks_Data")) Kendricks_Data else NULL,
      MG_Fragment_Data = if (exists("MG_DF")) MG_DF else NULL,
      Annotated_Loss_Data = if (exists("Fragmentation_Data")) Fragmentation_Data else NULL,
      DG_Fragment_Data = if (exists("FA_DF")) FA_DF else NULL,
      DB_Data = if (exists("Double_Bond_Data")) Double_Bond_Data else NULL,
      DB_After_Loss_Data = if (exists("Additional_DB_Series")) Additional_DB_Series else NULL
    )

    Plot_List <- list(
      Spectra_Kendrick_Plot = if (exists("MS_SPEC_KENDRICK")) MS_SPEC_KENDRICK else NULL,
      MS_Spec_Plot = if (exists("MS_Spec")) MS_Spec else NULL,
      MG_Plot = if (exists("MG_Plot")) MG_Plot else NULL,
      Kendrick_Plot = if (exists("Kendrick_gg_colored")) Kendrick_gg_colored else NULL,
      DB_from_FA_Plot = if (exists("FA_DB_Series")) FA_DB_Series else NULL,
      DG_Plot = if (exists("M_FA_Plot")) M_FA_Plot else NULL,
      DB_Plot = if (exists("Extracted_Plots")) Extracted_Plots else NULL,
      DB_After_Loss_Plot = if (exists("Alt_DB_Series")) Alt_DB_Series else NULL
    )

    # Combine data and plots into an output list for the current spectrum
    Output_List <- list(Data_List = Data_List, Plot_List = Plot_List)
    Output_List
  })  # End of lapply

  # Name the output based on whether Info_File exists
  names(Output) <- if (exists("Info_File")) Info_File$File.Name else Spectra_Files

  # Combine batch conditions with output results
  Output <- append(Batch_Conditions, Output)
  if (exists("Info_File")) {
    Output <- append(list(Info_File = Info_File), Output)
  }

  # Return the final output containing all results and settings
  Output
}  # End of function
