#===============================================================================
# Create Summary
#===============================================================================
#' Create a Summary Plot for Double Bonds from Batch Results Data
#'
#' This function generates a summary plot based on the input batch results
#' data, focusing on double bonds. It allows for optional parameters to
#' specify expected double bonds and the desired plot order.
#'
#' @param Data A data frame containing the batch results with necessary columns
#'   and structures, including annotated loss data.
#' @param Expected_DB An optional parameter indicating the number of double
#'   bonds column within the info.file. If NA, it calculates expected double
#'   bonds based on the data (default is NA).
#' @param Plot.Order A character string specifying the column within the
#'   info.file that indicates the desired plot order (default is "Order").
#'
#' @return A ggplot object representing the summary plot of double bonds.
#'
#' @examples
#' # Example usage of Create_Summary function
#' summary_plot <- Create_Summary(Batch_Results, Expected_DB = "Max_DB_Col")
#'
#' @export
Create_Summary <- function(Data, Expected_DB = NA, Plot.Order = "Order") {

  # Remove columns that are not needed for analysis
  Data.adj <- Data[!names(Data) %in% "Info_File"]
  Data.adj <- Data.adj[!names(Data.adj) %in% "Batch_Conditions"]

  # Retrieve batch conditions and relevant file and analyte columns
  Batch_Cond <- Data[["Batch_Conditions"]]
  Info_File <- Data["Info_File"]
  File.Column <- Batch_Cond[["File_Column"]]
  Analyte.Column <- Batch_Cond[["Analyte_Column"]]

  # Extract annotated loss data from each dataset
  D_list <- lapply(Data.adj, function(x) {
    x$Data_List$Annotated_Loss_Data
  })

  # Assign names to the list
  names(D_list) <- names(Data.adj)

  # Combine the annotated loss data into a single data frame
  D_df <- do.call(rbind, D_list)

  # Clean up file names by removing file extensions
  D_df[, File.Column] <- sub("\\.[^.]*$", "", rownames(D_df))

  # Merge with the info file based on file and analyte columns
  if (!Analyte.Column == F) {
    Merge_df <- Data$Info_File[c("File.Name", Analyte.Column)]
    Merged_DF <- merge(D_df, Merge_df, by = File.Column)
    colnames(Merged_DF)[colnames(Merged_DF) == Analyte.Column] <- "Plot.Name"
  } else {
    colnames(Merged_DF)[colnames(Merged_DF) == File.Column] <- "Plot.Name"
  }

  # Determine expected double bonds based on double bond equivalence, if not provided
  if (is.na(Expected_DB)) {
    Expected_DB <- lapply(Data.adj, function(x) {
      Output <- x$Data_List$Precursor_Data
      Output$DBE - (Output$O - 3)
    })

    Merged_DF2 <- lapply(names(Expected_DB), function(x) {
      sub <- Merged_DF[Merged_DF$File.Name %in% x,]
      DB_exp <- unlist(Expected_DB[names(Expected_DB) == x])

      # Filter based on delta_H and expected double bonds
      sub2 <- sub[sub$delta_H > -1 & sub$delta_H < (DB_exp * 2),]
      sub2[lapply(sub2$delta_H, "%%", 2) == 0,]
    })

    Merged_DF3 <- do.call(rbind, Merged_DF2)

  } else {

    if (is.numeric(Expected_DB)) {
      Merged_DF2 <- Merged_DF[Merged_DF$delta_H > -1 & Merged_DF$delta_H < Expected_DB * 2,]
      Merged_DF3 <- Merged_DF2[lapply(Merged_DF2$delta_H, "%%", 2) == 0,]

    } else if (is.character(Expected_DB)) {
      Merged_DF2 <- lapply(Info_File$Info_File$File.Name, function(x) {
        DB_Infos <- Info_File$Info_File[Info_File$Info_File$File.Name == x,]
        sub <- Merged_DF[Merged_DF$File.Name %in% x,]

        sub2 <- sub[sub$delta_H > -1 & sub$delta_H < (as.numeric(DB_Infos[, Expected_DB]) * 2),]
        sub2[lapply(sub2$delta_H, "%%", 2) == 0,]
      })

      Merged_DF3 <- do.call(rbind, Merged_DF2)
    }
  }

  if (is.character(Plot.Order)) {
    Merged_DF3$Plot.Name <- factor(Merged_DF3$Plot.Name, levels = unique(Merged_DF3$Plot.Name)[order(unique(Info_File$Info_File[, Plot.Order]), decreasing = F)])
  }

  Merged_DF3$delta_H_parsed <- ifelse(
    Merged_DF3$delta_H == 0,
    "M-CH3(CH2)n",
    paste("M-CH3(CH2)n ", Merged_DF3$delta_H, "H", sep = "")
  )

  # Create the summary plot using ggplot2
  ggSummary <- ggplot(Merged_DF3, aes(y = Plot.Name, x = Carbon_Position, fill = scaled_intensity)) +
    geom_tile(color = "black") +
    facet_wrap(~delta_H_parsed, ncol = 1, scales = "free_y")

  ggSummary +
    theme(text = element_text(colour = "black"),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          legend.position = "none",
          panel.background = element_rect(fill = "darkgrey", colour = "black"),
          panel.grid.minor.y = element_blank(),
          panel.border = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank(),
          strip.background = element_rect(fill = "grey90", colour = "black"),
          plot.title = element_blank(),
          axis.title.x = element_blank()) +
    ylab("") +
    xlab("")
}
