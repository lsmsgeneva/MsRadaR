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


#' Generate a Kendrick Plot
#'
#' This function creates a Kendrick plot using ggplot2, visualizing the relationship
#' between Kendrick Mass (Kendrick Defect) based on CH2. The plot differentiates data points
#' by color according to their oxygen content and adjusts point transparency based on
#' relative intensity thresholds.
#'
#' @param data A data frame containing at least the following columns:
#'   - `Height`: Signal intensity for the data points.
#'   - `Rel.Int`: Relative intensity of the data points.
#'   - `O`: Oxygen content.
#'   - `Ke_CH2`: Kendrick Mass values.
#'   - `KDM_CH2`: Kendrick Defect values.
#' @param toggle_alpha A numeric value indicating whether to adjust alpha levels for points
#'   below the relative intensity threshold. Defaults to 0.05.
#' @param Rel.Thres A numeric value representing the relative intensity threshold for
#'   adjusting alpha. Defaults to 0.01.
#' @param Abs.Thres A numeric value representing the absolute signal threshold for data points.
#' @param plot_legend A logical value indicating whether to display the legend on the plot.
#'   Defaults to TRUE.
#' @param point_size A numeric value for the size of the points in the plot. Defaults to 1.
#'
#' @return A ggplot object representing the Kendrick plot. The plot visualizes Kendrick Mass
#'   versus Kendrick Defect, with points colored by oxygen content and alpha transparency
#'   based on relative intensity.
#'
#' @examples
#' # Assuming 'kendrick_data' is a data frame containing the required columns
#' plot <- Plot_Kendricks_Plot(data = kendrick_data)
#' print(plot)
#'
#' @export
Plot_Kendricks_Plot <- function(data, toggle_alpha = 0.05, Rel.Thres = 1/100,
                                Abs.Thres = Signal_Noise, plot_legend = TRUE,
                                point_size = 1) {

  # Filter the input data to include only rows with Height greater than the absolute threshold
  data <- data[data$Height > Abs.Thres, ]

  # Create a new factor column indicating if the relative intensity exceeds the relative threshold
  data$alpha_thres <- as.factor(data$Rel.Int > Rel.Thres)

  # Set factor levels for the alpha threshold variable
  data$alpha_thres <- factor(data$Rel.Int > Rel.Thres, levels = c(FALSE, TRUE))

  # Determine the minimum and maximum values for the x-axis limits (Kendrick Mass)
  xmin <- min(data$Mass.Charge)
  xmax <- max(data$Mass.Charge)

  # Check if the toggle_alpha is greater than 0 for adjusted alpha levels in the plot
  if (toggle_alpha > 0) {

    # Create the ggplot object with points colored by oxygen content and transparency based on relative intensity
    Output <- ggplot(data, aes(x = Ke_CH2, y = KDM_CH2, colour = as.character(O), alpha = alpha_thres)) +
      geom_point(size = point_size * 5 / 14) +  # Adjust point size based on the provided parameter
      theme(
        panel.background = element_rect(fill = "white", colour = "black"),  # Set the background and border color
        panel.grid.minor.y = element_blank(),  # Remove minor grid lines on the y-axis
        panel.border = element_blank(),  # Remove panel border
        axis.ticks.x = element_blank(),  # Remove x-axis ticks
        panel.grid.minor.x = element_blank(),  # Remove minor grid lines on the x-axis
        panel.grid.major.x = element_blank(),  # Remove major grid lines on the x-axis
        panel.grid.major.y = element_blank(),  # Remove major grid lines on the y-axis
        strip.background = element_rect(fill = "grey90", colour = "black")  # Set strip background color
      ) +
      scale_colour_manual(values = c("6" = "darkred", "5" = "darkblue", "4" = "red2",
                                     "3" = "seagreen", "2" = "yellow3", "1" = "slateblue2",
                                     "0" = "darkorange2")) +  # Define colors for each oxygen level
      scale_alpha_manual(values = c("FALSE" = toggle_alpha, "TRUE" = 1)) +  # Set alpha levels for points
      xlab("Kendricks Mass") +  # X-axis label
      ylab("KMD based on CH2") +  # Y-axis label
      ggtitle(paste("Absolute Threshhold: ", Abs.Thres, ", Relative Threshold: ", Rel.Thres, "%", sep = "")) +  # Title with thresholds
      xlim(xmin, xmax) +  # Set x-axis limits
      theme(plot.title = element_text(size = 10))  # Adjust title font size

  } else {
    # If toggle_alpha is not greater than 0, filter data to keep only rows with relative intensity exceeding the threshold
    data <- data[data$Rel.Int > Rel.Thres, ]

    # Create the ggplot object without alpha transparency
    Output <- ggplot(data, aes(x = Ke_CH2, y = KDM_CH2, colour = as.character(O))) +
      geom_point(size = point_size * 5 / 14) +  # Adjust point size
      theme(
        panel.background = element_rect(fill = "white", colour = "black"),  # Set background and border color
        panel.grid.minor.y = element_blank(),  # Remove minor grid lines on the y-axis
        panel.border = element_blank(),  # Remove panel border
        axis.ticks.x = element_blank(),  # Remove x-axis ticks
        panel.grid.minor.x = element_blank(),  # Remove minor grid lines on the x-axis
        panel.grid.major.x = element_blank(),  # Remove major grid lines on the x-axis
        panel.grid.major.y = element_blank(),  # Remove major grid lines on the y-axis
        strip.background = element_rect(fill = "grey90", colour = "black")  # Set strip background color
      ) +
      scale_colour_manual(values = c("6" = "darkred", "5" = "darkblue", "4" = "red2",
                                     "3" = "seagreen", "2" = "yellow3", "1" = "slateblue2",
                                     "0" = "darkorange2")) +  # Define colors for each oxygen level
      xlab("Kendricks Mass") +  # X-axis label
      ylab("KMD based on CH2") +  # Y-axis label
      ggtitle(paste("Absolute Threshhold: ", Abs.Thres, ", Relative Threshold: ", Rel.Thres, "%", sep = "")) +  # Title with thresholds
      xlim(xmin, xmax) +  # Set x-axis limits
      theme(plot.title = element_text(size = 10))  # Adjust title font size
  }

  # Conditionally add a legend based on the plot_legend parameter
  if (plot_legend) {
    Output  # Return the plot with the legend
  } else {
    Output + theme(legend.position = "none")  # Return the plot without the legend
  }
}  # End of function



#' Generate a Bar Plot for Fatty Acid Losses
#'
#' This function creates a bar plot to visualize potential fatty acid attachments (FAs)
#' using delta hydrogen counts (`delta_H`) and height values, categorized by carbon
#' and oxygen counts. The plot distinguishes between different fragmentation mechanisms
#' (radical/neutral losses) using even/odd hydrogen counts.
#'
#' @param Data A data frame containing at least the following columns:
#'   - `H`: Hydrogen count.
#'   - `delta_H`: Delta hydrogen count representing changes due to losses.
#'   - `Height`: Height values for the data points.
#'   - `O`: Oxygen count.
#'   - `C`: Carbon count.
#' @param xlim1 A numeric value representing the lower limit for `delta_H`. Defaults to -5.
#' @param xlim2 A numeric value representing the upper limit for `delta_H`. Defaults to 10.
#' @param Oxygen_Range A numeric vector specifying the range of oxygen counts to include.
#'   Defaults to 0:10 (i.e., includes oxygen counts from 0 to 10).
#'
#' @return A ggplot object representing the bar plot of fatty acid losses.
#'   If no data matches the criteria, a message indicating this is returned.
#'
#' @examples
#' # Assuming 'fa_data' is a data frame containing the required columns
#' plot <- Plot_FA_Losses(Data = fa_data)
#' print(plot)
#'
#' @export
Plot_FA_Losses <- function(Data, xlim1 = -5, xlim2 = 10, Oxygen_Range = 0:10) {

  # Determine if the hydrogen count (H) is even and create a new column for it
  if(nrow(Data) > 0) {
    Data$Loss.Type <- NA

    Data$Loss.Type[Data$H %% 2 == 0] <- "Radical loss"
    Data$Loss.Type[!Data$H %% 2 == 0] <- "Neutral loss"

    # Create descriptive labels for oxygen and carbon counts
    Data$O_Paste <- paste(Data$O, "Oxygen lost")
    Data$C_Paste <- paste(Data$C, "C lost")

    # Filter the data based on delta_H limits and oxygen count range
    Data <- Data[Data$delta_H > xlim1 & Data$delta_H < xlim2, ]
    Data <- Data[Data$O %in% Oxygen_Range, ]

    # Group data by carbon count and scale height values
    Datax <- Data %>% group_by(C) %>% mutate(Scaled = rescale_max(Height))

    # Create the bar plot using ggplot2
    ggOut <- ggplot(Datax, aes(y = Height, x = delta_H, fill = Loss.Type, label = loss)) +
      geom_hline(yintercept = 0) +
      geom_vline(xintercept = 0, linetype = "dashed", size = 1, colour = "grey70") +
      geom_col(position = position_dodge2(preserve = "single"), color = "black") +
      facet_grid(C_Paste ~ O_Paste, scales = "free_y") +  # Create facets based on carbon and oxygen
      theme(
        panel.background = element_rect(fill = "white", colour = "black"),  # Set background and border color
        panel.grid.minor.y = element_blank(),  # Remove minor grid lines on the y-axis
        panel.border = element_blank(),  # Remove panel border
        axis.ticks.x = element_blank(),  # Remove x-axis ticks
        panel.grid.minor.x = element_blank(),  # Remove minor grid lines on the x-axis
        panel.grid.major.x = element_blank(),  # Remove major grid lines on the x-axis
        panel.grid.major.y = element_blank(),  # Remove major grid lines on the y-axis
        strip.background = element_rect(fill = "grey90", colour = "black"),  # Set strip background color
        plot.title = element_blank()  # No title
      ) +
      scale_fill_manual(values = c("steelblue", "darkgrey")) +  # Set custom colors for the fill
      xlab("Hydrogen difference from a saturated fatty acid loss")

    # Add text labels to the bars, adjusting position based on scaled height
    ggOut + geom_text(angle = 90, hjust = ifelse(Datax$Scaled > 0.5, 1.1, -0.1))
  } else {
    return("No matching fatty acids found")
  }
}


#' Generate a Bar Plot for Monoacylglycerol and Fatty Acid Fragments
#'
#' This function creates a bar plot to visualize monoacylglycerol and fatty acid
#' fragments based on their delta hydrogen counts (`delta_H`) and height values,
#' categorized by carbon and oxygen counts. It distinguishes between different
#' fragmentation mechanisms (radical/neutral losses) using even/odd hydrogen counts.
#'
#' @param Data A data frame containing at least the following columns:
#'   - `H`: Hydrogen count.
#'   - `delta_H`: Delta hydrogen count representing changes due to losses.
#'   - `Height`: Height values for the data points.
#'   - `O`: Oxygen count.
#'   - `C`: Carbon count.
#'   - `Mass.Charge`: Mass charge information to be displayed on the plot.
#' @param xlim1 A numeric value representing the lower limit for `delta_H`. Defaults to -10.
#' @param xlim2 A numeric value representing the upper limit for `delta_H`. Defaults to 5.
#' @param Oxygen_Range A numeric vector specifying the range of oxygen counts to include.
#'   Defaults to 0:10 (i.e., includes oxygen counts from 0 to 10).
#'
#' @return A ggplot object representing the bar plot of monoacylglycerol and fatty acid fragments.
#'   If no data matches the criteria, a message indicating this is returned.
#'
#' @examples
#' # Assuming 'fa_data' is a data frame containing the required columns
#' plot <- Plot_FAs(Data = fa_data)
#' print(plot)
#'
#' @export
Plot_FAs <- function(Data, xlim1 = -10, xlim2 = 5, Oxygen_Range = 0:10) {

  if(nrow(Data) > 0) {
    # Determine if the hydrogen count (H) is even and create a new column for it
    Data$Loss.Type <- NA
    Data$Loss.Type[Data$H %% 2 == 0] <- "Radical loss"
    Data$Loss.Type[!Data$H %% 2 == 0] <- "Neutral loss"

    # Create descriptive labels for oxygen and carbon counts
    Data$O_Paste <- paste(Data$O, "Oxygen")
    Data$C_Paste <- paste(Data$C, "C")

    # Invert delta_H values for plotting
    Data$delta_H <- Data$delta_H * -1

    # Filter the data based on delta_H limits and oxygen count range
    Data <- Data[Data$delta_H > xlim1 & Data$delta_H < xlim2, ]
    Data <- Data[Data$O %in% Oxygen_Range, ]

    # Group data by carbon count and scale height values
    Datax <- Data %>% group_by(C) %>% mutate(Scaled = rescale_max(Height))

    # Create the bar plot using ggplot2
    ggOut <- ggplot(Datax, aes(y = Height, x = delta_H, fill = Loss.Type, label = Mass.Charge)) +
      geom_hline(yintercept = 0) +
      geom_vline(xintercept = 0, linetype = "dashed", size = 1, colour = "grey70") +
      geom_col(position = position_dodge2(preserve = "single"), color = "black") +
      facet_grid(C_Paste ~ O_Paste, scales = "free_y") +  # Create facets based on carbon and oxygen
      theme(
        panel.background = element_rect(fill = "white", colour = "black"),  # Set background and border color
        panel.grid.minor.y = element_blank(),  # Remove minor grid lines on the y-axis
        panel.border = element_blank(),  # Remove panel border
        axis.ticks.x = element_blank(),  # Remove x-axis ticks
        panel.grid.minor.x = element_blank(),  # Remove minor grid lines on the x-axis
        panel.grid.major.x = element_blank(),  # Remove major grid lines on the x-axis
        panel.grid.major.y = element_blank(),  # Remove major grid lines on the y-axis
        strip.background = element_rect(fill = "grey90", colour = "black"),  # Set strip background color
        plot.title = element_blank()  # No title
      ) +
      scale_fill_manual(values = c("steelblue", "darkgrey")) +  # Set custom colors for the fill
      xlab("Hydrogen difference from a saturated fatty acid")

    # Add text labels to the bars, adjusting position based on scaled height
    ggOut + geom_text(angle = 90, hjust = ifelse(Datax$Scaled > 0.5, 1.1, -0.1))
  } else {
    return("No matching fatty acids found")
  }
}





#' Generate an Annotated MS/MS Plot
#'
#' This function creates an annotated MS/MS plot, visualizing relative intensities
#' of fragments based on mass-to-charge ratio (m/z) and oxygen count. The plot
#' automatically adjusts the y-axis to zoom in on lower intensity fragments when
#' there is a significant difference in intensities.
#'
#' @param Data A data frame containing at least the following columns:
#'   - `Rel.Int`: Relative intensity values of the fragments.
#'   - `Height`: Height values corresponding to the fragments (for reference).
#'   - `Mass.Charge`: Mass-to-charge ratio (m/z) for each fragment.
#'   - `O`: Oxygen count for each fragment.
#' @param Graph_Splitting_Ratio A numeric value that defines the ratio used to
#'   determine whether to split the graph for lower intensity fragments. Defaults to 1.1.
#' @param Axis_Split_Size A numeric value controlling the size of the axis split line.
#'   Defaults to 5.
#' @param Split_Text_Size A numeric value that sets the size of the split label text.
#'   Defaults to 1.
#' @param Label_Size A numeric value that specifies the size of the labels on the plot.
#'   Defaults to 8.
#' @param Line_Size A numeric value that determines the size of the line segments in the plot.
#'   Defaults to 1.
#' @param Rel_Threeshold A numeric threshold for relative intensity filtering.
#'   Defaults to 1.
#'
#' @return A ggplot object representing the annotated MS/MS plot. If significant intensity
#'   differences are not present, an unzoomed plot will be returned.
#'
#' @examples
#' # Assuming 'ms_data' is a data frame containing the required columns
#' plot <- Plot_MS2(Data = ms_data)
#' print(plot)
#'
#' @export
Plot_MS2 <- function(Data, Graph_Splitting_Ratio = 1.1, Axis_Split_Size = 5,
                     Split_Text_Size = 1, Label_Size = 8, Line_Size = 1,
                     Rel_Threeshold = 1) {

  # Sort relative intensity values and select the two highest
  Sorted_Values <- sort(Data$Rel.Int, decreasing = TRUE)[2:1]

  # Calculate the difference between the highest and second highest relative intensities
  Diff <- (diff(Sorted_Values) / Sorted_Values[1])

  # Determine the maximum height for the y-axis
  Maximum_Height <- round(max(Data$Height))

  # Check if the difference exceeds the defined graph splitting ratio
  if (Diff > Graph_Splitting_Ratio) {

    # Prepare the label for the difference
    Diffx <- paste("x", round(Diff, 1), sep = "")

    # Calculate the lower bound for visibility based on the highest relative intensity
    low_max <- ceiling(Sorted_Values[1])

    # Filter out values that won't be visible in the plot
    Data <- Data[Data$Relative_Intensity > low_max * Rel_Threeshold / 100, ]

    high_min <- 99  # Set a fixed upper limit for the relative intensity
    adjust <- high_min - low_max  # Calculate adjustment for the relative intensities

    # Adjust the relative intensity values for plotting
    Data <- Data %>%
      mutate(Rel.Int2 = as.numeric(Rel.Int),
             Rel.Int2 = case_when(Rel.Int < low_max ~ Rel.Int2,
                                  Rel.Int > high_min ~ Rel.Int2 - adjust,
                                  TRUE ~ NA_real_))

    # Get the mass-to-charge ratio (m/z) of the peak with the highest relative intensity
    M_Z_VAL <- Data$Mass.Charge[which.max(Data$Rel.Int)]

    # Create a ggplot with split ratios if significant differences are present
    ggplot(Data, aes(x = Mass.Charge, ymax = Rel.Int2, ymin = 0,
                     colour = as.factor(O), label = Mass.Charge, y = Rel.Int2)) +
      geom_linerange(size = Line_Size * 0.4699, key_glyph = draw_key_rect) +
      xlab("m/z") +
      annotate("segment", color = "grey80", size = Axis_Split_Size,
               x = -Inf, xend = Inf, y = low_max, yend = low_max) +
      geom_text_repel(box.padding = 0.1, max.overlaps = 5, point.size = NA,
                      nudge_y = 0.05, min.segment.length = 0, segment.size = 0.2,
                      angle = 90, force_pull = 0, hjust = 0, seed = 1,
                      show.legend = FALSE, size = Label_Size * 5 / 14, verbose = TRUE) +
      annotate(geom = "label", x = M_Z_VAL, y = low_max, label = Diffx,
               size = Label_Size * 5 / 14) +
      theme(axis.title.x = element_text(face = "italic"),
            legend.position = "top",
            panel.background = element_rect(fill = "white", colour = "black"),
            panel.grid.minor.y = element_blank(),
            panel.border = element_blank(),
            axis.ticks.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.major.y = element_blank(),
            strip.background = element_rect(fill = "grey90", colour = "black")) +
      scale_colour_manual(values = c("6" = "darkred", "5" = "darkblue",
                                     "4" = "red2", "3" = "seagreen",
                                     "2" = "yellow3", "1" = "slateblue2",
                                     "0" = "darkorange2")) +
      scale_y_continuous(expand = expansion(mult = c(0, 0)),
                         breaks = c(seq(0, 100, 2), 100.5),
                         label = function(x) {x + ifelse(x >= low_max, adjust, 0)}) +
      guides(color = guide_legend(nrow = 1, title = "Oxygen")) +
      ylab("Relative intensity") +
      ggtitle(paste("% of", Maximum_Height)) +
      theme(plot.title = element_text(hjust = 0, vjust = -17, size = 10))

  } else {
    # Create an unzoomed ggplot, in case no significant intensity differences are present
    ggplot(Data, aes(x = Mass.Charge, ymax = Rel.Int, ymin = 0,
                     colour = as.factor(O), label = Mass.Charge, y = Rel.Int)) +
      geom_text_repel(box.padding = 0.1, max.overlaps = 5, point.size = NA,
                      nudge_y = 0.05, min.segment.length = 0, segment.size = 0.2,
                      angle = 90, force_pull = 0, hjust = 0, seed = 1,
                      show.legend = FALSE, size = Label_Size * 5 / 14, verbose = TRUE) +
      geom_linerange(size = Line_Size * 0.469, key_glyph = draw_key_rect) +
      xlab("m/z") +
      theme(axis.title.x = element_text(face = "italic"),
            legend.position = "top",
            panel.background = element_rect(fill = "white", colour = "black"),
            panel.grid.minor.y = element_blank(),
            panel.border = element_blank(),
            axis.ticks.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.major.y = element_blank(),
            strip.background = element_rect(fill = "grey90", colour = "black")) +
      scale_colour_manual(values = c("6" = "darkred", "5" = "darkblue",
                                     "4" = "red2", "3" = "seagreen",
                                     "2" = "yellow3", "1" = "slateblue2",
                                     "0" = "darkorange2")) +
      scale_y_continuous(expand = expansion(mult = c(0, 0))) +
      guides(color = guide_legend(nrow = 1, title = "Oxygen")) +
      ylab("Relative intensity") +
      ggtitle(paste("% of", Maximum_Height)) +
      theme(plot.title = element_text(hjust = 0, vjust = -17, size = 10))
  }
}



#' Extract Double Bond Plot Objects
#'
#' This function extracts double bond plot objects created by the `Find_Double_Bonds`
#' functions. It filters and returns a series of plots based on specified criteria,
#' wrapped into a single output.
#'
#' @param Data_List A list of plot objects, potentially containing entries labeled
#'   as "No Peak". This list is generated by the `Find_Double_Bonds` function.
#' @param Figure_Pos An integer indicating the position of the plot to extract from the
#'   list. Defaults to 3.
#' @param delta_H A specification for filtering plots based on the number of double bonds.
#'   It can be:
#'   - A numeric value to extract specific double bond counts.
#'   - A string "even" to filter for even counts.
#'   - A string "odd" to filter for odd counts. Defaults to "even".
#'
#' @return A wrapped plot object containing the extracted double bond plots. If no valid
#'   data entries remain after filtering, a message indicating "No Data" will be returned.
#'
#' @examples
#' # Assuming 'data_list' is a list of plot objects obtained from Find_Double_Bonds
#' extracted_plots <- Extract_DB_Plots(Data_List = data_list, Figure_Pos = 3, delta_H = "even")
#' print(extracted_plots)
#'
#' @export
Extract_DB_Plots <- function(Data_List, Figure_Pos = 3, delta_H = "even") {

  # Filter out any entries labeled as "No Peak"
  Data_List <- Data_List[!Data_List %in% "No Peak"]

  # Check if there are any valid data entries left after filtering
  if (length(Data_List) > 0) {
    # Extract the relevant data based on the specified figure position
    DB_List <- Extract_List(Data_List, Figure_Pos)

    # Filtering Rule for extraction of series
    if (is.numeric(delta_H)) {
      Nos <- as.numeric(str_extract(names(DB_List), "[0-9]+"))
      DB_List <- DB_List[Nos %in% delta_H]
    } else {
      if (delta_H == "even") {
        Nos <- as.numeric(str_extract(names(DB_List), "[0-9]+"))
        DB_List <- DB_List[Nos %% 2 == 0]
      }
      if (delta_H == "odd") {
        Nos <- as.numeric(str_extract(names(DB_List), "[0-9]+"))
        DB_List <- DB_List[!Nos %% 2 == 0]
      }
    }

    # Order the extracted data list by names
    DB_List2 <- DB_List[order(names(DB_List))]

    # Wrap the plots into a single output with one column
    wrap_plots(DB_List2, ncol = 1)
  } else {
    # Return a message indicating that there is no data to plot
    "No Data"
  }
}


#' Generate a Heatmap of Double Bond Data
#'
#' This function creates a heatmap visualizing double bond data based on specified
#' parameters. It can fill tiles based on either scaled intensity or the logarithm
#' of height values.
#'
#' @param data A data frame containing the information to be plotted. It must include
#'   the columns `delta_H`, `C` (Carbon_Position), `Height`, and `O` (oxygen levels).
#' @param Oxygen_Range A vector specifying the range of oxygen values to include in
#'   the plot. Default is 0 to 6.
#' @param filling A string determining how the tiles in the heatmap are filled. Options
#'   are "scaled" for scaled intensity or "log10" of Height. Defaults to "scaled".
#' @param H_Count_min Minimum value for `delta_H` to include in the plot. Default is -10.
#' @param H_Count_Max Maximum value for `delta_H` to include in the plot. Default is 10.
#' @param Carbon_Pos_min Minimum value for `C` (Carbon_Position) to include in the plot.
#'   Default is 2.
#' @param Carbon_Pos_Max Maximum value for `C` (Carbon_Position) to include in the plot.
#'   Default is 20.
#'
#' @return A heatmap generated by ggplot2, showing the relationship between `delta_H`,
#'   `Carbon_Position`, and intensity values based on the specified filling method.
#'
#' @examples
#' # Assuming 'data_frame' is a data frame with the necessary columns
#' heatmap_plot <- DB_gg_func(data = data_frame, Oxygen_Range = 0:4, filling = "log10")
#' print(heatmap_plot)
#'
#' @export
DB_gg_func <- function(data, Oxygen_Range = c(0:6), filling = "scaled",
                       H_Count_min = -10, H_Count_Max = 10,
                       Carbon_Pos_min = 2, Carbon_Pos_Max = 20) {

  # Filter the data based on the specified delta_H range
  data <- data[data$delta_H >= H_Count_min & data$delta_H <= H_Count_Max,]

  # Filter the data based on the specified Carbon_Position range
  data <- data[data$C >= (Carbon_Pos_min - 1) & data$C <= (Carbon_Pos_Max - 1),]

  data$O_Label <- paste(data$O, "oxygen lost")

  # Create a heatmap using ggplot2
  if (filling == "scaled") {
    # If 'filling' is "scaled", use the 'scaled_intensity' column for fill
    ggplot(data[data$O %in% Oxygen_Range,], aes(y = Carbon_Position, x = delta_H, text = loss, fill = scaled_intensity)) +
      geom_tile() +
      xlab("Hydrogen difference") +
      facet_wrap(~O_Label, scales = "free")  # Create facets for each level of oxygen, allowing free scaling
  } else {
    # If 'filling' is not "scaled", use the log10 of 'Height' for fill
    ggplot(data[data$O %in% Oxygen_Range,], aes(y = Carbon_Position, x = delta_H, text = loss, fill = log10(Height))) +
      geom_tile() +
      xlab("Hydrogen difference") +
      facet_wrap(~O_Label, scales = "free")  # Create facets for each level of oxygen, allowing free scaling
  }
}
#' Generate a Plot of Carbon Counts vs Heights for Elemental Formulas
#'
#' This function creates a plot visualizing the relationship between carbon counts (C)
#' and heights of elemental formulas, conditioned on specific delta H values and other
#' parameters. It allows for filtering based on thresholds for height, delta H, and oxygen
#' levels.
#'
#' @param Data A data frame (default is `Elemental_Formula_df`) containing the elemental formula data.
#' @param column A string indicating which column in the data represents the delta H values.
#'   Default is "L_R".
#' @param Lower_H_Limit Minimum threshold for `H_Delta` to include in the plot. Default is -2.
#' @param Upper_H_Limit Maximum threshold for `H_Delta` to include in the plot. Default is 20.
#' @param Thres Minimum threshold for `Height` to include in the plot. Default is 5.
#' @param Max_C Maximum carbon count to include in the plot. Default is the maximum value in `FA_Peaks`.
#' @param Oxygen_Range A vector indicating the range of oxygen values to include. Default is 0 to 3.
#' @param point_size Size of the points in the plot. Default is 1.
#' @param alpha_v Transparency level for the lines in the plot. Default is 0.35.
#'
#' @return A ggplot object visualizing the relationship between carbon counts and heights
#'   of elemental formulas based on the specified filtering criteria. If no data meets the
#'   filtering criteria, a message is returned indicating this.
#'
#' @examples
#' # Assuming 'Elemental_Formula_df' is available and contains necessary columns
#' plot <- Plot_DB_from_FA(Data = Elemental_Formula_df, Lower_H_Limit = 0)
#' print(plot)
#'
#' @export
Plot_DB_from_FA <- function(Data = Elemental_Formula_df, column = "L_R",
                            Lower_H_Limit = -2, Upper_H_Limit = 20,
                            Thres = 5, Max_C = max(FA_DF$C),
                            Oxygen_Range = 0:3, point_size = 1,
                            alpha_v = 0.35) {

  # Rename the specified column to "H_Delta" for easier reference
  colnames(Data)[colnames(Data) == column] <- "H_Delta"

  # Filter the data based on the specified thresholds for H_Delta, Height, C, and O
  filtered_data <- Data[Data$H_Delta > Lower_H_Limit &
                          Data$H_Delta < Upper_H_Limit &
                          Data$Height > Thres &
                          Data$C < Max_C &
                          Data$O %in% Oxygen_Range, ]

  if(nrow(filtered_data) > 0) {
    filtered_data$H_Delta_label <- paste(filtered_data$H_Delta, "H", sep="")
    filtered_data$O_Label <- paste(filtered_data$O, "Oxygen", sep=" ")

    # Create the ggplot object
    Elemental_Formula_df_fig <- ggplot(filtered_data, aes(x = C, y = Height)) +
      geom_point(size = point_size) +  # Add points to the plot
      geom_line(linetype = "dashed", size = point_size, alpha = alpha_v) +  # Add dashed lines connecting the points
      facet_grid(H_Delta_label ~ O_Label, scales = "free_y") +  # Create a grid of facets based on H_Delta and O, allowing for free y-axis scaling
      scale_x_continuous(breaks = seq(1, Max_C, by = 2)) +  # Customize x-axis breaks
      xlab("Carbon Position")

    # Return the ggplot object
    Elemental_Formula_df_fig

  } else {
    "No fragments with given filtering criteria found"
  }
}

#' Extract Alternative Double Bond Series and Corresponding Plots
#'
#' This function extracts alternative double bond series and their corresponding plots
#' from a given data set of additional double bond series. It processes filtered double bonds
#' and loss heights, generating visualizations for each data entry.
#'
#' @param Data A data frame containing additional double bond series. Default is `Additional_DB_Series`.
#'
#' @return A list containing the extracted plots for starting points and double bond series,
#'   each annotated with appropriate titles based on the data names.
#'
#' @examples
#' # Assuming 'Additional_DB_Series' is available
#' result <- Extract_DB2(Data = Additional_DB_Series)
#'
#' @export
Extract_DB2 <- function(Data = Additional_DB_Series) {

  # Identify entries in the data corresponding to filtered double bonds and loss heights
  DB_List <- grepl("Sec_Data_DB_filtered", names(Data))
  FA_List <- grepl("Sec_Data_Loss_H3", names(Data))

  # Extract and plot starting points based on loss height data
  Starting_Point <- lapply(Data[FA_List], function(x) {
    Loss_Delta_H_gg <- Plot_FA_Losses(x, Oxygen_Range = 0)  # Plotting function for loss data
  })

  # Assign names to the starting points based on the data names
  names(Starting_Point) <- names(Data)[FA_List]

  # Add names to plot
  Starting_Point <- lapply(names(Starting_Point), function(x) {
    title_text <- paste("M-", str_extract(x, "(?<=\\s)\\S+(?=\\s)"), sep="")
    Starting_Point[[x]] +
      ggtitle(title_text) +
      theme(plot.title = element_text(size = 12, hjust = 0.5))
  })

  names(Starting_Point) <- names(Data)[FA_List]

  # Extract and process double bond series from the filtered data
  DB_Series <- lapply(Data[DB_List], function(x) {
    Sec_Double_Bond_Data <- Find_DB(x, delta_H = -5:8)  # Find double bonds
    Secondary_Double_Bonds_gg <- DB_gg_func(x,
                                            Oxygen_Range = 0,
                                            Carbon_Pos_Max = 40,
                                            Carbon_Pos_min = 2,
                                            H_Count_min = -10,
                                            H_Count_Max = 10)  # Generate double bond plots
  })

  Rename <- names(DB_Series)

  # Add names to plot
  DB_Series <- lapply(names(DB_Series), function(x) {
    title_text <- paste("M-", str_extract(x, "(?<=\\s)\\S+(?=\\s)"), sep="")
    DB_Series[[x]] +
      ggtitle(title_text) +
      theme(plot.title = element_text(size = 12, hjust = 0.5))
  })

  names(DB_Series) <- Rename

  # Combine starting points and double bond series results into a single output
  Output <- append(Starting_Point, DB_Series)

  # Return the combined output
  Output
}
