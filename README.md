# MsRadaR

**A Tool for De-Novo Structural Elucidation of Lipids from MS/MS Spectra (EDP-CID, EAD)**

MsRadaR contains functions to build a de-novo structural elucidation pipeline for determining double bond positions in lipid molecules from centroided MS/MS spectra. 
Elemental formula and elemental formula loss (EFL) assignment can be performed. The experimental determined EFL are compared to that of a saturated alkyl chain CH3(CH2)n
to determine the delta H difference related to the appearance of double bonds. Characteristic delta H differences are normalized and can be visualized as relative intensities
or intensity slopes. Intensity maxima peak picking is performed to spot putative double bond positions. Different m/z values can be used as starting points for CH3(CH2)n deltaH series.
Additionally, fragments related to the fatty acid composition can be visualized and indicate whether a neutral or radical loss occurred.

## Installation

To install the **MsRadaR** package, you can use the `devtools` package, which provides functions to install packages hosted on GitHub.

Make sure you have the `devtools` package installed. If you haven't installed it yet, you can do so by running:

   ```r
   install.packages("devtools")

   devtools::install_github("lsmsgeneva/MsRadaR")
```
## Usage

The following example demonstrates two typical workflows for using the **MsRadaR** package: single file processing and batch processing.
Example text files are provided within the directory. An example excel file for batch processing of these text files is also provided.

### Example Workflow

#### 1. Single File Processing
Single file processing requires a centroided MS/MS spectrum (currently as text file).
Additional information can be inputed within given functions (e.g. precursor mass), but are not required.
```r
# Load the MsRadaR library
library(MsRadaR)

#===============================================================================
# 1. Input Parameters
#===============================================================================

# Load Parameters required for processing
Load_Parameters()

#===============================================================================
# 2. Data Loading and Preprocessing
#===============================================================================

# Load and filter data from a single file
Data <- Load_Data(Dat.Name = "Lipid_Feature.txt", Precursor.Mass = FALSE)

#===============================================================================
# 3. Perform Elemental Formula Annotation of the MS/MS Spectra
#===============================================================================

Elemental_Formula_df <- Annotate_Fragments(Data)

#===============================================================================
# 4. Kendricks Transformation
#===============================================================================

# Perform Kendrick transformation on the annotated elemental formulas data frame
Kendricks_Data <- Kendricks_Transformation(Kendricks_Data = Elemental_Formula_df)

#===============================================================================
# 5. Annotate Losses
#===============================================================================

# Annotate elemental formula losses based on the original input file
Fragmentation_Data <- Annotate_Losses(Data = Data)

#===============================================================================
# 6. Find Fatty Acid Attachments
#===============================================================================

# From Elemental Formulas of FA and MG fragments
FA_Peaks <- Find_FAs(Elemental_Formula_df, Oxygen = 1:3)
FA_DF <- Extract_FAs(FA_Peaks = FA_Peaks, Precursor_Df = Precursor_DF, Data = Elemental_Formula_df)

# From Elemental Formula losses leading to DG Fragments
FA_Loss_Peaks <- Find_FA_Losses(Fragmentation_Data)
FA_Loss_DF <- Extract_FA_Losses(Data = Fragmentation_Data, FA_Peaks = FA_Loss_Peaks, Precursor_Df = Precursor_DF)

#===============================================================================
# 7. Annotate Double Bonds
#===============================================================================

Double_Bond_Data <- Find_DB(Fragmentation_Data, Threshold = Thres, delta_H = 0:9)

#===============================================================================
# 8. Find Double Bonds from Additional Starting Points
#===============================================================================

Additional_DB_Series <- Find_DB2(Data = Data)

#===============================================================================
# 9. Generate Graphical Outputs
#===============================================================================

# 9.1. Kendrick Plot
Kendrick_Plot <- Plot_Kendricks_Plot(Kendricks_Data, toggle_alpha = 0, Rel.Thres = 1/10, Abs.Thres = Thres, plot_legend = FALSE, point_size = 1)

# 9.2. MS Plot
MS_Spec <- Plot_MS2(Kendricks_Data, Axis_Split_Size = 1, Split_Text_Size = 8, Label_Size = 8, Line_Size = 1)

# 9.3. MS Plot with Kendrick Plot
MS_Kendrick_plot <- MS_Spec / Kendrick_Plot

# 9.4. Fatty Acid Plot
M_FA_Plot <- Plot_FA_Losses(FA_Loss_DF)

# 9.5. Double Bond Plot
Extracted_Plots <- Extract_DB_Plots(Double_Bond_Data, delta_H = "even")

# 9.6. Alternative Double Bond Fragmentation Series
Additional_DB_Series_Plot <- Extract_DB2(Data = Additional_DB_Series)

# 9.7. Fatty Acid Double Bond Fragments
FA_DB_Series_Plot <- Plot_DB_from_FA(Data = Elemental_Formula_df, column = "R_L")

```

#### 2. Batch Processing
Batch processing allows to utilize all the functions mentioned above to loop through a batch of MS/MS spectra (currently, single MS/MS text spectra).
Either all MS/MS spectra within the working directory are processed or an excel (.xlsx) can be used to specificy spectra locations, names and additional 
information to consider for data processing. For instance, precursor masses, double bond series after the loss of specific elemental formulas / m/z values can be specified
within the excel file. Inputing the column name within the function will consider these losses. Setting parameters to F will prevent certain functions from beeing executed.
```r
# Load the MsRadaR library
library(MsRadaR)

#===============================================================================
# 1. Input Parameters
#===============================================================================

# Load Parameters required for processing
Load_Parameters()

#===============================================================================
# 2. Batch Processing
#===============================================================================

Batch_Results <- Batch_Processing(Sample_Information = "Spectra_Informations.xlsx")

#===============================================================================
# 3. Create Summary
#===============================================================================

Batch_Summary <- Create_Summary(Data = Batch_Results, Expected_DB = "Max_DB_Col", Plot.Order = "Order")

```
## License
This project is licensed under the GPL-3 License.

