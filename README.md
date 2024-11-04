# MsRadaR

**A Tool for De-Novo Structural Elucidation of Lipids from MS/MS Spectra (EDP-CID, EAD)**

MsRadaR contains functions to build a de-novo structural elucidation pipeline for determining double bond positions in lipid molecules from centroided MS/MS spectra. 
Elemental formula and elemental formula loss (EFL) assignment can be performed. The experimental determined EFL are compared to that of a saturated alkyl chain CH3(CH2)n
to determine the delta H difference related to the appearance of double bonds. Characteristic delta H differences are normalized and can be visualized as relative intensities
or intensity slopes. Intensity maxima peak picking is performed to spot putative double bond positions. Different m/z values can be used as starting points for CH3(CH2)n deltaH series.
Additionally, fragments related to the fatty acid composition can be visualized and indicate whether a neutral or radical loss occurred.

Additional information about interpretation rules and applications can be found within the following preprint: https://chemrxiv.org/engage/chemrxiv/article-details/6655ac6e418a5379b0787191


## Installation
To install the **MsRadaR** R package, copy the following R code and run it in your console. The `devtools` package is required and will be installed if needed.

   ```r
   if (!requireNamespace("devtools", quietly = TRUE)) {install.packages("devtools")}

   devtools::install_github("lsmsgeneva/MsRadaR")
```

## Example Files
Example files of 7 acylglycerols are provided within the "Example Files" folder. 
An example excel file for batch processing of these text files is also provided.
The folder is automaticaly installed along the R package.
The following acylglycerols are available as examples:

| Lipid                                                       | Abbreviation |
|-------------------------------------------------------------|--------------|
| TG(18:2(n-6,n-9)/18:2(n-6,n-9)/18:2(n-6,n-9))               | LLL          |
| TG(18:1(n-9)/18:1(n-9)/18:1(n-9))                           | OOO          |
| TG(18:3(n-3,n-6,n-9)/18:3(n-3,n-6,n-9)/18:3(n-3,n-6,n-9))   | LnLnLn       |
| TG(18:3(n-6,n-9,n-12)/18:3(n-6,n-9,n-12)/18:3(n-6,n-9))     | gLngLngLn    |
| TG(18:1(n-9)/16:0/18:2(n-6,n-9))                            | OPL          |
| DG(18:2(n-6,n-9)/OH/18:2(n-6,n-9))                          | 1,3LL        |
| MG(OH/20:4(n-6,n-9,n-12,n-15)/OH)                           | 2Ar          |



## Usage

The following example workflow demonstrates two different ways to use the MsRadaR package: single-file processing and batch processing.

### Example Workflow

#### 1. Single File Processing
Single file processing requires a centroided MS/MS spectrum (currently as text file).
Additional information can be inputed within given functions (e.g. precursor mass), but are not required.
```r
# Load the MsRadaR library
library(MsRadaR)

#===============================================================================
#1. Input Parameters
#===============================================================================

#Load Parameters required for processing
Load_Parameters()

#===============================================================================
#2. Data loading and preprocessing
#===============================================================================

#Specify path. e.g. for example data. replace with path of your data.
Example_Spectra_Path <- Load_Example_Path(Sample = "OPL")

#Load and filter data from a single file
Data <- Load_Data(Dat.Name = Example_Spectra_Path, Precursor.Mass = FALSE)

#===============================================================================
#3. Perform elemental formula annotation of the MS/MS spectra
#===============================================================================

Elemental_Formula_df <- Annotate_Fragments(Data)

#===============================================================================
#4. Kendricks Transformation
#===============================================================================

#Perform Kendrick transformation on the annotated elemental formulas data.frame
Kendricks_Data <- Kendricks_Transformation(Kendricks_Data = Elemental_Formula_df)

#===============================================================================
#5. Annotate Losses
#===============================================================================

#Annotate elemental formula losses based on the original input file
Fragmentation_Data <- Annotate_Losses(Data = Data)

#===============================================================================
#6. Find FA attachments
#===============================================================================

#From Elemental Formulas of FA and MG fragments
FA_Peaks <- Find_FAs(Elemental_Formula_df, Oxygen = 1:3,median_Factor = 4)
FA_DF <- Extract_FAs(FA_Peaks = FA_Peaks, Precursor_Df = Precursor_DF, Data = Elemental_Formula_df)


#From Elemental Formula losses leading to DG Fragments
FA_Loss_Peaks <- Find_FA_Losses(Fragmentation_Data)
FA_Loss_DF <- Extract_FA_Losses(Data= Fragmentation_Data,FA_Peaks = FA_Loss_Peaks, Precursor_Df = Precursor_DF )
#Find numbers of double bonds


#===============================================================================
#7. Annotate Double bonds
#===============================================================================

Double_Bond_Data <- Find_DB(Fragmentation_Data,Threeshold = Thres, delta_H = 0:9)

#===============================================================================
#7. Find DB from additional starting points
#===============================================================================

Additional_DB_Series <- Find_DB2(Data = Data)

#===============================================================================
#8. Generate Graphical Outputs
#===============================================================================
#Various graphical outputs can be generated for data visualization and interpretation
#===============================================================================
#8.1.1 Kendricks plot
#===============================================================================
Kendrick_Plot <- Plot_Kendricks_Plot(Kendricks_Data, toggle_alpha = 0, Rel.Thres = 1/10, Abs.Thres = Thres, plot_legend = F, point_size = 1)

#===============================================================================
#8.1.2 MS Plot
#===============================================================================
MS_Spec <- Plot_MS2(Kendricks_Data, Axis_Split_Size = 1, Split_Text_Size = 8, Label_Size = 8, Line_Size = 1)
#===============================================================================
#8.1.3 MS Plot with kendricks plot
#===============================================================================
MS_Kendrick_plot <- MS_Spec/Kendrick_Plot
#===============================================================================
#8.2 Fatty acid plot
#===============================================================================
#Lookg for correct function
M_FA_Plot <- Plot_FA_Losses(FA_Loss_DF)

MG_Plot <- Plot_FAs(FA_DF)
#===============================================================================
#8.3 DB Plot
#===============================================================================
Extracted_Plots <- Extract_DB_Plots(Double_Bond_Data, delta_H = "even")
#===============================================================================
#8.4 Alternative DB Fragmentation Series
#===============================================================================

Additional_DB_Series_Plot <- Extract_DB2(Data = Additional_DB_Series)

#===============================================================================
#8.5 Fatty Acid DB Frags
#===============================================================================

FA_DB_Series_PLot <- Plot_DB_from_FA(Data = Elemental_Formula_df, column = "R_L")

```

#### 2. Batch Processing
Batch processing allows to utilize all the functions mentioned above to loop through a batch of MS/MS spectra (currently, single MS/MS text spectra).
Either all MS/MS spectra within the working directory are processed or an excel (.xlsx) can be used to specificy spectra locations, names and additional 
information to consider for data processing. For instance, precursor masses, double bond series after the loss of specific elemental formulas / m/z values can be specified
within the excel file. Inputing the column name within the function will consider these losses. Setting parameters to F will prevent certain functions from beeing executed.
The results of the batch processing can be summarized within an heatmap.
```r
# Load the MsRadaR library
library(MsRadaR)

#===============================================================================
# 1. Input Parameters
#===============================================================================

# Load Parameters required for processing
Load_Parameters()

#Load Example Info for processing example files
Load_Example_Info()

#Process data
Batch_Results <- Batch_Processing(Spectra_Location=  Example_Folder, Sample_Information = Example_Info)

Batch_Summary <- Create_Summary(Data = Batch_Results, Expected_DB = "Max_DB_Col", Plot.Order = "Order")

```
## License
This project is licensed under the GPL-3 License.

