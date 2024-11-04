# Citations
This project uses the following R packages. Proper citations for each package are provided to acknowledge the authors.
In addition, licenses for each package are listed. The peakwindow function from the cardidates package has been modified
as doucmented. 

## Rdisop
  License: GPL-2

  Sebastian Boecker, Matthias Letzel, Zsuzsanna Liptak and Anton Pervukhin SIRIUS: Decomposing isotope patterns for
  metabolite identification. Bioinformatics, 25(2):218-224, 2009.

  Sebastian Boecker, Zsuzsanna Liptak, Marcel Martin, Anton Pervukhin and Henner Sudek DECOMP---from interpreting Mass
  Spectrometry peaks to solving the Money Changing Problem. Bioinformatics, 24(4):591-593, 2008.

  Sebastian Boecker, Matthias Letzel, Zsuzsanna Liptak and Anton Pervukhin Decomposing metabolomic isotope patterns. In
  Proc. of Workshop on Algorithms in Bioinformatics (WABI 2006), volume 4175 of Lect. Notes Comput. Sci., pages 12-23.
  Springer, 2006.

  Sebastian Boecker and Zsuzsanna Liptak A fast and simple algorithm for the Money Changing Problem. Algorithmica,
  48(4):413-432, 2007.

## CHNOSZ
  License: GPL-3

  Dick JM (2019). “CHNOSZ: Thermodynamic calculations and diagrams for geochemistry.” _Frontiers in Earth Science_, *7*,
  180. doi:10.3389/feart.2019.00180 <https://doi.org/10.3389/feart.2019.00180>.

  Dick JM (2021). “Diagrams with multiple metals in CHNOSZ.” _Applied Computing and Geosciences_, *10*, 100059.
  doi:10.1016/j.acags.2021.100059 <https://doi.org/10.1016/j.acags.2021.100059>.

  Dick JM (2008). “Calculation of the relative metastabilities of proteins using the CHNOSZ software package.”
  _Geochemical Transactions_, *9*, 10. doi:10.1186/1467-4866-9-10 <https://doi.org/10.1186/1467-4866-9-10>.

## dplyr
  License: MIT

  Wickham H, François R, Henry L, Müller K, Vaughan D (2023). _dplyr: A Grammar of Data Manipulation_. R package version
  1.1.3, <https://CRAN.R-project.org/package=dplyr>.

## ggplot2
  License: MIT
  
  H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.

## scales
  License: MIT
  
  Wickham H, Pedersen T, Seidel D (2023). _scales: Scale Functions for Visualization_. R package version 1.3.0,
  <https://CRAN.R-project.org/package=scales>.

## stringr
  License: MIT
  
  Wickham H (2023). _stringr: Simple, Consistent Wrappers for Common String Operations_. R package version 1.5.1,
  <https://CRAN.R-project.org/package=stringr>.

## readxl
  License: MIT
  
  Wickham H, Bryan J (2023). _readxl: Read Excel Files_. R package version 1.4.3,
  <https://CRAN.R-project.org/package=readxl>.

## rlist
  License: MIT
  
  Ren K (2021). _rlist: A Toolbox for Non-Tabular Data Manipulation_. R package version 0.4.6.2,
  <https://CRAN.R-project.org/package=rlist>

## tidyr
  License: MIT
  
  Wickham H, Vaughan D, Girlich M (2024). _tidyr: Tidy Messy Data_. R package version 1.3.1,
  <https://CRAN.R-project.org/package=tidyr>.

## ggrepel
  License: GPL-3
  
  Slowikowski K (2024). _ggrepel: Automatically Position Non-Overlapping Text Labels with 'ggplot2'_. R package version
  0.9.5, <https://CRAN.R-project.org/package=ggrepel>.

## patchwork
  License: MIT
  
  Pedersen T (2024). _patchwork: The Composer of Plots_. R package version 1.3.0,
  <https://CRAN.R-project.org/package=patchwork>.

## pastecs
  License: GPL-2
  
  Grosjean P, Ibanez F (2024). _pastecs: Package for Analysis of Space-Time Ecological Series_. R package version 1.4.2,
  <https://CRAN.R-project.org/package=pastecs>.

## cardidates
  Licenses: GPL-2 and GPL-3
  
  Rolinski, S., Horn, H., Petzoldt, T., Paul, L. (2007). Identification of cardinal dates in phytoplankton time series to
  enable the analysis of long-term trends. Oecologia 153, 997--1008. doi:10.1007/s00442-007-0783-2

### Modifications

- **Function Modified:** `peakwindow`  
  *Description of Changes:*  
  - The `peakwindow` function from the `cardidates` package was adapted to ignore errors steming from incomplete data series.  
  - Changes are documented within the function `peakwindow_adapt` within this R package. 
    This new function prevents errors steming from insufficient datapoints.
    As insufficient datapoints prevent the evaluation of an intensity maxima, no intensity maxima will be printed.
