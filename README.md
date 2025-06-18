# SPLain
SPLain is a Shiny app designed to interactively explore and visualize the results of Differential Transcript Usage (DTU) analyses. 


## Installation

### Clone the repository
To run the SPLain app, you need to clone the repository to your local machine. You can do this using the following command in your terminal or command prompt:

```bash
   git clone https://github.com/beamimc/SPLain.git
```

### Install dependencies
To ensure that all the necessary dependencies for the SPLain app are installed, you can use the `renv` package along with `BiocManager`. This will help manage the R package environment and install any missing packages.

1. Install the `renv` and `BiocManager` packages if you haven't already:
```r
install.packages("renv","BiocManager") # if not installed
```
2. To install all of SPLain’s dependencies, first set your working directory to the root of the cloned SPLain repository, then run the following R code:
```r
library(renv)
library(BiocManager)

# 1. Make sure you’re in the SPLain directory
setwd("~/SPLain")
getwd()                 # should show ".../SPLain"

# 2. Check dependencies that are not installed in your R library
# and install them using BiocManager
pkg_list <- renv::dependencies("app")$Package
pkgs_already_installed <- find.package(pkg_list, quiet=TRUE) |> basename()

pkgs_to_install <- setdiff(pkg_list, pkgs_already_installed)
BiocManager::install(pkgs_to_install)
```

## Data Requirements
SPLain requires two main data inputs: the Differential Transcript Usage (DTU) results and exon definitions.
1. The DTU results should be in a SummarizedExperiment format, output from running satuRn.

2. The exon is a GRangesList object containing exon definitions for all transcripts present in the SummarizeExperiment object. 
Currently, SPLain only supports sequence data for the human genome (hg38) via the `BSgenome.Hsapiens.UCSC.hg38` package.

We provide a sample of a formated dataset in the `data` directory that can be used to test the app. A tutorial on how to generate these files is available in https://txomics.github.io/tapir/testing.html. Here, its explained how to run the `satuRn` pipeline to generate the DTU results and how to create the exon definitions from a GTF file.

### Data Format

## Usage
To run the app with the sample data, we provide the example code in run_app.R file at the repository root. This file contains the necessary commands needed to:

1. Set the working directory to the SPLain root  
2. Load the example DTU results and exon definitions from `data/`  
3. Source the main app code in `app/`  
4. Start the Shiny application

Once you’re in the SPLain directory, in the R terminal run:

```r
source("run_app.R")
```

```r
# run_app.R
# 1. Set your working directory to the root of the SPLain repo
setwd(here::here())

# 2. Load required libraries
library(here)
library(shiny)

# 3. Read in the example DTU results and exon definitions
se    <- readRDS(here::here("data", "glinos_saturn_dtu.rds")) # you may change this to your own SummarizedExperiment object
exons <- readRDS(here::here("data", "glinos_exons.rds")) # you may change this to your own GRangesList object with exon definitions
#    (Default genome: hg38 via BSgenome.Hsapiens.UCSC.hg38)

# 4. Point to the app directory and source the main app code
#    (Default genome: hg38 via BSgenome.Hsapiens.UCSC.hg38)
app_dir <- here::here("app")
source(file.path(app_dir, "app.R"))

# 5. Run the Shiny app with your data and app directory
shiny::runApp(
  app(
    se      = se,
    exons   = exons,
    app_dir = app_dir
  )
)
```
This will start the Shiny app, and you can access it in your web browser in the automatically opened tab. 
You can also open the app manually by navigating to shinys output that will be printed in the R console, which will look like this:
Once you run `source("run_app.R")`, the Shiny app will launch and automatically open in your default web browser.  

If it doesn’t open automatically, look for a message in the R console and copy the URL provided:

```r
Listening on http://127.0.0.1:xxxx 
```



## Feedback
We would love to hear your feedback and welcome requests to SPLain! If you have suggestions for improvements or new features request, please  post an Issue on GitHub.

## License
SPLain is licensed under the MIT License. See the [LICENSE](LICENSE) file for more details.

