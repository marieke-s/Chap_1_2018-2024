Template from Intro2R (https://alexd106.github.io/intro2R/project_setup.html)

# Directory Structure
      Root 
           |
           |__data
           |    |_raw_data
           |    |_processed_data
           |    |_metadata
           |
           |_R
           |
           |
           |_scripts/analyses
           |
           |_output
           |
           |__figures

In my working directory I have the following directories:

Root - This is your project directory containing your .Rproj file.
  - data - I store all my data in this directory.
        - The subdirectory called raw_data contains raw data files and only raw data files. These files should be treated as read only and should not be changed in any way. If you need to process/clean/modify your data do this in R (not MS Excel) as you can document (and justify) any changes made.
        - Any processed data should be saved to a separate file and stored in the processed_data subdirectory.
        - Information about data collection methods, details of data download and any other useful metadata should be saved in a text document (see README text files below) in the metadata subdirectory.

- R - This is an optional directory where I save all of my custom R functions I have written for the current analysis. These can then be sourced into R using the source() function.

- scripts/analyses - All of the main R scripts/Rmd I have written for the current project are saved here.

- output - Outputs from my R scripts such as plots, HTML files and data summaries are saved in this directory. This helps me and my collaborators distinguish what files are outputs and which are source files.
- figures 


(R Style guide to make clean scripts : http://adv-r.had.co.nz/Style.html)
test2
