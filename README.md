# R1R2R3_data_codes
Data and python codes for "Locating the largest event observed on Mars with multi-orbit surface waves"

Python codes for reading and processing data and generating figures 1--3 are included as `figure1_rev.py`, `figure2_mod.py` and `figure3.py`, respectively. Figure 4 is a map figure produced by Geraldine Zenhausern based on the Cartopy plotting package, but the code is not included here.

The data subdirectory contains the raw and deglitched data in miniseed format (see [IRIS website](http://ds.iris.edu/ds/nodes/dmc/data/formats/miniseed/)  for more details) with instrument response removed and rotated into ZNE components, as well as a Microsoft Excel file showing all the raw picks provided by different piucking methods and the derived location info and a csv file with just the summary information sheet on the locations using the different surface wave picking methods.
