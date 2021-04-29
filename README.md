# GeoPEARL_Parameteriser
Code to parameterise GeoPEARL using GIS data

MarsFileExaminer.py
The MarsFileExaminer file was created to read, extract and export to an input file MARS climate data. The class, upon initialisation, 
reads a file (in this case a MARS climate file) and stores the information in a pandas dataframe. This data is then clipped to the region, 
and dates, of interest. The data is formatted to appear in the correct columns in the output, then added to a text file with the extension 
.met with a predefined header.

SoilProfileCreator.py
The SoilProfileCreator file contains the Soil_profile_creator class which stores the soil and met station information. The class reads 
information from the following files as generated above:

SoilSTUtoProfileName.csv:        Links between the soil profile names and the STUs
ProfileName_to_int.csv:        Link between soil profile names and the int value

In addition the following file is read by the class from SPADE14:

FR_EST_HOR.xls:        Soil properties from the different horizons

The mode of operation required is also supplied to the class - this determines which ptf rules to use. Currently there are two options: 
‘HYPRES’ or ‘Toth’ which use the ptfs they refer to.

Once the data have been imported, the class methods to create the soil and plot information can be invoked. The sol_creator() method 
determines which ptfs to use, and then proceeds to create the relevant data structures relevant for outputting as a file for the soil information. 
The sol_export() method exports the relevant data, including a predefined header, to a text file ending in the .sol extension. The plo_creator() 
method similarly creates a data structure containing all the information required to create the plot information. The plo_export() method when 
invoked writes this information, again including a predefined header, to a text file with the .plo extension. These files are suitable for use 
as inputs into the GeoPEARL model.

Finally, the map_export() method can be used to output an asc file that contains the soil profile for each raster cell in the study area. This 
is used when creating the maps using GeoPEARL for import into a GIS.
