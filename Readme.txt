Created by Jeffrey D. Munzar
Published February 22, 2017

- - - - - -

TO RUN AN ACE SCANNING ANALYSIS

The raw microarray data can be analyzed using Matlab by editing lines at the top of the following m-file and running the script:
ACE_Scanning_Analyzer_v1_0.m

1. Importantly, the working directory (the directory this file is in) must be added to the matlab path. This is provided as one of the first lines in the script.

2. The initialization of some script variables is done by importing a .mat database that has been created for each ACE scanning microarray experiment. This is done by uncommenting the desired line of code at the start of the script.

3. Additional functions and software bundles required to view the heat maps are included in the <Code> folder. Inkscape can be used to label and truncate the Matlab-generated heat maps.

4. Raw data for each ACE scanning dataset is found in subfolders in the <RawDataFiles> folder.


Note 1: The ACE scanning analysis will save figures and analyzed results as dataset-specific subfolders under the <CompiledResults> folder.

Note 2: A script to generate ACEs varying in length, location and/or incorporating single mismatches to an aptamer sequence is included in the <Code> folder. This file can be used to design microarrays for ACE scanning:
MYcroArray_Generator.m


- - - - - -

Anyone is free to use or modify this software for personal, non-commercial use, so long as reference is made back to the current work and accompanying publication.