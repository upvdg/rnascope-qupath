# RNAScope QuPath
Script to detect RNAScope Spots in QuPath.

# Code Description
This script will apply a Laplacian of Gaussian filtered followed by a local maxima detector in order to locate 
RNAScope spots in fluorescent images. The script works with an arbitrary number of channels, provided as a list.
It is intended to run in QuPath[1] 0.4.3 or later. 

# Materials & Methods
RNAScope spots are identified using a Laplacian of Gaussian[2] filter (LoG) followed by a 2D local maximum finder.
In summary, annotated regions are imported from QuPath into ImageJ[3] and run through a median filter (radius=1.5 px) 
in order to remove Poisson noise.
The diffraction-limited spots are then enhanced using a LoG filter (sigma=1.0 px). The resulting image is run through 
ImageJ's local maximum finder algorithm with a user-selected channel-dependent tolerance value[4]. 
The local maxima points are then reimported as spots into QuPath for each channel.

# References
1. [Bankhead, P. et al. QuPath: Open source software for digital pathology image analysis](https://doi.org/10.1038/s41598-017-17204-5). Scientific Reports (2017)
doi:10.1038/s41598-017-17204-5
2. https://imagescience.org/meijering/software/featurej/laplacian/
3. [Schneider, C. A., Rasband, W. S., & Eliceiri, K. W. (2012). NIH Image to ImageJ: 25 years of image analysis](doi.org/10.1038/nmeth.2089). Nature Methods, 9(7), 671–675.
doi:10.1038/nmeth.2089
4. https://github.com/imagej/ImageJ/blob/master/ij/plugin/filter/MaximumFinder.java


# Inputs
A single or multichannel fluorescence image with annotations. This script will run on all annotations
The user should specify the channel names to use, as well as the prominence values for each channel.
Higher prominence values will find fewer spots, as it represents how much of a difference the peaks should have with their 
local background. 

# Outputs
After running the code, mRNA spots will be shown as Points, in the same color as the channel.
Two new measurements `RNAScope CHANNEL Spots` and `RNASctop CHANNEL Density` are also appended to each annotation.

# Dependencies
This script makes use of the ImageScience library at https://imagescience.org/meijering/software/imagescience/

# Installtion
You must download 'imagescience.jar' and place it into your QuPath extensions directory before running this script. Alternatively, you can simply drag and drop the jar file into QuPath. 
https://imagescience.org/meijering/software/imagescience/

# Running the script 
To run the script, simply drag and drop it into the QuPath interface, which will open the script editor, and hit "Run"

# Author Information
Code written by Olivier Burri, EPFL - SV - PTECH - BIOP
for Lucie Bracq, Van Der Goot Lab
Last update: 20230314
