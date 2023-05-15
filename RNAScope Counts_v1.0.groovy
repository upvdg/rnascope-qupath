/* = CODE DESCRIPTION =
 * This script will apply a Laplacian of Gaussian filtered followed by a local maxima detector in order to locate 
 * RNAScope spots in fluorescent images. The script works with an arbitrary number of channels, provided as a list.
 * 
 * == MATERIALS & METHODS == 
 * RNAScope spots are identified using a Laplacian of Gaussian[1] filter (LoG) followed by a 2D local maximum finder.
 * In summary, annotated regions are imported from QuPath into ImageJ[2] and run through a median filter (radius=1.5 px) 
 * in order to remove Poisson noise.
 * The diffraction-limited spots are then enhanced using a LoG filter (sigma=1.0 px). The resulting image is run through 
 * ImageJ's local maximum finder algorithm with a channel-dependent tolerance value[3]. 
 * The local maxima points are then reimported into QuPath for each channel.
 * 
 * REFERENCES
 * ----------
 * [1] https://imagescience.org/meijering/software/featurej/laplacian/
 * [2] Schneider, C. A., Rasband, W. S., & Eliceiri, K. W. (2012). NIH Image to ImageJ: 25 years of image analysis. Nature Methods, 9(7), 671–675. doi:10.1038/nmeth.2089
 * [3] https://github.com/imagej/ImageJ/blob/master/ij/plugin/filter/MaximumFinder.java
 *
 * == INPUTS ==
 * A single or multichannel fluorescence image with annotations. This script will run on all annotations
 * The user should specify the channel names to use, as well as the prominence values for each channel.
 * Higher prominence values will find fewer spots, as it represents how much of a difference the peaks should have with their 
 * local background. 
 *
 * == OUTPUTS ==
 * After running the code, mRNA spots will be shown Points, in the color of the channel.
 * Two new measurements "RNAScope CHANNEL Spots" and RNASctop CHANNEL Density are also appended to each annotation
 *
 * = DEPENDENCIES =
 * This script makes use of the ImageScience library at https://imagescience.org/meijering/software/imagescience/
 *
 * = INSTALLATION =
 * You must download 'imagescience.jar' and place it into your QuPath extensions directory
 * https://imagescience.org/meijering/software/imagescience/
 *
 * = AUTHOR INFORMATION =
 * Code written by Olivier Burri, EPFL - SV - PTECH - BIOP
 * for Lucie Bracq, Van Der Goot Lab
 * Last update: 20230314
 *
 * = COPYRIGHT =
 * © All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, BioImaging and Optics Platform (BIOP), 2023
 *
 * Licensed under the BSD-3-Clause License:
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided
 * that the following conditions are met:
 * 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer
 *    in the documentation and/or other materials provided with the distribution.
 * 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products
 *    derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
 * BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

// USER SETTINGS

// Set prominence and channels
def rnaScopeChannelNames = ["CY5"]
def prominencePerChannel = [80]


// START OF SCRIPT

// Need pixel calibration for subsequent steps
px = getCurrentServer().getPixelCalibration().getAveragedPixelSizeMicrons()

// This could be a parameter for the user, but as they are diffraction-limited, we pick the smallest sigma
def laplacianSigma = prominencePerChannel.collect{ px }

def annotations = getAnnotationObjects()

// Get the channel names
def channelNames = getCurrentServer().getMetadata().getChannels().collect{ it.getName() }

// Clean previous detections on the image
clearDetections()

annotations.eachWithIndex{ annotation, i ->
    
    println "Processing $annotation - (${i+1}/${annotations.size()})"
    
    // Get the image into ImageJ, use a downsample of 1 to keep the max resolution
    def request = RegionRequest.createInstance( getCurrentServerPath(), 1, annotation.getROI() )
    def pathImage = IJTools.convertToImagePlus( getCurrentServer(), request )
    
    def allChannels = ChannelSplitter.split( pathImage.getImage() )
    
    [rnaScopeChannelNames, prominencePerChannel, laplacianSigma].transpose().each{ channel, prominence, lapSigma -> 
        
        // Find the channel
        def channelIndex = channelNames.indexOf( channel )
        def image = allChannels[ channelIndex ]
        
        // Remove noisy pixels
        IJ.run( image, "Median...", "radius=1.5" )
        
        // Use channel color
        def color = getCurrentServer().getChannel( channelIndex ).getColor()
        
        // Make LoG and get the ImageJ spots
        def spots_roi = getSpots( image, lapSigma, prominence )
    
        image.close()
        
        // convert it to QuPath ROI
        def pathroi = IJTools.convertToROI( spots_roi, pathImage )
        
        def detection = PathObjects.createDetectionObject( pathroi )
        def detectionClass = getDerivedPathClass( annotation.getPathClass(), "RNAScope "+channel )
        detectionClass.setColor( color )
        detection.setPathClass( detectionClass )
        
        annotation.addPathObject( detection )
        
        // Add statistics
        def area = annotation.getROI().getArea() / px / px
        annotation.measurements["RNAScope "+channel+" Spots"] = spots_roi.size() as double
        annotation.measurements["RNAScope "+channel+" Density"] = spots_roi.size() / area
    }
}

fireHierarchyUpdate()

// END OF SCRIPT

// Helper method to get spots roi using a LoG filter
def getSpots( def image, def sigma, def prominence ) {
    
    def laplace = new Laplacian()
    
    // ImageScience needs the output to be defined
    def output =  new FloatImage( Image.wrap( image ) )
    
    // Compute LoG
    output = laplace.run( output, sigma )
    def lap = output.imageplus()
    lap.setRoi( image.getRoi() )
    
    // Finally get the local maxima
    IJ.run( lap, "Find Maxima...", "prominence=$prominence light output=[Point Selection]" )
        
    def roi = lap.getRoi()
    return roi
}

// IMPORTS
import qupath.lib.roi.*
import ch.epfl.biop.qupath.utils.*
import imagescience.feature.Laplacian
import imagescience.image.Image
import imagescience.image.FloatImage
import ij.IJ
import qupath.imagej.objects.*
import qupath.lib.objects.*
import ij.plugin.ChannelSplitter