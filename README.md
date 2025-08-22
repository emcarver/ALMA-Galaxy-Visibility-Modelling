# ALMA Visibility Modelling
This repository contains code and other documentation I wrote during the course of my co-op at HAA over the summer of 2025. My project allowed me to learn methods of image synthesis, and compare that to analysis performed directly with the data by testing visibility modelling on galaxies observed by ALMA. 

When visibilities are imaged for visualization and analysis, the process of imaging introduces artifacts such as the resolution limit imposed by the synthesized beam. Rather than bringing visibilities into the image plane for analysis, it is possible to bring a model into the visibility space and perform analysis there. Visibility modelling techniques can allow for sub-beam structure to be recovered, without the sensitivity trade-off that comes with increased resolution when imaging. Visibility modelling has been utilized in the study of protoplanetary disks, where it has been used to investigate disk features such as gaps, rings, and inner cavities in great detail. However, this same technique has not been applied to structures in other systems, such as galaxies. 

This code was written to test visibility modelling on the starburst ring at the center of the nearby galaxy NGC 3351. See Sun et al., 2024, *ApJ*, **967**, 133 for the data used in this example and for an example of how the measurements provided by visibility modelling can allow us to study the nature of star formation. 

This repo contains two folders, which contain code used in two different implementations of visibility modelling. 

The radial_version folder contains code used to fit axisymmetric models defined by a radial brightness distribution. This method was tested with observations of the inner star-forming ring in NGC 3351, and produced fit models as well as products useful in the visualization of those fits. A README file in that folder goes into greater detail on that tool. 

The 2d_model_version folder contains code designed to fit models defined in a two-dimensional plane, thus allowing a greater freedom in the complexity of the model. During my co-op term I was not able to get this version of the code running as expected, but it is provided to show generally what this use case would look like. A README file in that folder goes into greater detail on what is included. 
