This tutorial explains how to perform whole-brain tractography parcellation of a single subject using the whitematteranalysis (WMA) software and the anatomically curated ORG white matter atlas. On this page, we provide step-by-step instructions to guide a user to run the entire tractography parcellation pipeline.

> **_Note_**: The commands used on this page outline the major steps in the wrapper script “wm_apply_ORG_atlas_to_subject.sh” in WMA. This wrapper script is provided to run the whole pipeline in one command, and it is useful for batch processing of data from multiple subjects (see here for a tutorial for batch processing).

1. Software prerequisites
   - Install the whitematteranalysis software
   - Install 3D Slicer and SlicerDMRI
    
1. Download the tutorial data
   - Download the tutorial data package (WMA_tutorial_data.zip, ~2.5GB)
   - Decompress the downlaoded zip file to *Desktop* of your computer
   - Files in the decompressed data folder should be organized as below, including:
      - An example tractography dataset (computed from one example Human Connectome Project (HCP) subject using the two-tensor Unscented Kalman Filter (UKF) fiber tracking method)
      - The ORG white matter atlas (see here for details)
      
      ![test image](../wma_small.jpg)

1. Prepare terminal environment to run related commands
    - MacOS: Open */Applications/Utilities/Terminal.app*
    - Linux (e.g. Red Hat): Open */Applications/System Tools/Terminal*
    - Windows: Open *XX*
    
      > **_Note_**: The tutorial on this page is based on MacOS. All commands in this tutorial can be directly used on Linux. For Windows, users need to change the commands by using Windows system separator “\”.

1. Go to the tutorial data folder from terminal
    - From the terminal, type the following command. (Make sure that you have decompressed the tutorial data to your desktop, as mentioned above.)
    
        ```bash
        cd /Users/YOUR_USER_NAME/Desktop/WMA_tutorial_data
        ```
    - You terminal should look like the below image (the username should change according to your computer). Type ls to list the files in the tutorial data folder on the terminal.
    
      ![test image](../wma_small.jpg)
    
1. Initial tractography data quality control

    This QC step is important to: 1) verify correct appearance of tract anatomy, and 2) verify tractography is stored using the same spatial coordinate system as the atlas tractography data.
   
   1. Run QC using “wm_quality_control_tractography.py”
   
      - This script outputs rendered images of the input tractography data (as well as other information such as fiber length distributions and diffusion data fields stored along the tracts). Here, this script to verify correct appearance of tract anatomy of the input tractography dataset. From your terminal, type the following command:
      
        ```bash
        wm_quality_control_tractography.py ./ ./QC-TractVisualization/
        ```
        
        > **_Note_**: This script also allows for tractography visualization and QC of multiple tractography datasets at the same time. This is useful when you perform a study (e.g. between-population analysis) where you have tractography data from multiple subjects. Please see here for details of batch processing using WMA.
        
      - A new “QC-TractVisualization” folder is generated, including several types of outputs:
        - Multiple HTML files to enable tractography visualization from different views. Click on one of them, e.g., “view_left.html”, the HTML file will be opened in your browser to show the tractography from the left view (as displayed below). Click on the “sample_UKF_HCP” image, another HTML page will be opened to show the tractography from 6 different views (as displayed below). These images show that the input tractography has a correct appearance of white matter tract anatomy.
        - “fiber_length_histograms.pdf” and “quality_control_fibers.txt” to show fiber length distributions. 
        - “quality_control_data.txt” to show the diffusion data fields stored along the tracts.
        
          ![test image](../wma_small.jpg)

   1. Run QC using “wm_quality_control_tract_overlap.py”
     
      - This script outputs rendered images of two input tractography datasets together to show tract overlap. Here, this script is used to check if the input tractography data and the atlas tractography data are stored using the same spatial coordinate system. From your terminal, type the following command:
       
        ```bash
        wm_quality_control_tract_overlap.py ./ORG-Atlases-1.1.1/ORG-800FC-100HCP/atlas.vtp ./example-UKF-data.vtk ./QC-TractOverlap-Initial/
        ```
      - A new “QC-TractOverlap-Initial” folder is generated, including multiple JPG files to enable visualization of tract overlap from different views. Open one of them, e.g., “view_left_tract_overlap.jpg”, where the different colors represent the different tractography data (as displayed below). This image shows that the two tractography files are in the same coordinate system, but they are not aligned together.
       
        ![test image](../wma_small.jpg)
        
        
      

        
