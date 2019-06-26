#!/bin/bash

function Usage {
    cat <<USAGE
Usage:
`basename $0` -i InputTractography -o OutputDirectory -a ORGAtlasFolder -s PathTo3DSlicer
Compulsory arguments:
     -i:  Input tractography data stored in VTK (.vtk or .vtp). Note: fiber streamlines need to be in the RAS coordinates.
     -o:  Output directory to save all fiber clustering outputs.
     -a:  Path to the ORG atlas (the anatomically curated atlas ORG-800FC-100HCP), within which there should be 
          two folders named: ORG-RegAtlas-100HCP and ORG-800FC-100HCP 
     -s:  Path to 3D Slicer, e.g., under MacOS, it is /Applications/Slicer.app/Contents/MacOS/Slicer
Optional arguments:
     -r:  whole brain tractography registration mode (default = 'rig')
            rig: rigid_affine_fast : this enables a rough tractography registraion. This mode in general 
                                   applicable to tractography data generated from different dMRI 
                                   acquisions and different populations (such as babies)
            nonrig: affine + nonrigid (2 stages) : this enables nonrigid deformations of the fibers. This mode
                                              needs the input tractography to be similar to the atals tractography, 
                                              e.g. two-tensor UKF tractography + HCP dMRI acquisiton. 
     -n:  Number of threads (default = 1). If mutiple cores are available, recommended setting is 4. 
          Increasing the number of threads does not speed up too much, but it requires more computational resources.) 
     -x:  If the job is being run in an environment without a X client, a virtual X server environment is needed (for 
           transforming fiber clusters in the atlas space back to the tractography space using 3D Slicer). Use value 1 
           to indicate the usage of a virtual X server (default 0).
     -d: Export diffusion tensor measurements for each fiber cluster (or fiber tract). Note that diffusion tensor 
         (or other diffusion measurements) must be stored in the input VTK file. If this is specified, -m needs to be provided.
     -m: Path to the FiberTractMeasurements module in SlicerDMRI. For example, in 3D Slicer 4.11 stable release under MacOS, 
         the CLI module path is:
           /Applications/Slicer.app/Contents/Extensions-28264/SlicerDMRI/lib/Slicer-4.11/cli-modules/FiberTractMeasurements
     -c Clean the internal temporary files (default : 0)
            0 : keep all files
            1 : mininal removal : initial bilateral clusters, transformed bilateral clusters
            2 : maximal removal : remove registration temp files, initial bilateral clusters, outlier removed bilateral clusters, transformed bilateral clusters

Example:
`basename $0` -i UKF-tract.vtk -o ./FiberClusteringTest -a ./ORG-Atlases-1.1.1 -s /Applications/Slicer.app/Contents/MacOS/Slicer
--------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------
script by Fan Zhang (fzhang@bwh.harvard.edu)
USAGE
    exit 1
}

function reportMappingParameters {
    cat <<REPORTMAPPINGPARAMETERS
======================================================================================
 Mapping parameters
--------------------------------------------------------------------------------------
 Input tractography:          $InputTractography
 Output directory:            $OutputDir
 ORG atlas folder:            $AtlasBaseFolder
 3D Slier:                    $SlicerPath
 Registraion mode:            $RegMode
 Number of threads:           $NumThreads
 virtual X server:            $VX
 Export dMRI measures:        $DiffMeasure
 FiberTractMeasurementsCLI:   $FiberTractMeasurementsCLI
 Clean Files:                 $CleanFiles
======================================================================================

REPORTMAPPINGPARAMETERS
}

while getopts ":hi:i:o:a:s:n:r:x:d:m:c:" opt; do
  case $opt in
  	h) Usage; exit 0
    ;;
    i) InputTractography="$OPTARG"
    ;;
    o) OutputDir="$OPTARG"
    ;;
    a) AtlasBaseFolder="$OPTARG"
    ;;
    s) SlicerPath="$OPTARG"
    ;;
    r) RegMode="$OPTARG"
	;;
    n) NumThreads="$OPTARG"
	;;
    x) VX="$OPTARG"
    ;;
    d) DiffMeasure="$OPTARG"
    ;;
    m) FiberTractMeasurementsCLI="$OPTARG"
    ;;
    c) CleanFiles="$OPTARG"
    ;;
    \?) echo "\nERROR: Invalid option -$OPTARG"; echo ""; Usage
    ;;
  esac
done

if ! [[ ! "$RegMode" ]] ; then
	if ! [[ "$RegMode" == "rig" || "$RegMode" == "nonrig" ]]; then
		echo "ERROR: invalid tractography registration mode" $RegMode
		echo ""
		Usage
	fi
else
	RegMode='rig'
fi

if [[ $NumThreads -lt 1 ]]; then
    NumThreads=1
fi

if [ -z "$VX" ] ; then
	VX=0
else
	if [[ $VX -lt 1 ]]; then
	    VX=0
	else
		VX=1
	fi
fi

if [ -z "$DiffMeasure" ] ; then
	DiffMeasure=0
	FiberTractMeasurementsCLI=None
else
	if [[ $DiffMeasure -lt 1 ]]; then
	    DiffMeasure=0
	else
		DiffMeasure=1
	fi
fi

if [ $DiffMeasure == 1 ]; then
	if [ -z "$FiberTractMeasurementsCLI" ] ; then
		echo "ERROR: -m path to FiberTractMeasurements Module must be provided."
		echo ""
		Usage
	elif [ ! -f $FiberTractMeasurementsCLI ]; then
		echo "ERROR: FiberTractMeasurements Module does not exist."
		echo ""
		Usage
	fi
fi

if [ -z "$CleanFiles" ] ; then
	CleanFiles=0
else
	if ! [[ $CleanFiles -lt 3 ]]; then
	    echo "ERROR: invalid clean file mode."
		echo ""
		Usage
	fi
fi

reportMappingParameters

# Get CaseID 
ext=${InputTractography#*.}
fn=$(basename -- $InputTractography)
caseID=${fn//.$ext/}

echo "<WMA_batch> Working on input tractography:" $InputTractography
echo "<WMA_batch> Case ID is assigned to be the file name of the input: [" ${caseID} "]"
echo ""

# Setup output
OutputCaseFolder=$OutputDir/${caseID}
echo "<WMA_batch> Fiber clustering result will be stored at:" $OutputCaseFolder
if [ ! -d $OutputCaseFolder ]; then
	echo ' - create an output folder.' 
	mkdir -p $OutputCaseFolder
else
	
	numfiles=($OutputCaseFolder/AnatomicalTracts/T*vtp)
	numfiles=${#numfiles[@]}
	if [ $numfiles -gt 1 ] ;then
		echo ""
		echo "** Final anantomical tracts ($numfiles tracts) exist in the output folder. Manually remove all files to rerun."
		echo ""
		exit
	fi
	echo ' - output folder exists.'

fi
echo ""

# Setup white matter parcellation atlas
RegAtlasFolder=$AtlasBaseFolder/ORG-RegAtlas-100HCP
FCAtlasFolder=$AtlasBaseFolder/ORG-800FC-100HCP
echo "<WMA_batch> White matter atlas: " $AtlasBaseFolder
echo " - tractography registration atlas:" $RegAtlasFolder
echo " - fiber clustering atlas:" $FCAtlasFolder
echo ""


echo "<WMA_batch> Tractography registraion with mode [" $RegMode "]"
RegistrationFolder=$OutputCaseFolder/TractRegistration
if [ "$RegMode" == "rig" ]; then
	RegTractography=$RegistrationFolder/${caseID}/output_tractography/${caseID}_reg.vtk
	
	if [ ! -f $RegTractography ]; then
		wm_register_to_atlas_new.py -l 40 -mode rigid_affine_fast \
			$InputTractography $RegAtlasFolder/registration_atlas.vtk $RegistrationFolder/
	else
		echo " - registration has been done."
	fi
elif [ "$RegMode" == "nonrig" ]; then
	RegTractography=$RegistrationFolder/${caseID}_reg/output_tractography/${caseID}_reg_reg.vtk
	
	if [ ! -f $RegTractography ]; then
		wm_register_to_atlas_new.py -l 40 -mode affine \
			$InputTractography $RegAtlasFolder/registration_atlas.vtk $RegistrationFolder/

		affineRegTract=$RegistrationFolder/${caseID}/output_tractography/${caseID}_reg.vtk
		
		wm_register_to_atlas_new.py -l 40 -mode nonrigid \
			$affineRegTract $RegAtlasFolder/registration_atlas.vtk $RegistrationFolder/
	else
		echo " - registration has been done."
	fi
fi
echo ""

echo $RegTractography
# Get the case ID for fiber clustering
ext=${RegTractography#*.}
fn=$(basename -- $RegTractography)
FCcaseID=${fn//.$ext/}

echo "<WMA_batch> Fiber clustering for whole-brain 800 fiber cluster parcellation."
FiberClusteringInitialFolder=$OutputCaseFolder/FiberClustering/InitialClusters
if [ ! -f $FiberClusteringInitialFolder/$FCcaseID/cluster_00800.vtp ]; then
	wm_cluster_from_atlas.py -l 40 -j $NumThreads \
		$RegTractography $FCAtlasFolder $FiberClusteringInitialFolder
else
	echo " - initial fiber clustering has been done."
fi
echo ""

echo "<WMA_batch> Outlier fiber removal."
FiberClusteringOutlierRemFolder=$OutputCaseFolder/FiberClustering/OutlierRemovedClusters
if [ ! -f $FiberClusteringOutlierRemFolder/${FCcaseID}_outlier_removed/cluster_00800.vtp ]; then
	wm_cluster_remove_outliers.py \
		$FiberClusteringInitialFolder/$FCcaseID $FCAtlasFolder $FiberClusteringOutlierRemFolder
else
	echo " - outlier fiber removal has been done."
fi
echo ""

echo "<WMA_batch> Hemisphere location assessment in the atlas space."
if [ ! -f $FiberClusteringOutlierRemFolder/${FCcaseID}_outlier_removed/cluster_location_by_hemisphere.log ]; then
	wm_assess_cluster_location_by_hemisphere.py \
		-clusterLocationFile $FCAtlasFolder/cluster_hemisphere_location.txt \
		$FiberClusteringOutlierRemFolder/${FCcaseID}_outlier_removed 
else
	echo " - hemisphere location assessment has been done."
fi
echo ""

echo "<WMA_batch> Transform fiber clusters back to input tractography space."
FiberClustersInTractographySpace=$OutputCaseFolder/FiberClustering/TransformedClusters/${caseID}
if [ "$RegMode" == "rig" ]; then
	
	tfm=$RegistrationFolder/${caseID}/output_tractography/itk_txform_${caseID}.tfm
	
	if [ ! -f $FiberClustersInTractographySpace/cluster_00800.vtp ]; then
		if [ $VX == 0 ]; then
			wm_harden_transform.py -i -t $tfm \
				$FiberClusteringOutlierRemFolder/${FCcaseID}_outlier_removed $FiberClustersInTractographySpace $SlicerPath
		else
			xvfb-run wm_harden_transform.py -i -t $tfm \
				$FiberClusteringOutlierRemFolder/${FCcaseID}_outlier_removed $FiberClustersInTractographySpace $SlicerPath
		fi
	else
		echo " - transform has been done."
	fi
elif [ "$RegMode" == "nonrig" ]; then

	tfm=$RegistrationFolder/${caseID}_reg/output_tractography/itk_txform_${caseID}_reg.tfm
	
	if [ ! -f $FiberClustersInTractographySpace/tmp/cluster_00800.vtp ]; then
		if [ $VX == 0 ]; then
			wm_harden_transform.py -i -t $tfm \
				$FiberClusteringOutlierRemFolder/${FCcaseID}_outlier_removed $FiberClustersInTractographySpace/tmp $SlicerPath
		else
			xvfb-run wm_harden_transform.py -i -t $tfm \
				$FiberClusteringOutlierRemFolder/${FCcaseID}_outlier_removed $FiberClustersInTractographySpace/tmp $SlicerPath
		fi
	else
		echo " - transform has been done."
	fi

	tfm=$RegistrationFolder/${caseID}/output_tractography/itk_txform_${caseID}.tfm
	
	if [ ! -f $FiberClustersInTractographySpace/cluster_00800.vtp ]; then
		if [ $VX == 0 ]; then
			wm_harden_transform.py -i -t $tfm \
				$FiberClustersInTractographySpace/tmp $FiberClustersInTractographySpace $SlicerPath
		else
			xvfb-run wm_harden_transform.py -i -t $tfm \
				$FiberClustersInTractographySpace/tmp $FiberClustersInTractographySpace $SlicerPath
		fi
	else
		echo " - transform has been done."
	fi
fi
echo ""

echo "<WMA_batch> Separate fiber clusters by hemisphere."
SeparatedClustersFolder=$OutputCaseFolder/FiberClustering/SeparatedClusters
if [ ! -f $SeparatedClustersFolder/tracts_commissural/cluster_00800.vtp ]; then

	wm_separate_clusters_by_hemisphere.py $FiberClustersInTractographySpace $SeparatedClustersFolder

else
	echo " - separation has been done."
fi
echo ""

if [ $DiffMeasure == 1 ]; then

	echo "<WMA_batch> Report diffusion measurements of the fiber clusters."
	if [ ! -f $SeparatedClustersFolder/tracts_commissural/DiffusionMeasurements.csv ]; then
		wm_diffusion_measurements.py \
			$SeparatedClustersFolder/tracts_commissural $SeparatedClustersFolder/tracts_commissural $FiberTractMeasurementsCLI
	else
		echo " - diffusion measurements of commissural clusters has been done."
	fi
	if [ ! -f $SeparatedClustersFolder/tracts_left_hemisphere/DiffusionMeasurements.csv ]; then
		wm_diffusion_measurements.py \
			$SeparatedClustersFolder/tracts_left_hemisphere $SeparatedClustersFolder/tracts_left_hemisphere $FiberTractMeasurementsCLI
	else
		echo " - diffusion measurements of left hemisphere clusters has been done."
	fi
	if [ ! -f $SeparatedClustersFolder/tracts_right_hemisphere/DiffusionMeasurements.csv ]; then
		wm_diffusion_measurements.py \
			$SeparatedClustersFolder/tracts_right_hemisphere $SeparatedClustersFolder/tracts_right_hemisphere $FiberTractMeasurementsCLI
	else
		echo " - diffusion measurements of right hemisphere clusters has been done."
	fi
fi
echo ""

echo "<WMA_batch> Append clusters into anatomical tracts."

hemispheric_tracts=(T_AF T_CB T_MdLF T_PLIC T_SLF-I T_SLF-II T_SLF-III T_EC T_EmC T_ICP T_ILF T_IOFF T_UF \
	T_Intra-CBLM-I\&P T_Intra-CBLM-PaT \
	T_CPC T_CR-F T_CR-P T_CST \
	T_SF T_SO T_SP T_TF T_TO T_TP\
	T_Sup-F T_Sup-FP T_Sup-O T_Sup-OT T_Sup-P T_Sup-PO T_Sup-PT T_Sup-T )

commissural_tracts=(T_CC1 T_CC2 T_CC3 T_CC4 T_CC5 T_CC6 T_CC7 T_MCP )

AnatomicalTractsFolder=$OutputCaseFolder/AnatomicalTracts

for T_mrml in "${hemispheric_tracts[@]}"
do 
	tract_name=$(basename $T_mrml)
	# echo " - Appending tract:" $tract_name
	if [ ! -f $AnatomicalTractsFolder/${tract_name}_left.vtp ]; then
		wm_append_clusters.py $SeparatedClustersFolder/tracts_left_hemisphere $AnatomicalTractsFolder \
			-appendedTractName ${tract_name}_left -tractMRML $FCAtlasFolder/${tract_name}.mrml
	fi
	if [ ! -f $AnatomicalTractsFolder/${tract_name}_right.vtp ]; then
		wm_append_clusters.py $SeparatedClustersFolder/tracts_right_hemisphere $AnatomicalTractsFolder \
			-appendedTractName ${tract_name}_right -tractMRML $FCAtlasFolder/${tract_name}.mrml
	fi
done

for T_mrml in "${commissural_tracts[@]}"
do 
	tract_name=$(basename $T_mrml)
	# echo " - Appending tract:" $tract_name
	if [ ! -f $AnatomicalTractsFolder/${tract_name}.vtp ]; then
		wm_append_clusters.py $SeparatedClustersFolder/tracts_commissural $AnatomicalTractsFolder \
			-appendedTractName ${tract_name} -tractMRML $FCAtlasFolder/${tract_name}.mrml
	fi
done
echo ""

if [ $DiffMeasure == 1 ]; then

	echo "<WMA_batch> Report diffusion measurements of the anatomical tracts."
	if [ ! -f $AnatomicalTractsFolder/DiffusionMeasurements.csv ]; then
		wm_diffusion_measurements.py \
			$AnatomicalTractsFolder $AnatomicalTractsFolder $FiberTractMeasurementsCLI
	else
		echo " - diffusion measurements of anatomical tracts has been done."
	fi
fi
echo ""

if [ $CleanFiles == 1 ]; then
	echo "<WMA_batch> Clean files using mininal removal."
	rm -rf $OutputCaseFolder/FiberClustering/InitialClusters
	rm -rf $OutputCaseFolder/FiberClustering/TransformedClusters

elif [ $CleanFiles == 2 ]; then
	echo "<WMA_batch> Clean files using maximal removal."
	rm -rf $OutputCaseFolder/TractRegistration/*/output_tractography/*vtk
	rm -rf $OutputCaseFolder/TractRegistration/*/iteration*
	rm -rf $OutputCaseFolder/FiberClustering/InitialClusters/*
	rm -rf $OutputCaseFolder/FiberClustering/OutlierRemovedClusters/*
	rm -rf $OutputCaseFolder/FiberClustering/TransformedClusters/*
fi

