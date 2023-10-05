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
     -s:  Path to 3D Slicer, e.g., under macOS, it is /Applications/Slicer.app/Contents/MacOS/Slicer
Optional arguments:
     -t:  apply a transformation file to match the data to the adult brain size of the atlas.
     -r:  whole brain tractography registration mode (default = 'rig')
            rig: rigid_affine_fast : this enables a rough tractography registration. This mode in general
                                   applicable to tractography data generated from different dMRI 
                                   acquisitions and different populations (such as babies)
            nonrig: affine + nonrigid (2 stages) : this enables nonrigid deformations of the fibers. This mode
                                              needs the input tractography to be similar to the atlas tractography,
                                              e.g. two-tensor UKF tractography + HCP dMRI acquisition.
     -n:  Number of threads (default = 1). If multiple cores are available, recommended setting is 4.
          Increasing the number of threads does not speed up too much, but it requires more computational resources.) 
     -x:  If the job is being run in an environment without a X client, a virtual X server environment is needed (for 
           transforming fiber clusters in the atlas space back to the tractography space using 3D Slicer). Use value 1 
           to indicate the usage of a virtual X server (default 0).
     -d: Export diffusion tensor measurements for each fiber cluster (or fiber tract). Note that diffusion tensor 
         (or other diffusion measurements) must be stored in the input VTK file. If this is specified, -m needs to be provided.
     -m: Path to the FiberTractMeasurements module in SlicerDMRI. For example, in 3D Slicer 4.11 stable release under macOS,
         the CLI module path is:
           /Applications/Slicer.app/Contents/Extensions-28264/SlicerDMRI/lib/Slicer-4.11/cli-modules/FiberTractMeasurements
     -c Clean the internal temporary files (default : 0)
            0 : keep all files
            1 : minimal removal : initial bilateral clusters, transformed bilateral clusters
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
 Case ID (file name):         $caseID
 Input tractography:          $InputTractography
 Output directory:            $OutputDir
 ORG atlas folder:            $AtlasBaseFolder
 3D Slicer:                   $SlicerPath
 Registration mode:           $RegMode
 Number of threads:           $NumThreads
 virtual X server:            $VX
 Export dMRI measures:        $DiffMeasure
 FiberTractMeasurementsCLI:   $FiberTractMeasurementsCLI
 Clean Files:                 $CleanFiles
======================================================================================

REPORTMAPPINGPARAMETERS
}

while getopts ":hi:i:o:a:s:n:r:t:x:d:m:c:" opt; do
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
		t) TfmFile="$OPTARG"
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
	tmp=($FiberTractMeasurementsCLI)
	if [ -z "$FiberTractMeasurementsCLI" ] ; then
		echo "ERROR: -m path to FiberTractMeasurements Module must be provided."
		echo ""
		Usage
	
	# check existence of the last string in $FiberTractMeasurementsCLI
	elif [ ! -f ${tmp[-1]} ]; then
		echo "WARNING: FiberTractMeasurements could not be found, program may fail."
		echo ""
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

# Get CaseID 
fn=$(basename -- $InputTractography)
ext=${fn#*.}
caseID=${fn//.$ext/}

reportMappingParameters

echo "<wm_apply_ORG_atlas_to_subject> Working on input tractography:" $InputTractography
echo "<wm_apply_ORG_atlas_to_subject> Case ID is assigned to be the file name of the input: [" ${caseID} "]"
echo ""

# Setup output
OutputCaseFolder=$OutputDir/${caseID}
echo "<wm_apply_ORG_atlas_to_subject> Fiber clustering result will be stored at:" $OutputCaseFolder
if [ ! -d $OutputCaseFolder ]; then
	echo ' - create an output folder.' 
	mkdir -p $OutputCaseFolder
else
	
	numfiles=($OutputCaseFolder/AnatomicalTracts/T*vtp)
	numfiles=${#numfiles[@]}
	if [ $numfiles -gt 1 ] ;then
		echo ""
		echo "** Anatomical tracts ($numfiles tracts) are detected in the output folder. Manually remove all files to rerun."
		echo ""
		#exit
	fi

	echo ' - output folder exists.'

fi
echo ""

# Setup white matter parcellation atlas
RegAtlasFolder=$AtlasBaseFolder/ORG-RegAtlas-100HCP
FCAtlasFolder=$AtlasBaseFolder/ORG-800FC-100HCP
echo "<wm_apply_ORG_atlas_to_subject> White matter atlas: " $AtlasBaseFolder
echo " - tractography registration atlas:" $RegAtlasFolder
echo " - fiber clustering atlas:" $FCAtlasFolder
echo ""

# Apply transformation if provided
InTractographyDirname="$(dirname "$InputTractography")"
OutTfmTractography=$OutputCaseFolder/TransformedTracts
if [ "$TfmFile" ]; then
  if [ ! -f "$TfmFile" ]; then
    echo ""
    echo "ERROR: Transformation file not found."
    echo ""
    exit
  else
    echo "<wm_harden_transform.py> Apply transformation with file: " "$TfmFile"
    if [ $VX == 0 ]; then
      wm_harden_transform.py -t "$TfmFile" \
        "$InTractographyDirname" "$OutTfmTractography" "$SlicerPath"
    else
      xvfb-run wm_harden_transform.py -t "$TfmFile" \
        "$InTractographyDirname" "$OutTfmTractography" "$SlicerPath"
    fi
    echo ""
  fi
fi

echo "<wm_apply_ORG_atlas_to_subject> Tractography registration with mode [" $RegMode "]"
RegistrationFolder=$OutputCaseFolder/TractRegistration
if [ "$TfmFile" ]; then
  TractographyData=$(find $OutTfmTractography -type f)
else
  TractographyData=$InputTractography
fi

if [ "$RegMode" == "rig" ]; then
	RegTractography=$RegistrationFolder/${caseID}/output_tractography/${caseID}_reg.vtk
	
	if [ ! -f $RegTractography ]; then
		wm_register_to_atlas_new.py -mode rigid_affine_fast \
			$TractographyData $RegAtlasFolder/registration_atlas.vtk $RegistrationFolder/
	else
		echo " - registration has been done."
	fi
elif [ "$RegMode" == "nonrig" ]; then
	RegTractography=$RegistrationFolder/${caseID}_reg/output_tractography/${caseID}_reg_reg.vtk
	
	if [ ! -f $RegTractography ]; then
		wm_register_to_atlas_new.py -mode affine \
			$TractographyData $RegAtlasFolder/registration_atlas.vtk $RegistrationFolder/

		affineRegTract=$RegistrationFolder/${caseID}/output_tractography/${caseID}_reg.vtk
		
		wm_register_to_atlas_new.py -mode nonrigid \
			$affineRegTract $RegAtlasFolder/registration_atlas.vtk $RegistrationFolder/
	else
		echo " - registration has been done."
	fi
fi
echo ""

if [ ! -f $RegTractography ]; then
	echo ""
	echo "ERROR: Tractography registration failed. The output registered tractography data can not be found."
	echo ""
	exit
fi

# Get the case ID for fiber clustering
fn=$(basename -- $RegTractography)
ext=${fn#*.}
FCcaseID=${fn//.$ext/}

echo "<wm_apply_ORG_atlas_to_subject> Fiber clustering for whole-brain 800 fiber cluster parcellation."
FiberClusteringInitialFolder=$OutputCaseFolder/FiberClustering/InitialClusters
if [ ! -f $FiberClusteringInitialFolder/$FCcaseID/cluster_00800.vtp ]; then
	wm_cluster_from_atlas.py -j $NumThreads \
		$RegTractography $FCAtlasFolder $FiberClusteringInitialFolder -norender
else
	echo " - initial fiber clustering has been done."
fi
echo ""

numfiles=($FiberClusteringInitialFolder/$FCcaseID/*vtp)
numfiles=${#numfiles[@]}
if [ $numfiles -lt 800 ]; then
	echo ""
	echo "ERROR: Initial fiber clustering failed. There should be 800 resulting fiber clusters, but only $numfiles generated."
	echo ""
	exit
fi
 
echo "<wm_apply_ORG_atlas_to_subject> Outlier fiber removal."
FiberClusteringOutlierRemFolder=$OutputCaseFolder/FiberClustering/OutlierRemovedClusters
if [ ! -f $FiberClusteringOutlierRemFolder/${FCcaseID}_outlier_removed/cluster_00800.vtp ]; then
	wm_cluster_remove_outliers.py -j $NumThreads \
		$FiberClusteringInitialFolder/$FCcaseID $FCAtlasFolder $FiberClusteringOutlierRemFolder
else
	echo " - outlier fiber removal has been done."
fi
echo ""

numfiles=($FiberClusteringOutlierRemFolder/${FCcaseID}_outlier_removed/*vtp)
numfiles=${#numfiles[@]}
if [ $numfiles -lt 800 ]; then
	echo ""
	echo "ERROR: Outlier removal failed. There should be 800 resulting fiber clusters, but only $numfiles generated."
	echo ""
	exit
fi

echo "<wm_apply_ORG_atlas_to_subject> Hemisphere location assessment in the atlas space."
if [ ! -f $FiberClusteringOutlierRemFolder/${FCcaseID}_outlier_removed/cluster_location_by_hemisphere.log ]; then
	wm_assess_cluster_location_by_hemisphere.py \
		-clusterLocationFile $FCAtlasFolder/cluster_hemisphere_location.txt \
		$FiberClusteringOutlierRemFolder/${FCcaseID}_outlier_removed 
else
	echo " - hemisphere location assessment has been done."
fi
echo ""

if [ ! -f $FiberClusteringOutlierRemFolder/${FCcaseID}_outlier_removed/cluster_location_by_hemisphere.log ]; then
	echo ""
	echo "ERROR: Hemisphere location assessment failed. There should a cluster_location_by_hemisphere.log file, stating: \"<wm_assess_cluster_location_by_hemisphere.py> Done!!!\" "
	echo ""
	exit
fi

echo "<wm_apply_ORG_atlas_to_subject> Transform fiber clusters back to input tractography space."
FiberClustersInTractographySpace=$OutputCaseFolder/FiberClustering/TransformedClusters/${caseID}
if [ "$RegMode" == "rig" ]; then
	
	tfm=$RegistrationFolder/${caseID}/output_tractography/itk_txform_${caseID}.tfm
	
	if [ ! -f $FiberClustersInTractographySpace/cluster_00800.vtp ]; then
		if [ $VX == 0 ]; then
			wm_harden_transform.py -i -t $tfm \
				$FiberClusteringOutlierRemFolder/${FCcaseID}_outlier_removed $FiberClustersInTractographySpace "$SlicerPath"
		else
			xvfb-run wm_harden_transform.py -i -t $tfm \
				$FiberClusteringOutlierRemFolder/${FCcaseID}_outlier_removed $FiberClustersInTractographySpace "$SlicerPath"
		fi
	else
		echo " - transform has been done."
	fi
elif [ "$RegMode" == "nonrig" ]; then

	tfm=$RegistrationFolder/${caseID}_reg/output_tractography/itk_txform_${caseID}_reg.tfm
	
	if [ ! -f $FiberClustersInTractographySpace/tmp/cluster_00800.vtp ]; then
		if [ $VX == 0 ]; then
			wm_harden_transform.py -i -t $tfm \
				$FiberClusteringOutlierRemFolder/${FCcaseID}_outlier_removed $FiberClustersInTractographySpace/tmp "$SlicerPath"
		else
			xvfb-run wm_harden_transform.py -i -t $tfm \
				$FiberClusteringOutlierRemFolder/${FCcaseID}_outlier_removed $FiberClustersInTractographySpace/tmp "$SlicerPath"
		fi
	else
		echo " - transform has been done."
	fi

	tfm=$RegistrationFolder/${caseID}/output_tractography/itk_txform_${caseID}.tfm
	
	if [ ! -f $FiberClustersInTractographySpace/cluster_00800.vtp ]; then
		if [ $VX == 0 ]; then
			wm_harden_transform.py -i -t $tfm \
				$FiberClustersInTractographySpace/tmp $FiberClustersInTractographySpace "$SlicerPath"
		else
			xvfb-run wm_harden_transform.py -i -t $tfm \
				$FiberClustersInTractographySpace/tmp $FiberClustersInTractographySpace "$SlicerPath"
		fi
	else
		echo " - transform has been done."
	fi
fi
echo ""

numfiles=($FiberClustersInTractographySpace/*vtp)
numfiles=${#numfiles[@]}
if [ $numfiles -lt 800 ]; then
	echo ""
	echo "ERROR: Transforming fiber clusters failed. There should be 800 resulting fiber clusters, but only $numfiles generated."
	echo ""
	exit
fi

echo "<wm_apply_ORG_atlas_to_subject> Separate fiber clusters by hemisphere."
SeparatedClustersFolder=$OutputCaseFolder/FiberClustering/SeparatedClusters
if [ ! -f $SeparatedClustersFolder/tracts_commissural/cluster_00800.vtp ]; then

	wm_separate_clusters_by_hemisphere.py $FiberClustersInTractographySpace $SeparatedClustersFolder

else
	echo " - separation has been done."
fi
echo ""

numfiles=($SeparatedClustersFolder/tracts_commissural/*vtp)
numfiles=${#numfiles[@]}
if [ $numfiles -lt 800 ]; then
	echo ""
	echo "ERROR: Separating fiber clusters failed. There should be 800 resulting fiber clusters in each folder, but only $numfiles generated."
	echo ""
	exit
fi

# Apply the inverse transformation if provided
ClusterDirnames=($(find "$SeparatedClustersFolder" -mindepth 1 -type d | sort))
OutInvTfmTractographyDirnameRoot=$OutputCaseFolder/InvTransformedTracts
if [ "$TfmFile" ]; then
  if [ ! -f "$TfmFile" ]; then
    echo ""
    echo "ERROR: Transformation file not found."
    echo ""
    exit
  else
    echo "<wm_harden_transform.py> Apply inverse transformation with file: " "$TfmFile"
    mkdir $OutInvTfmTractographyDirnameRoot
    for ClusterDirname in "${ClusterDirnames[@]}"; do
      GroupName=($(basename $ClusterDirname))
      OutInvTfmTractographyDirname=$OutInvTfmTractographyDirnameRoot/$GroupName
      mkdir $OutInvTfmTractographyDirname
      if [ $VX == 0 ]; then
        wm_harden_transform.py -i -t "$TfmFile" \
          "$ClusterDirname" "$OutInvTfmTractographyDirname" "$SlicerPath"
      else
        xvfb-run wm_harden_transform.py -i -t "$TfmFile" \
          "$ClusterDirname" "$OutInvTfmTractographyDirname" "$SlicerPath"
      fi
    done
    echo ""
  fi
fi

echo "<wm_apply_ORG_atlas_to_subject> Append clusters into anatomical tracts."
AnatomicalTractsFolder=$OutputCaseFolder/AnatomicalTracts
if [ "$TfmFile" ]; then
  TractographyData=$OutInvTfmTractographyDirnameRoot
else
  TractographyData=$SeparatedClustersFolder
fi

if [ ! -f $AnatomicalTractsFolder/T_UF_right.vtp ]; then
	
	wm_append_clusters_to_anatomical_tracts.py $TractographyData/ $FCAtlasFolder $AnatomicalTractsFolder

else
	echo " - Appending clusters into anatomical tracts has been done."
fi
echo ""

numfiles=($AnatomicalTractsFolder/*vtp)
numfiles=${#numfiles[@]}
if [ $numfiles -lt 73 ]; then
	echo ""
	echo "ERROR: Appending clusters into anatomical tracts failed. There should be 73 resulting fiber clusters, but only $numfiles generated."
	echo ""
	exit
fi

if [ $DiffMeasure == 1 ]; then

  os=$(uname)

  if  [[ "$os" == 'Linux' ]]; then
    measurement_cli_cmd=$SlicerPath" --launch "$FiberTractMeasurementsCLI
  else
    measurement_cli_cmd=$FiberTractMeasurementsCLI
  fi

	echo "<wm_apply_ORG_atlas_to_subject> Report diffusion measurements of fiber clusters."
	if [ ! -f $SeparatedClustersFolder/diffusion_measurements_commissural.csv ]; then
		wm_diffusion_measurements.py \
			$SeparatedClustersFolder/tracts_commissural $SeparatedClustersFolder/diffusion_measurements_commissural.csv "$measurement_cli_cmd"
	else
		echo " - diffusion measurements of commissural clusters has been done."
	fi
	if [ ! -f $SeparatedClustersFolder/diffusion_measurements_left_hemisphere.csv ]; then
		wm_diffusion_measurements.py \
			$SeparatedClustersFolder/tracts_left_hemisphere $SeparatedClustersFolder/diffusion_measurements_left_hemisphere.csv "$measurement_cli_cmd"
	else
		echo " - diffusion measurements of left hemisphere clusters has been done."
	fi
	if [ ! -f $SeparatedClustersFolder/diffusion_measurements_right_hemisphere.csv ]; then
		wm_diffusion_measurements.py \
			$SeparatedClustersFolder/tracts_right_hemisphere $SeparatedClustersFolder/diffusion_measurements_right_hemisphere.csv "$measurement_cli_cmd"
	else
		echo " - diffusion measurements of right hemisphere clusters has been done."
	fi

	if [ ! -f $SeparatedClustersFolder/diffusion_measurements_right_hemisphere.csv ]; then
		echo ""
		echo "ERROR: Reporting diffusion measurements of fiber clusters. failed. No diffusion measurement (.csv) files generated."
		echo ""
		exit
	fi

fi
echo ""


if [ $DiffMeasure == 1 ]; then

  os=$(uname)

  if  [[ "$os" == 'Linux' ]]; then
    measurement_cli_cmd=$SlicerPath" --launch "$FiberTractMeasurementsCLI
  else
    measurement_cli_cmd=$FiberTractMeasurementsCLI
  fi

	echo "<wm_apply_ORG_atlas_to_subject> Report diffusion measurements of the anatomical tracts."
	if [ ! -f $AnatomicalTractsFolder/diffusion_measurements_anatomical_tracts.csv ]; then
		wm_diffusion_measurements.py \
			$AnatomicalTractsFolder $AnatomicalTractsFolder/diffusion_measurements_anatomical_tracts.csv "$measurement_cli_cmd"
	else
		echo " - diffusion measurements of anatomical tracts has been done."
	fi

	if [ ! -f $AnatomicalTractsFolder/diffusion_measurements_anatomical_tracts.csv ]; then
		echo ""
		echo "ERROR: Reporting diffusion measurements of fiber clusters. failed. No diffusion measurement (.csv) files generated."
		echo ""
		exit
	fi

fi
echo ""

if [ $CleanFiles == 1 ]; then
	echo "<wm_apply_ORG_atlas_to_subject> Clean files using minimal removal."
	rm -rf $OutputCaseFolder/FiberClustering/InitialClusters
	rm -rf $OutputCaseFolder/FiberClustering/TransformedClusters

elif [ $CleanFiles == 2 ]; then
	echo "<wm_apply_ORG_atlas_to_subject> Clean files using maximal removal."
	rm -rf $OutputCaseFolder/TractRegistration/*/output_tractography/*vtk
	rm -rf $OutputCaseFolder/TractRegistration/*/iteration*
	rm -rf $OutputCaseFolder/FiberClustering/InitialClusters/*
	rm -rf $OutputCaseFolder/FiberClustering/OutlierRemovedClusters/*
	rm -rf $OutputCaseFolder/FiberClustering/TransformedClusters/*
fi

exit
 
