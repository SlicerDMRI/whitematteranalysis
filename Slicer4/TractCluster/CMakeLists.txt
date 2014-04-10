cmake_minimum_required(VERSION 2.8.9)

project(TractCluster)

#-----------------------------------------------------------------------------
set(EXTENSION_HOMEPAGE "http://www.slicer.org/slicerWiki/index.php/Documentation/4.3/Extensions/TractCluster")
set(EXTENSION_CATEGORY "Examples")
set(EXTENSION_CONTRIBUTORS "Lauren O'Donnell BWH")
set(EXTENSION_DESCRIPTION "This is an example of extension bundling a scripted loadable module")
set(EXTENSION_ICONURL "http://viewvc.slicer.org/viewvc.cgi/Slicer4/trunk/Extensions/Testing/TractCluster/TractCluster.png?revision=21747&view=co")
set(EXTENSION_SCREENSHOTURLS "http://wiki.slicer.org/slicerWiki/images/e/e2/Slicer-r19441-TractCluster-screenshot.png")

#-----------------------------------------------------------------------------
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

#-----------------------------------------------------------------------------
add_subdirectory(PickTracts)

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_CPACK})
