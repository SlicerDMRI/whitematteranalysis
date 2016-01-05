
import os
import argparse
import vtk
import numpy

def read_keypoint_file(filename):
    """ Read Matt Toews keypoint file """

    ## # featExtract 1.1
    ## # Extraction Voxel Resolution (ijk) : 176 208 176
    ## # Extraction Voxel Size (mm)  (ijk) : 1.000000 1.000000 1.000000
    ## # Feature Coordinate Space: millimeters (gto_xyz)
    ## Features: 829
    ## Scale-space location[x y z scale] orientation[o11 o12 o13 o21 o22 o23 o31 o32 o32] 2nd moment eigenvalues[e1 e2 e3] info flag[i1] descriptor[d1 .. d64]
    
    print filename

    # ignore first 6 rows for now and assume the data is as listed above
    location = []
    scale = []
    orientation1 = []
    orientation2 = []
    orientation3 = []
    eigenvalues = []
    flag = []
    descriptor = []
    f = open(filename, 'r')
    line_number = 0
    feature_count = 0
    for line in f:
        if line_number < 6:
            print line
        else:
            #print line_number
            feature_count += 1
            #print "FEATURE:", line
            # convert line to floating point list of numbers
            data = numpy.fromstring(line, dtype=float, sep='\t')
            #print data.shape
            location.append(data[0:3])
            scale.append(data[3:4])
            orientation1.append(data[4:7])
            orientation2.append(data[7:10])
            orientation3.append(data[10:13])
            eigenvalues.append(data[13:16])
            flag.append(data[16:17])
            descriptor.append(data[17:81])

        line_number += 1
    print "Read in", feature_count, "features."
    #print line
    return location, scale, orientation1, orientation2, orientation3, eigenvalues, flag, descriptor
    
def convert_keypoint_to_polydata(location, scale, orientation1, orientation2, orientation3, eigenvalues, flag, descriptor):

    number_of_keypoints = len(location)
    print number_of_keypoints
    
    # set up vtk objects
    pd = vtk.vtkPolyData()
    points = vtk.vtkPoints()
    points.SetNumberOfPoints(number_of_keypoints)
    pd.SetPoints(points)
    point_data = pd.GetPointData()
    
    array_scale = vtk.vtkFloatArray()
    array_scale.SetNumberOfComponents(1)
    array_scale.SetNumberOfTuples(number_of_keypoints)
    array_scale.SetName('scale')
    
    array_orientation1 = vtk.vtkFloatArray()
    array_orientation1.SetNumberOfComponents(3)
    array_orientation1.SetNumberOfTuples(number_of_keypoints)    
    array_orientation1.SetName('orientation1')

    array_orientation2 = vtk.vtkFloatArray()
    array_orientation2.SetNumberOfComponents(3)
    array_orientation2.SetNumberOfTuples(number_of_keypoints)    
    array_orientation2.SetName('orientation2')

    array_orientation3 = vtk.vtkFloatArray()
    array_orientation3.SetNumberOfComponents(3)
    array_orientation3.SetNumberOfTuples(number_of_keypoints)    
    array_orientation3.SetName('orientation3')

    array_eigenvalues = vtk.vtkFloatArray()
    array_eigenvalues.SetNumberOfComponents(3)
    array_eigenvalues.SetNumberOfTuples(number_of_keypoints)
    array_eigenvalues.SetName('eigenvalues')
    array_flag = vtk.vtkIntArray()
    array_flag.SetNumberOfComponents(1)
    array_flag.SetNumberOfTuples(number_of_keypoints)    
    array_flag.SetName('flag')
    array_descriptor = vtk.vtkFloatArray()
    array_descriptor.SetNumberOfComponents(64)
    array_descriptor.SetNumberOfTuples(number_of_keypoints)
    array_descriptor.SetName('descriptor')
    
    point_data.AddArray(array_scale)
    point_data.AddArray(array_orientation1)
    point_data.AddArray(array_orientation2)
    point_data.AddArray(array_orientation3)
    point_data.AddArray(array_eigenvalues)
    point_data.AddArray(array_flag)
    point_data.AddArray(array_descriptor)
    
    # fill data into vtk object
    point_idx = 0
    for [loc, sc, ori1, ori2, ori3, eig, fl, desc] in zip(location, scale, orientation1, orientation2, orientation3, eigenvalues, flag, descriptor):
        ## print loc
        ## print sc
        ## print ori
        ## print eig
        ## print fl
        ## print desc
        
        points.SetPoint(point_idx, loc)
        array_scale.SetTuple(point_idx, sc)
        array_orientation1.SetTuple(point_idx, ori1)
        array_orientation2.SetTuple(point_idx, ori2)
        array_orientation3.SetTuple(point_idx, ori3)
        array_eigenvalues.SetTuple(point_idx, eig)
        array_flag.SetTuple(point_idx, fl)
        array_descriptor.SetTuple(point_idx, desc)
        
        point_idx += 1
    return pd
    
#-----------------
# Parse arguments
#-----------------
parser = argparse.ArgumentParser(
    description="Convert key feature file to vtk.",
    epilog="Written by Lauren O\'Donnell, odonnell@bwh.harvard.edu",
    version='1.0')

parser.add_argument(
    'inputFile',
    help='Contains feature keypoint data in Matt Toews format.')

args = parser.parse_args()

subject_id = os.path.splitext(os.path.basename(args.inputFile))[0]
output_filename = subject_id + '.vtk'
print "Converting keypoint to polydata:\n", args.inputFile, "==>>", output_filename

location, scale, orientation1, orientation2, orientation3, eigenvalues, flag, descriptor = read_keypoint_file(args.inputFile)

pd = convert_keypoint_to_polydata(location, scale, orientation1, orientation2, orientation3, eigenvalues, flag, descriptor)

print "Polydata has", pd.GetNumberOfPoints(), "keypoints."

#print pd

# Save the polydata
writer = vtk.vtkPolyDataWriter()
writer.SetInput(pd)
writer.SetFileName(output_filename)
writer.Write()
del writer

print "Saved file:", output_filename

# Also make a file to visualize to test
source1 = vtk.vtkSphereSource()
source1.SetRadius(1)
source1.Update()
source2 = vtk.vtkSphereSource()
source2.SetRadius(1.2)
source2.Update()
appender = vtk.vtkAppendPolyData()
#appender.AddInputData(source1.GetOutput())
#appender.AddInputData(source2.GetOutput())
appender.AddInput(source1.GetOutput())
appender.AddInput(source2.GetOutput())
appender.Update()
# This scalar will be used for the scaling of the glyphs
pd.GetPointData().SetActiveScalars('scale')
# glyph filter to create the glyphs at each point
glyph = vtk.vtkGlyph3D()
glyph.ScalingOn()
#glyph.SetInputData(pd)
#glyph.SetSourceConnection(0, appender.GetOutputPort())
glyph.SetInput(pd)
glyph.SetSource(appender.GetOutput())
glyph.Update()
pd_glyph = glyph.GetOutput()
print pd_glyph
# Tell the polydata to treat the orientation as vectors and create another polydata for those
pd.GetPointData().SetActiveVectors('orientation1')
glyph1 = vtk.vtkGlyph3D()
glyph1.SetVectorModeToUseVector()
glyph1.Update()
sourcec1 = vtk.vtkCylinderSource()
# make it larger than the sphere
sourcec1.SetHeight(3)
glyph1.SetInput(pd)
glyph1.SetSource(sourcec1.GetOutput())
# turn scaling off to just see orientation
#glyph1.ScalingOff()
pd_glyph1 = glyph1.GetOutput()
print pd_glyph1
# By duplicating the above code it's possible to make a polydata for the other vectors too.

# Save the polydata
output_filename_glyph = subject_id + 'glyphs.vtk'
writer = vtk.vtkPolyDataWriter()
writer.SetInput(pd_glyph)
writer.SetFileName(output_filename_glyph)
writer.Write()


# Save the polydata
output_filename_glyph1 = subject_id + 'glyphs1.vtk'
writer = vtk.vtkPolyDataWriter()
writer.SetInput(pd_glyph1)
writer.SetFileName(output_filename_glyph1)
writer.Write()





