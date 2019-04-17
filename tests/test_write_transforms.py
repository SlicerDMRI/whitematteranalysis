from __future__ import print_function
import whitematteranalysis as wma
import vtk
import numpy
import os

outdir = 'test_write_transforms'
fname = 'test_transform_input/brain_0001.vtk'
slicer_command = "/Applications/Slicer.app/Contents/MacOS/Slicer"

# make output directory and read test data
if not os.path.isdir(outdir):
    os.mkdir(outdir)
pd = wma.io.read_polydata(fname)

print("================================================================================")
print("Creating synthetic transforms, saving to disk, and applying to polydatas.")
print("Output directory is:", outdir)
print("Input polydata is:", fname)
print("================================================================================")

# create synthetic affine transform and save it
affine_trans = [0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0]
vtktrans = wma.register_two_subjects.convert_transform_to_vtk(affine_trans)
transform_list = []
transform_list.append(vtktrans)
subject_ids = ['test_affine']
wma.io.write_transforms_to_itk_format(transform_list, outdir, subject_ids)

# transform the polydata directly using this transform
transformer = vtk.vtkTransformPolyDataFilter()
if (vtk.vtkVersion().GetVTKMajorVersion() >= 6.0):
    transformer.SetInputData(pd)
else:
    transformer.SetInput(pd)
transformer.SetTransform(vtktrans)
transformer.Update()
wma.io.write_polydata(transformer.GetOutput(),  os.path.join(outdir, 'test_affine.vtk'))
del transformer

# create synthetic bspline transform and save it
res = 8
#bspline_trans = numpy.zeros(res*res*res*3)
# move all points by the vector 10,10,10
#bspline_trans = bspline_trans + 10.0
# random movement of all grid points up to 10mm
print("================================================================================")
print("Creating random displacement field up to 10mm for bspline test.")
print("================================================================================")
bspline_trans = numpy.random.random(res*res*res*3) * 10

vtktrans = wma.register_two_subjects_nonrigid_bsplines.convert_transform_to_vtk(bspline_trans)
transform_list = []
transform_list.append(vtktrans)
subject_ids = ['test_bspline']
wma.io.write_transforms_to_itk_format(transform_list, outdir, subject_ids)

# transform the polydata directly using this transform
transformer = vtk.vtkTransformPolyDataFilter()
if (vtk.vtkVersion().GetVTKMajorVersion() >= 6.0):
    transformer.SetInputData(pd)
else:
    transformer.SetInput(pd)
transformer.SetTransform(vtktrans)
transformer.Update()
wma.io.write_polydata(transformer.GetOutput(),  os.path.join(outdir, 'test_bspline.vtk'))


# transform the polydata using Slicer and compare
print("\n================================================================================")
print("Run these commands to transform the polydata using Slicer and compare:")
print("================================================================================")
print("wm_harden_transform.py test_transform_input -t test_write_transforms/itk_txform_test_affine.tfm test_affine /Applications/Slicer.app/Contents/MacOS/Slicer")
print("python ../bin-test/wm_compare_vtks.py test_affine/brain_0001_trans.vtk test_write_transforms/test_affine.vtk")

print("wm_harden_transform.py test_transform_input -t test_write_transforms/itk_txform_test_bspline.tfm test_bspline /Applications/Slicer.app/Contents/MacOS/Slicer")
print("python ../bin-test/wm_compare_vtks.py test_bspline/brain_0001_trans.vtk test_write_transforms/test_bspline.vtk")

# Actually run the bspline test. This depends on Slicer path being correct,
# and also on this code being run from within the test directory.
print("\n================================================================================")
print("Testing B Spline transform application with Slicer (wm_harden_transform.py).")
print("================================================================================")
if os.path.isdir("test_bspline"):
    if os.path.exists("test_bspline/brain_0001_trans.vtk"):
        os.remove("test_bspline/brain_0001_trans.vtk")
    os.rmdir("test_bspline")
os.system("wm_harden_transform.py test_transform_input -t test_write_transforms/itk_txform_test_bspline.tfm test_bspline /Applications/Slicer.app/Contents/MacOS/Slicer")
print("\n================================================================================")
print("Comparing transformed polydatas (wm_compare_vtks.py):")
print("python-vtk RAS forward transform versus saved ITK IJK inverse transform, applied using Slicer:")
print("================================================================================")
os.system("python ../bin-test/wm_compare_vtks.py test_bspline/brain_0001_trans.vtk test_write_transforms/test_bspline.vtk")
