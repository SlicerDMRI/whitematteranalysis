#!/usr/bin/env python
#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/ipython

# Run registration on the test dataset.

import argparse
import os
import numpy
import vtk
import time

try:
    import whitematteranalysis as wma
except:
    print("<wm_register.py> Error importing white matter analysis package\n")
    raise

def main():
    #-----------------
    # Parse arguments
    #-----------------
    parser = argparse.ArgumentParser(
        description="Compare two tractography files to see how much registration has changed the points. The files must have the exact same size and lines and points, e.g. the original file and the file that results from applying a transform to a file.",
        epilog="Written by Lauren O\'Donnell, odonnell@bwh.harvard.edu.  Please reference \"O'Donnell, Lauren J., and C-F. Westin. Automatic tractography segmentation using a high-dimensional white matter atlas. Medical Imaging, IEEE Transactions on 26.11 (2007): 1562-1575.\"",
        version='1.0')
    
    parser.add_argument(
        'inputFile1',
        help='A file of whole-brain tractography as vtkPolyData (.vtk or .vtp).')
    parser.add_argument(
        'inputFile2',
        help='A file of whole-brain tractography as vtkPolyData (.vtk or .vtp).')
    
    args = parser.parse_args()
    
    pd1 = wma.io.read_polydata(args.inputFile1)
    pd2 = wma.io.read_polydata(args.inputFile2)
    
    # Go through the whole polydata and compute the differences
    # Estimate this by converting to arrays and computing distances
    # figure out hemisphere labels for this cluster
    farray1 = wma.fibers.FiberArray()
    farray1.convert_from_polydata(pd1, points_per_fiber=30)
    farray2 = wma.fibers.FiberArray()
    farray2.convert_from_polydata(pd2, points_per_fiber=30)
    
    # compute the distances between the matching points
    ddx = farray1.fiber_array_r - farray2.fiber_array_r
    ddy = farray1.fiber_array_a - farray2.fiber_array_a
    ddz = farray1.fiber_array_s - farray2.fiber_array_s
    
    dx = numpy.square(ddx)
    dy = numpy.square(ddy)
    dz = numpy.square(ddz)
    
    distance = numpy.sqrt(dx + dy + dz)
    
    # print information
    print("MIN\tMAX\tMEAN\tSTD\tMEDIAN")
    print(numpy.min(distance), '\t', numpy.max(distance), '\t', numpy.mean(distance), '\t', numpy.std(distance), '\t', numpy.median(distance))
     
            
    # figure out what percentage is less than 1mm, 2mm, 3mm, 4mm, 5mm etc.
    sz = farray1.fiber_array_r.shape
    number_of_points = float(sz[0]*sz[1])
    
    print(100 * numpy.sum(distance < 0.001) / number_of_points, 'percent < 0.001 mm')
    print(100 * numpy.sum(distance < 0.01) / number_of_points, 'percent < 0.01 mm')
    print(100 * numpy.sum(distance < 0.1) / number_of_points, 'percent < 0.1 mm')
    print(100 * numpy.sum(distance < 0.5) / number_of_points, 'percent < 0.5 mm')
    print(100 * numpy.sum(distance < 1) / number_of_points, 'percent < 1 mm')
    print(100 * numpy.sum(distance < 1.2) / number_of_points, 'percent < 1.2 mm')
    print(100 * numpy.sum(distance < 1.3) / number_of_points, 'percent < 1.3 mm')
    print(100 * numpy.sum(distance < 1.4) / number_of_points, 'percent < 1.4 mm')
    print(100 * numpy.sum(distance < 1.5) / number_of_points, 'percent < 1.5 mm')
    print(100 * numpy.sum(distance < 2) / number_of_points, 'percent < 2 mm')
    print(100 * numpy.sum(distance < 3) / number_of_points, 'percent < 3 mm')
    print(100 * numpy.sum(distance < 4) / number_of_points, 'percent < 4 mm')
    print(100 * numpy.sum(distance < 5) / number_of_points, 'percent < 5 mm')
    
    # output a new polydata showing what the distances are in case there is a pattern

if __name__ == '__main__':
    main()
