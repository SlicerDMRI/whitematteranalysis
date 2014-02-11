""" io.py

This module provides input of vtk polydata tractography files (vtk/vtp).
It also provides a class for I/O of laterality results.

read_polydata

Function to read vtkPolyData in .vtk or .vtp form

write_laterality_results

Function to write laterality indices, histograms, polydata to summarize
laterality output

read_laterality_results

This function reads in the laterality data for further analysis.

"""

import os
import pickle
import glob

import numpy
import vtk
try:
    import matplotlib.pyplot
    USE_MATPLOTLIB = 1
except ImportError:
    USE_MATPLOTLIB = 0
    print "<io.py> Failed to import matplotlib, cannot save histograms."
    print "<io.py> Please install matplotlib for this functionality."


import render
import filter

VERBOSE = 0


def read_polydata(filename):
    """Read whole-brain tractography as vtkPolyData format."""

    if VERBOSE:
        print "Reading in data from", filename, "..."

    basename, extension = os.path.splitext(filename)

    if   (extension == '.vtk'):
        reader = vtk.vtkPolyDataReader()
    elif (extension == '.vtp'):
        reader = vtk.vtkXMLPolyDataReader()
    else:
        print 'Cannot recognize model file format'
        return None

    reader.SetFileName(filename)
    reader.Update()
    outpd = reader.GetOutput()
    del reader
    if VERBOSE:
        print "Done reading in data from", filename
        print "Number of lines found:", outpd.GetNumberOfLines()

    return outpd

def list_vtk_files(input_dir):
    # Find input files
    input_mask = "{0}/*.vtk".format(input_dir)
    input_mask2 = "{0}/*.vtp".format(input_dir)
    input_pd_fnames = glob.glob(input_mask) + glob.glob(input_mask2)
    return(input_pd_fnames)
    
def read_and_preprocess_polydata_directory(input_dir, fiber_length, number_of_fibers):
    """ Find and read all .vtk and .vtp files in the given directory
    input_dir. Preprocess with fiber length threshold and downsample
    to desired number of fibers."""
    
    input_pd_fnames = list_vtk_files(input_dir)
    num_pd = len(input_pd_fnames)
    
    print "<io.py> ======================================="
    print "<io.py> Reading vtk and vtp files from directory: ", input_dir
    print "<io.py> Total number of files found: ", num_pd
    print "<io.py> ======================================="

    input_pds = list()
    subject_ids = list()
    sidx = 0

    for fname in input_pd_fnames:
        subject_id = os.path.splitext(os.path.basename(fname))[0]
        subject_ids.append(subject_id)
        print "<io.py>  ", sidx + 1, "/",  num_pd, subject_id, " Reading ", fname, "..."
        pd = read_polydata(fname)
        print "<io.py>  ", sidx + 1, "/",  num_pd, subject_id, " Input number of fibers:", pd.GetNumberOfLines()
        pd2 = filter.preprocess(pd, fiber_length)
        print "<io.py>  ", sidx + 1, "/",  num_pd, subject_id, " Length threshold", fiber_length, "mm. Number of fibers retained:", pd2.GetNumberOfLines()
        pd3 = filter.downsample(pd2, number_of_fibers)
        print "<io.py>  ", sidx + 1, "/",  num_pd, subject_id, " Downsample to", number_of_fibers, "fibers. Number of fibers retained:", pd3.GetNumberOfLines()        
        input_pds.append(pd3)
        sidx += 1
        print "<io.py> ======================================="

    print "<io.py> ======================================="
    print "<io.py> Done reading vtk and vtp files from directory: ", input_dir
    print "<io.py> Total number of files read: ", len(input_pds)
    print "<io.py> ======================================="
        
    return input_pds, subject_ids

    
                            
def write_polydata(polydata, filename):
    """Write polydata as vtkPolyData format, according to extension."""

    if VERBOSE:
        print "Writing ", filename, "..."

    basename, extension = os.path.splitext(filename)

    if   (extension == '.vtk'):
        writer = vtk.vtkPolyDataWriter()
        writer.SetFileTypeToBinary()
    elif (extension == '.vtp'):
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetDataModeToBinary()
    else:
        print 'Cannot recognize model file format'
        return None

    writer.SetFileName(filename)
    if (vtk.vtkVersion().GetVTKMajorVersion() >= 6.0):
        writer.SetInputData(polydata)
    else:
        writer.SetInput(polydata)
    writer.Update()

    del writer

    if VERBOSE:
        print "Done writing ", filename
        print "Number of lines found:", outpd.GetNumberOfLines()


class LateralityResults:

    """Results of laterality computation for a subject.

    This class defines the structure returned by ComputeWhiteMatterLaterality.
    I/O functions included in this class are read and write.

    """

    def __init__(self):
        # I/O parameters
        self.directory = ''
        # results data storage
        self.polydata = None
        self.laterality_index = None
        self.right_hem_distance = None
        self.left_hem_distance = None
        # computation parameters
        self.sigma = None
        self.points_per_fiber = None
        self.threshold = None
        self.left_hem_similarity = None
        self.right_hem_similarity = None
        self.hemisphere = None
        
    def write(self, dirname, savedist=False):
        """Write output laterality results for one subject."""
        print "a"
        if not os.path.isdir(dirname):
            os.mkdir(dirname)
        print "b"
        # output polydata
        writer = vtk.vtkPolyDataWriter()
        if (vtk.vtkVersion().GetVTKMajorVersion() >= 6.0):
            writer.SetInputData(self.polydata)
        else:
            writer.SetInput(self.polydata)
        
        writer.SetFileName(os.path.join(dirname, 'tractography_with_LI.vtk'))
        writer.Write()
        print "c"
        # output LI and other data values to text file
        # first pickle everything for later python processing
        fid = open(os.path.join(dirname, 'pickle_laterality_index.txt'), 'w')
        pickle.dump(self.laterality_index, fid)
        fid.close()

        fid = open(os.path.join(dirname, 'pickle_left_hem_similarity.txt'), 'w')
        pickle.dump(self.left_hem_similarity, fid)
        fid.close()
        fid = open(os.path.join(dirname, 'pickle_right_hem_similarity.txt'), 'w')
        pickle.dump(self.right_hem_similarity, fid)
        fid.close()
        fid = open(os.path.join(dirname, 'pickle_hemisphere.txt'), 'w')
        pickle.dump(self.hemisphere, fid)
        fid.close()

        if savedist:
            fid = open(os.path.join(dirname, 'pickle_right_hem_distance.txt'), 'w')
            pickle.dump(self.right_hem_distance, fid)
            fid.close()
            fid = open(os.path.join(dirname, 'pickle_left_hem_distance.txt'), 'w')
            pickle.dump(self.left_hem_distance, fid)
            fid.close()
            print "d"
        # now output human-readable LI values
        fid = open(os.path.join(dirname, 'laterality_index_values.txt'), 'w')
        for idx in range(0, len(self.laterality_index)):
            fid.write(str(self.laterality_index[idx]))
            fid.write('\n')
        fid.close()
        print "e"
        # generate histogram (needs matplotlib)
        li_stats = self.laterality_index[numpy.nonzero(self.laterality_index)]
        if USE_MATPLOTLIB and 0:
            try:
                print "f"
                matplotlib.pyplot.hist(li_stats, bins=25, color=[.6, .6, .6])
                print "f1"
                matplotlib.pyplot.savefig(
                    os.path.join(dirname, 'LI_histogram.pdf'))
                print "f2"
                matplotlib.pyplot.close()
                print "g"
            except Exception:
                print "<io.py> matplotlib was unable to write histogram."
                raise
        print "1"
        # generate fiber visualization
        #try:
            # pd3_li = filter.maskFibers(pd3_out,abs(li)>max(abs(li))*.01, li)
            #ren = render.render(self.polydata)
            #ren.save_views(dirname)
        #except Exception:
        #    print "<io.py> vtk or rendering issue. Failed to save views."
        #    print "<io.py> polydata was saved to disk so you can re-render."
        #    raise

        #print "IMPLEMENT SAVING OF PARAMETERS TOO"

    def read(self, dirname, readpd=False, readdist=False):
        """Read output (class laterality.LateralityResults) for one subject."""

        if not os.path.isdir(dirname):
            print "<io.py> error: directory does not exist.", dirname

        if readpd:
            # input polydata
            reader = vtk.vtkPolyDataReader()
            fname = os.path.join(dirname, 'tractography_with_LI.vtk')
            reader.SetFileName(fname)
            reader.Update()
            self.polydata = reader.GetOutput()

        # input LI and other data values using pickle
        fname = os.path.join(dirname, 'pickle_laterality_index.txt')
        fid = open(fname, 'r')
        self.laterality_index = pickle.load(fid)
        fid.close()

        fname = os.path.join(dirname, 'pickle_left_hem_similarity.txt')
        fid = open(fname, 'r')
        self.left_hem_similarity = pickle.load(fid)
        fid.close()
        fname = os.path.join(dirname, 'pickle_right_hem_similarity.txt')
        fid = open(fname, 'r')
        self.right_hem_similarity = pickle.load(fid)
        fid.close()
        fname = os.path.join(dirname, 'pickle_hemisphere.txt')
        fid = open(fname, 'r')
        self.hemisphere = pickle.load(fid)
        fid.close()

        if readdist:
            fid = open(os.path.join(dirname, 'pickle_right_hem_distance.txt'), 'r')
            self.right_hem_distance = pickle.load(fid)
            fid.close()
            fid = open(os.path.join(dirname, 'pickle_left_hem_distance.txt'), 'r')
            self.left_hem_distance = pickle.load(fid)
            fid.close()

        self.directory = dirname

        #print "IMPLEMENT READING OF PARAMETERS TOO"
