#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/ipython

print "LAUREN should mean rotation be 0?"

import os
import glob

import whitematteranalysis as wma
import vtk
import numpy
print 'Read and preprocess'

minimum_length = 60

#indir = '/Users/odonnell/Desktop/OUTPUT/laterality_test_mp_allfibers'
indir = '/Users/odonnell/Data/TBI/TenMatchedTracts'
outdir = '/Users/odonnell/Dropbox/Coding/OUTPUTS/TBI/registration'

#number_of_fibers = 1000
#number_of_fibers = 150
#number_of_fibers_step_one = 30
#number_of_fibers_step_two = 40
#number_of_fibers_step_three = 60

number_of_fibers = 300
number_of_fibers_step_one = 50
number_of_fibers_step_two = 75
number_of_fibers_step_three = 100

import multiprocessing
number_of_jobs = multiprocessing.cpu_count()
print 'CPUs detected:', number_of_jobs
#number_of_jobs = 12
#number_of_jobs = 16

inputMask = "{0}/*.vtp".format(indir)
inputPolyDatas = glob.glob(inputMask)
input_pds = list()

print inputPolyDatas

for fname in inputPolyDatas:
    print fname
    pd = wma.io.read_polydata(fname)
    pd2 = wma.filter.preprocess(pd, minimum_length)
    pd3 = wma.filter.downsample(pd2, number_of_fibers)
    input_pds.append(pd3)

def view_polydatas(polydata_list):
    appender = vtk.vtkAppendPolyData()
    idx = 0
    for pd in polydata_list:
        nf = pd.GetNumberOfLines()
        print idx
        print nf
        mask = numpy.ones(nf)
        colors = numpy.multiply(mask, idx-1)
        pd2 = wma.filter.mask(pd, mask, colors)
        appender.AddInput(pd2)
        idx = idx + 1
    appender.Update()
    pd3 = appender.GetOutput()
    ren = wma.render.render(pd3)
    return ren

def transform_polydatas(input_pds, register):
    transforms = register.convert_transforms_to_vtk()
    idx = 0
    output_pds = list()
    for transform in transforms:
        transformer = vtk.vtkTransformPolyDataFilter()
        transformer.SetInput(input_pds[idx])
        transformer.SetTransform(transform)
        transformer.Update()
        pd = transformer.GetOutput()
        output_pds.append(pd)
        idx = idx + 1
    return output_pds


def transform_polydatas_from_disk(input_pd_fnames, register, outdir):
    transforms = register.convert_transforms_to_vtk()
    for idx in range(0, len(input_pd_fnames)):
        transformer = vtk.vtkTransformPolyDataFilter()
        fname = input_pd_fnames[idx]
        print fname
        pd = wma.io.read_polydata(fname)
        transformer.SetInput(pd)
        transformer.SetTransform(transforms[idx])
        transformer.Update()
        pd2 = transformer.GetOutput()
        writer = vtk.vtkPolyDataWriter()
        writer.SetInput(pd2)
        fname = 'white_matter_{:04}.vtk'.format(idx)
        writer.SetFileName(os.path.join(outdir, fname))
        #writer.SetFileName(os.path.join(outdir, str(idx) + '.vtk'))
        writer.SetFileTypeToBinary()
        writer.Write()
        del writer
        del transformer
        del pd2
        del pd
    
# view input data
#ren = view_polydatas(input_pds)
#outdir_subdir = os.path.join(outdir, 'input')
#os.mkdir(outdir_subdir)
#ren.save_views(outdir_subdir)
#del ren

register = wma.congeal.CongealTractography()
register.parallel_jobs = number_of_jobs
#register.parallel_jobs = 10
#register.threshold = 5
register.threshold = 0
register.sigma = 10
register.points_per_fiber = 5
register.fiber_sample_size = number_of_fibers_step_one

# inputs are fixed, moving
for pd in input_pds:
    register.add_subject(pd)

model_pds=list()
for subj in register._subjects:
    model_pds.append(subj._original_fibers.convert_to_polydata())

#ren = view_polydatas(model_pds)

# RUN registration (initial)
register.translate_only()
register.maxfun = 150
register.compute()

# RUN registration (initial)
register.rotate_only()
register.compute()

transform_polydatas_from_disk(inputPolyDatas, register, outdir)


# view output data
#output_pds = transform_polydatas(input_pds, register)
#ren = view_polydatas(output_pds)
#outdir_subdir = os.path.join(outdir, 'output')
#os.mkdir(outdir_subdir)
#ren.save_views(outdir_subdir)
#del ren


# narrower rho parameters, solution should move less
inc_rot = (3.0 / 180.0) * numpy.pi
inc_trans = 5.0
inc_scale = .2 
register.set_rhobeg(inc_rot, inc_trans, inc_scale)
inc_rot = (1.0 / 180.0) * numpy.pi
inc_trans = 2.0
inc_scale = .01
register.set_rhoend(inc_rot, inc_trans, inc_scale)

# RUN registration (second)
register.translate_only()
register.maxfun = 250
register.compute()

register.rotate_only()
register.compute()


# test more...
inc_rot = (1 / 180.0) * numpy.pi
inc_trans = 1.0
inc_scale = .05 
register.set_rhobeg(inc_rot, inc_trans, inc_scale)
inc_rot = (.5 / 180.0) * numpy.pi
inc_trans = .5
inc_scale = .01
register.set_rhoend(inc_rot, inc_trans, inc_scale)

register.translate_and_rotate_and_scale()
#register.maxfun = 200
register.fiber_sample_size = number_of_fibers_step_three
register.maxfun = 1000
register.compute()
