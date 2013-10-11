import os
import glob
import matplotlib.pyplot as plt
import numpy
import scipy.stats

import vtk

import whitematteranalysis as wma

import multiprocessing


def add_array_to_polydata(pd, array, array_name='Test', array_type='Cell'):
    out_array = vtk.vtkFloatArray()
    for idx in range(len(array)):
        out_array.InsertNextTuple1(array[idx])
    out_array.SetName(array_name)
    ret = pd.GetCellData().AddArray(out_array)
    print ret
    pd.GetCellData().SetActiveScalars(array_name)
    return(pd)


number_of_fibers_per_subject = 3000
points_per_fiber = 30
fiber_length = 30

indir = '/Users/odonnell/Dropbox/Coding/OUTPUTS/MICCAI2012/controls_with_scalars'

input_mask = "{0}/*.vtk".format(indir)
input_poly_datas = glob.glob(input_mask)

input_poly_datas = input_poly_datas[0:1]
print input_poly_datas

# now visualize as-is, and after conversion to fixed length n=30 and back


# read in ones with scalars already
input_pds = list()
for fname in input_poly_datas:
    print fname
    pd = wma.io.read_polydata(fname)
    input_pds.append(pd)

# downsample pds
input_pds_downsampled = list()
downsample_indices = list()
for pd in input_pds:
    pd_input, fiber_indices1 = wma.filter.preprocess(pd, fiber_length, return_indices=True)
    pd_input, fiber_indices2 = wma.filter.downsample(pd_input, number_of_fibers_per_subject,return_indices=True)
    fiber_indices = fiber_indices1[fiber_indices2]
    #pd_input, fiber_indices = wma.filter.downsample(pd, number_of_fibers_per_subject,return_indices=True)
    downsample_indices.append(fiber_indices)
    input_pds_downsampled.append(pd_input)

# convert to fixed length and back

# convert to array representation
print 'Converting fibers to array representation'
fiber_array = wma.fibers.FiberArray()
fiber_array.convert_from_polydata(pd_input, points_per_fiber)
print 'Done converting fibers to array representation'
pd_fixed = fiber_array.convert_to_polydata()
print 'Done converting fibers back to pd'

# view both with FA (?)
# grab scalars of interest
input_mean_perp_per_subject = list()
input_mean_para_per_subject = list()
input_mean_md_per_subject = list()
input_mean_fa_per_subject = list()
pidx = 0
for pd in input_pds:
    print pidx, '----------------------------------------------------------'
    mean_perp = pd.GetCellData().GetArray('mean_perpendicular_diffusivity')
    mean_para = pd.GetCellData().GetArray('mean_parallel_diffusivity')
    mean_md = pd.GetCellData().GetArray('mean_MD')
    mean_fa = pd.GetCellData().GetArray('mean_FA')
    mean_perp_subj = list()
    mean_para_subj = list()
    mean_md_subj = list()
    mean_fa_subj = list()
    fiber_indices = downsample_indices[pidx]
    for idx in fiber_indices:
        mean_perp_subj.append(mean_perp.GetTuple1(idx))
        mean_para_subj.append(mean_para.GetTuple1(idx))
        mean_md_subj.append(mean_md.GetTuple1(idx))
        mean_fa_subj.append(mean_fa.GetTuple1(idx))
    input_mean_perp_per_subject.append(mean_perp_subj)    
    input_mean_para_per_subject.append(mean_para_subj)    
    input_mean_md_per_subject.append(mean_md_subj)    
    input_mean_fa_per_subject.append(mean_fa_subj)    
    pidx += 1

color_by = input_mean_fa_per_subject[0]
color_name = 'MeanFA'

pd_input = add_array_to_polydata(pd_input, color_by, array_name=color_name)
pd_fixed = add_array_to_polydata(pd_fixed, color_by, array_name=color_name)

ren_input = wma.render.render(pd_input, scalar_bar=True, scalar_range=[0.35,0.65])
ren_fixed = wma.render.render(pd_fixed, scalar_bar=True, scalar_range=[0.35,0.65])


