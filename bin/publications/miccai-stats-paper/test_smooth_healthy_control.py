import glob
import os

import whitematteranalysis as wma

import multiprocessing

parallel_jobs = multiprocessing.cpu_count()
print 'CPUs detected:', parallel_jobs
parallel_jobs = 10
print 'Using N jobs:', parallel_jobs

indir = '/Users/odonnell/Data/CONTROLS_50_PNL/Tracts/tract_data/good'
outdir = '/Users/odonnell/Dropbox/Coding/OUTPUTS/MICCAI2012/test_smooth_healthy_control'

input_mask = "{0}/*.vtp".format(indir)
input_poly_datas = glob.glob(input_mask)

input_poly_datas = input_poly_datas[0:1]
print input_poly_datas

number_of_fibers = 4000
points_per_fiber = 30
#fiber_length = 35
fiber_length = 50
#kernel_sigmas = [2.0, 5.0, 10.0, 20.0, 30.0]
kernel_sigmas = [0.0000001, 5.0, 10.0, 15.0, 20.0]

for fname in input_poly_datas:

    pd = wma.io.read_polydata(fname)
    print 'processing'
    pd = wma.filter.preprocess(pd,fiber_length)
    print 'random downsampling'
    pd = wma.filter.downsample(pd, number_of_fibers)
    print 'smoothing'

    density_list = list()
    pd_list = list()
    
    for sigma in kernel_sigmas:
        pd2, weights = wma.filter.smooth(pd, fiber_distance_sigma = sigma, points_per_fiber=points_per_fiber, n_jobs=parallel_jobs, upper_thresh=40)

        density_list.append(weights)
        pd_list.append(pd2)
        # make an output subdirectory for views and data
        outdir_current =  os.path.join(outdir, 'smooth_{0:03d}mm'.format(int(sigma)))
        if not os.path.exists(outdir_current):
            os.makedirs(outdir_current)

        # write new one in current directory
        fname_initial = 'smooth_{0:03d}mm_'.format(int(sigma))
        fname_out = os.path.join(outdir_current, fname_initial + os.path.basename(fname))
        wma.io.write_polydata(pd2, fname_out)
        
        # make an output subdirectory for views and data
        outdir_current =  os.path.join(outdir, 'wgr2_smooth_{0:03d}mm_wgr2'.format(int(sigma)))
        if not os.path.exists(outdir_current):
            os.makedirs(outdir_current)

        # write one without junk fibers also
        fname_initial = 'smooth_{0:03d}mm'.format(int(sigma))
        fname_out = os.path.join(outdir_current, fname_initial + os.path.basename(fname))        
        pdA = wma.filter.mask(pd2, weights >= 2, weights)
        wma.io.write_polydata(pdA, fname_out)


    # save views with all the same colors.
    idx = 0
    for sigma in kernel_sigmas:
        pd = pd_list[idx]
        density = density_list[3]
        pd_to_render = wma.filter.mask(pd, density >0, density)
         
        # make an output subdirectory for views and data
        outdir_current =  os.path.join(outdir, 'smooth_{0:03d}mm'.format(int(sigma)))
        if not os.path.exists(outdir_current):
            os.makedirs(outdir_current)
        #ren = wma.render.render(pd_to_render, scalar_bar=True)
        #ren = wma.render.render(pd_to_render, scalar_range=[1,15])
        # len 40
        ren = wma.render.render(pd_to_render, scalar_range=[3,35])
        # len 50, n 4000
        ren = wma.render.render(pd_to_render, scalar_range=[20, 60])        
        ren.save_views(outdir_current)
        idx += 1
        

        
