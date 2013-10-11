import whitematteranalysis as wma
import glob
import os

parallel_jobs = 2
number_of_fibers = 10000
fiber_length  = 20


#sigma = 10
input_directory = '/Users/odonnell/Dropbox/Data/TBI_FE_PNL/controls'
input_mask = "{0}/*.vtp".format(input_directory)
input_polydatas = glob.glob(input_mask)

#input_polydatas = input_polydatas[0:1]

wma.registration_functions.smooth_polydatas_from_disk(input_polydatas, sigma, number_of_fibers, fiber_length, parallel_jobs)
