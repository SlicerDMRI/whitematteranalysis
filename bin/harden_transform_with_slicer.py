#!/usr/bin/env python
import argparse
import slicer
import os
import glob

parser = argparse.ArgumentParser(
    description="Harden transform with Slicer.",
    epilog="Written by Fan Zhang, fzhang@bwh.harvard.edu")

parser.add_argument("-v", "--version",
    action="version", default=argparse.SUPPRESS,
    version='1.0',
    help="Show program's version number and exit")
parser.add_argument(
    'polydata',
    help='')
parser.add_argument(
    'transform',
    help='')
parser.add_argument(
    'inverse',
    help='')
parser.add_argument(
    'outdir',
    help='')

args = parser.parse_args()

def harden_transform(polydata, transform, inverse, outdir):

    check_load, polydata_node = slicer.util.loadModel(str(polydata), 1)
    if not check_load:
        print 'Could not load polydata file:', polydata
        return

    check_load, transform_node = slicer.util.loadTransform(str(transform), 1)
    if not check_load:
        print 'Could not load transform file:', transform
        return

    if inverse == "1":
        transform_node.Inverse()

    logic = slicer.vtkSlicerTransformLogic()
    t_node_id = transform_node.GetID()

    # harden transform
    polydata_node.SetAndObserveTransformNodeID(t_node_id)
    logic.hardenTransform(polydata_node)

    polydata_base_path, polydata_name = os.path.split(polydata)
    output_name = polydata_name.replace(".", "_trans.")
    if inverse == "1":
        output_name = polydata_name.replace(".", "_inv_trans.")
    slicer.util.saveNode(polydata_node, os.path.join(outdir, output_name))

if os.path.isfile(args.polydata):
    harden_transform(args.polydata, args.transform, args.inverse, args.outdir)
elif os.path.isdir(args.polydata):
    input_mask = "{0}/*.vtk".format(args.polydata)
    input_mask2 = "{0}/*.vtp".format(args.polydata)
    input_pd_fnames = glob.glob(input_mask) + glob.glob(input_mask2)
    input_polydatas = sorted(input_pd_fnames)

    for polydata in input_polydatas:
        harden_transform(polydata, args.transform, args.inverse, args.outdir)