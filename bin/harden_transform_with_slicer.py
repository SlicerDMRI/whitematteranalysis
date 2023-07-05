#!/usr/bin/env python
import argparse
import slicer
import os
import glob
import shutil

def harden_transform(polydata, transform, inverse, outdir):
    
    polydata_base_path, polydata_name = os.path.split(polydata)
    output_name = os.path.join(outdir, polydata_name)
    
    if os.path.exists(output_name):
        return
    
    check_load, polydata_node = slicer.util.loadFiberBundle(str(polydata), 1)
    if not check_load:
        print('Could not load polydata file:', polydata)
        return

    check_load, transform_node = slicer.util.loadTransform(str(transform), 1)
    if not check_load:
        print('Could not load transform file:', transform)
        return

    if polydata_node.GetPolyData().GetNumberOfCells() == 0:
        print('Empty cluster:', polydata)
        shutil.copyfile(polydata, output_name)
        return

    if inverse == "1":
        transform_node.Inverse()
    
    logic = slicer.vtkSlicerTransformLogic()
    t_node_id = transform_node.GetID()

    # harden transform
    polydata_node.SetAndObserveTransformNodeID(t_node_id)
    logic.hardenTransform(polydata_node)

    slicer.util.saveNode(polydata_node, output_name)

def main():
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
    
    if os.path.isfile(args.polydata):
        harden_transform(args.polydata, args.transform, args.inverse, args.outdir)
    elif os.path.isdir(args.polydata):
        input_mask = f"{args.polydata}/*.vtk"
        input_mask2 = f"{args.polydata}/*.vtp"
        input_pd_fnames = glob.glob(input_mask) + glob.glob(input_mask2)
        input_polydatas = sorted(input_pd_fnames)
    
        for polydata in input_polydatas:
            harden_transform(polydata, args.transform, args.inverse, args.outdir)

if __name__ == '__main__':
    main()
