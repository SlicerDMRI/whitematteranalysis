#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import glob
import multiprocessing
import os

import numpy as np
import vtk
from joblib import Parallel, delayed

import whitematteranalysis as wma


def _build_arg_parser():

    parser = argparse.ArgumentParser(
        description="Make sign change of the gradient direction",
        epilog="Written by Fan Zhang, fzhang@bwh.harvard.edu")
    parser.add_argument(
        'inputnrrd',
        help='Nrrd header in .nhdr format')
    parser.add_argument(
        'outputnrrd',
        help='New Nrrd header')
    parser.add_argument(
        '-d', action="store", dest="dim", type=str,
        help='The dimension to change: x, y, or z.')

    return parser


def _parse_args(parser):

    return parser.parse_args()


def main():

    parser = _build_arg_parser()
    args = _parse_args(parser)

    inputnrrd = args.inputnrrd
    outputnrrd = args.outputnrrd

    with open(inputnrrd) as f:
        contents = f.readlines()

    outstr = ''
    for line in contents:
        tmpstr = line
        if line.startswith('DWMRI_gradient'):
            strs = line.split(':=')
            dirs = strs[1].split(' ')

            if args.dim == 'x':
                val = float(dirs[0])
                if val != 0:
                    val = -val
                    dirs[0] = str(val)
            elif args.dim == 'y':
                val = float(dirs[1])
                if val != 0:
                    val = -val
                    dirs[1] = str(val)
            elif args.dim == 'z':
                val = float(dirs[2])
                if val != 0:
                    val = -val
                    dirs[2] = str(val) + '\n'

            tmpstr = strs[0] + ':=' + dirs[0] + ' ' + dirs[1] + ' ' + dirs[2]

        outstr += tmpstr

    with open(outputnrrd, 'w') as o:
        o.write(outstr)

    print(f'Done! New nrrd header {outputnrrd}')


if __name__ == "__main__":
    main()
