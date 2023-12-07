#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os
import re

import numpy as np
import pandas


def _build_arg_parser():

    parser = argparse.ArgumentParser(
        description="Compute averaged statistics of the anatomical tracts. For diffusion measure with multiple statistics, only the Mean will be outputted.",
        epilog="Written by Fan Zhang, fzhang@bwh.harvard.edu.")
    parser.add_argument(
        'inputmeasure',
        help='A .csv or .txt measurement files.')
    parser.add_argument(
        'outoutmeasure',
        help='Output csv file.')
    parser.add_argument(
        '-appendedTractName', action="store", type=str, default="",
        help="Name of the appended tracts.")
    parser.add_argument(
        '-tractList', action="store", type=str, nargs='+',
        help='A list of tracts to be appended, e.g., AF_left AF_right.')

    return parser


def _parse_args(parser):

    return parser.parse_args()


def main():

    parser = _build_arg_parser()
    args = _parse_args(parser)

    stats = pandas.read_table(args.inputmeasure, delimiter=',')

    fields = []
    for col in stats.columns:
        fields.append(col)
    fields = np.array(fields)

    print(fields)

    print('All available tracts:')
    all_tracts = [f_name.replace('.Num_Points', '') for f_name in fields if f_name.endswith('.Num_Points')]
    print(f'N = {len(all_tracts)} :')
    print(all_tracts)

    append_list = []
    print('Tracts and related measures to be appended.')
    for t_idx, tract in enumerate(args.tractList):
        print(f'* {t_idx} - {tract}')
        indices = [index for index in range(len(fields)) if fields[index].startswith(tract+'.Num_') or re.search(tract+".*.Mean", fields[index])]
        if len(indices) == 0:
            print("Error: tract not founded in the input file.")
            exit()
        print(fields[indices])
        append_list.append(indices)

    append_list = np.array(append_list)
    append_measures = [field.replace(args.tractList[-1], '') for field in fields[indices]]

    print(f"Output measures: {append_measures}")

    avg_stats = [stats.to_numpy()[:, 0]]
    for m_idx, m_name in enumerate(append_measures):

        if m_name == '.Num_Points':

            avg_stat = np.sum(stats.to_numpy()[:, append_list[:, m_idx]], axis=1)

        elif m_name == '.Num_Fibers':

            avg_stat = np.sum(stats.to_numpy()[:, append_list[:, m_idx]], axis=1)

        elif m_name == '.Mean_Length': # weighted by NoS

            weight = stats.to_numpy()[:, append_list[:, 1]]
            val = stats.to_numpy()[:, append_list[:, m_idx]]

            val_weighted_sum = np.sum(val * weight, axis=1)
            weight_sum = np.sum(weight, axis=1)

            empty_indices = np.where(weight_sum == 0)[0]
            weight_sum[empty_indices] = 1

            avg_stat = val_weighted_sum / weight_sum
            avg_stat[empty_indices] = 0

        else: # weighted by NoP

            weight = stats.to_numpy()[:, append_list[:, 0]]
            val = stats.to_numpy()[:, append_list[:, m_idx]]

            val_weighted_sum = np.sum(val.astype(np.double) * weight.astype(np.double), axis=1)
            weight_sum = np.sum(weight, axis=1)

            empty_indices = np.where(weight_sum == 0)[0]
            weight_sum[empty_indices] = 1

            avg_stat = val_weighted_sum / weight_sum
            avg_stat[empty_indices] = -1

        avg_stats.append(avg_stat)

    avg_stats = np.array(avg_stats)
    avg_stats = avg_stats.transpose()

    column_names = ['subjectkey']
    for m in append_measures:
        column_names.append(args.appendedTractName+m)

    df = pandas.DataFrame(avg_stats, columns=column_names)

    print('Averaged tract measures:')
    print(df)

    df.to_csv(args.outoutmeasure, index=False)
    print()
    print(f'Output file at: {os.path.abspath(args.outoutmeasure)}')

    exit()


if __name__ == "__main__":
    main()
