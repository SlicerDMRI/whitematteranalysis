import argparse
import pstats



parser = argparse.ArgumentParser(
    description="Print info from a cProfile text file to the terminal.",
    epilog="Written by Lauren O\'Donnell, odonnell@bwh.harvard.edu.",
    version='1.0')
parser.add_argument(
    'input_file',
    help='Text file containing cProfile output.')
parser.add_argument(
    '-n', action="store", dest="num_lines", type=int, default=10,  
    help='Lines of profile output to display (number of functions in list).')
parser.add_argument(
    '-o', action="store", dest="order", type=str, default='time',  
    help='Sort by time (default) or number of calls (calls) or cumulative time (cumulative)')
parser.add_argument(
    '-f', action="store", dest="input_file_2", type=str,
    help='File to compare times')
args = parser.parse_args()


p = pstats.Stats(args.input_file)
p.strip_dirs().sort_stats(args.order).print_stats(args.num_lines)

if args.input_file_2 is not None:
    p2 = pstats.Stats(args.input_file_2)
    p2.strip_dirs().sort_stats(args.order).print_stats(args.num_lines)






