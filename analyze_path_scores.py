#!/usr/bin/python

program_description = "Script to analyze output of score_paths_vs_fragments"

import numpy as np
import os
import sys
import argparse
import io
import gzip
import json
import score_paths_vs_fragments

def parse_results_file(filepath):
    if filepath.endswith('.gz'):
        f = io.TextIOWrapper(io.BufferedReader(gzip.open(filepath)))
    else:
        f = open(filepath, 'r')

    raw_data = json.load(f)
    f.close()

    data = {}
    for str_key in raw_data:
        data[int(str_key)] = raw_data[str_key]

    return data

def path_is_contiguous(path):
    last_x = path[0]
    for x in path[1:]:
        if x != last_x+1:
            return False
        last_x += 1
    return True

def main():
    parser = argparse.ArgumentParser(description=program_description)

    parser.add_argument('results_json_file',
                        help = 'Results file location')
    parser.add_argument('-p', '--path_file',
                        required = False,
                        default = None,
                        help = 'Path file to read from')

    args = parser.parse_args()

    assert( os.path.isfile( args.results_json_file) )
    assert( os.path.isfile( args.path_file) )

    results_data = parse_results_file(args.results_json_file)
    
    rmsds = []
    position_coverage_counts = {}
    for path_number in sorted(results_data.keys()):
        rmsds.append( [r[0] for r in results_data[path_number]] )
        for r in results_data[path_number]:
            starting_position_number = r[1]
            frag_length = 9 # Change this hard code
            for x in xrange(frag_length):
                position_number = starting_position_number + x
                if position_number not in position_coverage_counts:
                    position_coverage_counts[position_number] = 0
                position_coverage_counts[position_number] += 1

    mean_rmsds = np.mean(rmsds, axis=1)
    # print mean_rmsds
    # print position_coverage_counts

    if args.path_file:
        # Get easy paths that are contiguous
        paths, path_length = score_paths_vs_fragments.parse_path_data( args.path_file )
        contiguous_paths = []
        for i, path in enumerate(paths):
            if path_is_contiguous(path):
                contiguous_paths.append( (i, path) )

        for path_index, path in sorted(contiguous_paths):
            path_rank = 0
            for i, r in enumerate(results_data[path_index]):
                if r[1] == path[0]:
                    path_rank = i+1
                    break
                else:
                    path_rank = i
            print path[0], path_index, path_rank

if __name__ == "__main__":
    main()
