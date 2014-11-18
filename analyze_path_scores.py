#!/usr/bin/python

program_description = "Script to analyze output of score_paths_vs_fragments"

import numpy as np
import os
import sys
import argparse
import io
import gzip
import json

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

def main():
    parser = argparse.ArgumentParser(description=program_description)

    parser.add_argument('results_json_file',
                        help = 'Results file location')

    args = parser.parse_args()

    assert( os.path.isfile( args.results_json_file) )

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
    print mean_rmsds
    print position_coverage_counts

if __name__ == "__main__":
    main()
