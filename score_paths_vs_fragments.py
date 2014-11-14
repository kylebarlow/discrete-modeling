#!/usr/bin/python

program_description = "Script to compare paths to fragments"

import numpy as np
import os
import sys
import argparse

class FragmentData:
    def __init__ (self, split_lines):
        self.resnum_list = []
        self.aa_list = []
        self.secondary_struct_list = []
        self.phi_list = []
        self.psi_list = []
        self.omega_list = []
        self.x_list = []
        self.y_list = []
        self.z_list = []
        for i, data in enumerate(split_lines):
            if i == 0:
                self.pdb_id = data[0]
                self.chain = data[1]
            else:
                assert( self.pdb_id == data[0] )
                assert( self.chain == data[1] )

            self.resnum_list.append( int(data[2]) )
            assert( len(data[3]) == 1 )
            self.aa_list.append( data[3] )
            assert( len(data[4]) == 1 )
            self.secondary_struct_list.append( data[4] )
            self.phi_list.append( float(data[5]) )
            self.psi_list.append( float(data[6]) )
            self.omega_list.append( float(data[7]) )
            self.x_list.append( float(data[8]) )
            self.y_list.append( float(data[9]) )
            self.z_list.append( float(data[10]) )

        self.length = len(split_lines)

    def get_coords(self):
        return [(x, y, z) for x, y, z in zip(self.x_list, self.y_list, self.z_list)]

    def __len__ (self):
        return self.length

def parse_fragment_data(frag_file):
    fragment_data = []
    with open(frag_file, 'r') as f:
        split_lines = []
        for line in f:
            split_line = line.strip().split()
            if len(split_line) >= 11:
                split_lines.append(split_line)
            elif len(split_line) == 0 and len(split_lines) > 0:
                fragment_data.append( FragmentData(split_lines) )
                split_lines = []

    return fragment_data

def parse_path_data( path_file ):
    path_data = []
    with open(path_file, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                path_data.append( [int(x) for x in line.strip().split()] )
    return path_data

def parse_anchor_data( anchor_file ):
    anchor_points = {}
    with open(anchor_file, 'r') as f:
        for line in f:
            split_line = line.strip().split('|')
            if len( split_line ) == 6:
                anchor_points[int(split_line[1])] = (float(split_line[2]), float(split_line[3]), float(split_line[4]))
    return anchor_points

def main():
    parser = argparse.ArgumentParser(description=program_description)

    parser.add_argument('-f', '--fragment_file',
                        required = True,
                        help = 'Fragment file to read from')
    parser.add_argument('-p', '--path_file',
                        required = True,
                        help = 'Path file to read from')
    parser.add_argument('-a', '--anchor_file',
                        required = True,
                        help = 'Anchor file to read from')

    args = parser.parse_args()

    assert( os.path.isfile( args.fragment_file ) )
    assert( os.path.isfile( args.path_file) )
    assert( os.path.isfile( args.anchor_file) )

    fragment_data = parse_fragment_data( args.fragment_file )
    all_coords = [x.get_coords() for x in fragment_data]
    path_data = parse_path_data( args.path_file )
    anchor_points = parse_anchor_data( args.anchor_file )
    for path_coords in [[anchor_points[x] for x in path] for path in path_data]:
        print path_coords

if __name__ == "__main__":
    main()
