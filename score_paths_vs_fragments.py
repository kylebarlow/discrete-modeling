#!/usr/bin/python

program_description = "Script to compare paths to fragments"

import numpy as np
import os
import sys
import argparse
import gzip
import io
import time
import math

try:
    import pyRMSD
    from pyRMSD.matrixHandler import MatrixHandler
    import pyRMSD.RMSDCalculator
except ImportError:
    pyRMSD = None
    import Bio.PDB
    import Bio.PDB.Atom as Atom

class Reporter:
    def __init__(self,task,entries='files',print_output=True):
        self.print_output=print_output
        self.start=time.time()
        self.entries=entries
        self.lastreport=self.start
        self.task=task
        self.report_interval=1 # Interval to print progress (seconds)
        self.n=0
        self.completion_time = None
        if self.print_output:
            print '\nStarting '+task
        self.total_count = None # Total tasks to be processed
    def set_total_count(self, x):
        self.total_count = x
    def report(self,n):
        self.n=n
        t=time.time()
        if self.print_output and self.lastreport<(t-self.report_interval):
            self.lastreport=t
            if self.total_count:
                percent_done = float(self.n)/float(self.total_count)
                time_now = time.time()
                est_total_time = (time_now-self.start) * (1.0/percent_done)
                time_remaining = est_total_time - (time_now-self.start)
                minutes_remaining = math.floor(time_remaining/60.0)
                seconds_remaining = int(time_remaining-(60*minutes_remaining))
                sys.stdout.write("  Processed: "+str(n)+" "+self.entries+" (%.1f%%) %02d:%02d\r"%(percent_done*100.0, minutes_remaining, seconds_remaining) )
            else:
                sys.stdout.write("  Processed: "+str(n)+" "+self.entries+"\r" )
            sys.stdout.flush()
    def increment_report(self):
        self.report(self.n+1)
    def decrement_report(self):
        self.report(self.n-1)
    def add_to_report(self,x):
        self.report(self.n+x)
    def done(self):
        self.completion_time = time.time()
        if self.print_output:
            print 'Done %s, processed %d %s, took %.3f seconds\n' % (self.task,self.n,self.entries,self.completion_time-self.start)
    def elapsed_time(self):
        if self.completion_time:
            return self.completion_time - self.start
        else:
            return time.time() - self.start

class Fragment:
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
        self.coords_cached = False

    def get_coords(self):
        if not self.coords_cached:
            self.coords = [(x, y, z) for x, y, z in zip(self.x_list, self.y_list, self.z_list)]
            self.coords_cached = True
        return self.coords

    def __len__ (self):
        return self.length
    
    def __repr__(self):
        return str(self.get_coords())

class FragmentPositionData:
    def __init__ (self, start_position):
        self.start_position = start_position
        self.length = None
        self.fragment_list = []
        self.end_position = None

    def add_fragment(self, new_fragment, start_position):
        assert( self.start_position == start_position )

        if self.length:
            assert( self.length == len(new_fragment) )
        else:
            self.length = len(new_fragment)
            self.end_position = self.start_position + self.length
        self.fragment_list.append( new_fragment )

    def __iter__(self):
        return iter(self.fragment_list)

    def __getitem__(self, key):
        return self.fragment_list[key]

    def __len__(self):
        return len(self.fragment_list)

def parse_fragment_data(frag_file):
    fragment_data = {}

    if frag_file.endswith('.gz'):
        f = io.TextIOWrapper(io.BufferedReader(gzip.open(frag_file)))
    else:
        f = open(frag_file, 'r')

    split_lines = []
    current_start_position = None
    for line in f:
        split_line = line.strip().split()
        if len(split_line) >= 11:
            split_lines.append(split_line)
        elif len(split_line) == 4:
            if 'position' in split_line[0]:
                current_start_position = int(split_line[1])
        elif len(split_line) == 0 and len(split_lines) > 0:
            fd = Fragment(split_lines)
            assert( current_start_position )
            if current_start_position not in fragment_data:
                fragment_data[current_start_position] = FragmentPositionData(current_start_position)
            fragment_data[current_start_position].add_fragment( fd, current_start_position )
            split_lines = []

    f.close()

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

def calc_rms(ref_coords, alt_coords):
    assert( len(ref_coords) == len(alt_coords) )
    if pyRMSD:
        calculator = pyRMSD.RMSDCalculator.RMSDCalculator("QCP_SERIAL_CALCULATOR", np.array([ref_coords, alt_coords]))
        return calculator.pairwiseRMSDMatrix()[0]
    else:
        super_imposer = Bio.PDB.Superimposer()
        ref_atoms = [Atom.Atom('CA', coords, 0.0, 1.0, '', ' CA ', i+1, element='C') for i, coords in enumerate(ref_coords)]
        alt_atoms = [Atom.Atom('CA', coords, 0.0, 1.0, '', ' CA ', i+1, element='C') for i, coords in enumerate(alt_coords)]
        super_imposer.set_atoms(ref_atoms, alt_atoms)
        return super_imposer.rms

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
    # all_coords = [x.get_coords() for x in fragment_data.values()]
    path_data = parse_path_data( args.path_file )
    anchor_points = parse_anchor_data( args.anchor_file )
    r = Reporter('calculating RMS for all paths vs. all first fragments for each position', entries='paths')
    r.total_count = len(path_data[:50])
    print 'total count:', r.total_count
    for i, path_coords in enumerate([[anchor_points[x] for x in path] for path in path_data[:50]]):
        path_length = len(path_coords)
        best_rms = float("inf")
        rms_results = []
        for fragment_position in fragment_data.values():
            for j, fragment in enumerate(fragment_position):
                fragment_coords = fragment.get_coords()
                rms = calc_rms(fragment_coords, path_coords)
                rms_results.append( (rms, fragment_position.start_position, j) )
        rms_results.sort()
        r.increment_report()

    r.done()

if __name__ == "__main__":
    main()
