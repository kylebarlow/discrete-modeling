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
import gc
import traceback
import json

import ctypes
import multiprocessing
from multiprocessing.sharedctypes import Value, Array, RawArray

# Use pyRMSD for RMSD calculations if available as it's much faster than BioPython
try:
    import pyRMSD
    from pyRMSD.matrixHandler import MatrixHandler
    import pyRMSD.RMSDCalculator
    from pyRMSD.availableCalculators import availableCalculators
    # Use CUDA for GPU calculations, if avialable
    if 'QCP_CUDA_MEM_CALCULATOR' in availableCalculators():
        pyrmsd_calc = 'QCP_CUDA_MEM_CALCULATOR'
    else:
        pyrmsd_calc = 'QCP_SERIAL_CALCULATOR'
except ImportError:
    print 'Warning: could not import faster RMSD module. Falling back to slower biopython...'
    pyRMSD = None
    import Bio.PDB
    import Bio.PDB.Atom as Atom

# Use MDAnalysis for distance calculations if available
try:
    import MDAnalysis
    import MDAnalysis.analysis
    # from MDAnalysis.analysis import distances
except ImportError:
    print 'Warning: could not import faster distance module. Falling back...'
    MDAnalysis = None
    sys.exit(0)

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
        return [(x, y, z) for x, y, z in zip(self.x_list, self.y_list, self.z_list)]

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
    fragment_position_data = {}

    if frag_file.endswith('.gz'):
        f = io.TextIOWrapper(io.BufferedReader(gzip.open(frag_file)))
    else:
        f = open(frag_file, 'r')

    split_lines = []
    current_start_position = None
    coordinate_count = 0
    for line in f:
        split_line = line.strip().split()
        if len(split_line) >= 11:
            split_lines.append(split_line)
        elif len(split_line) == 4:
            if 'position' in split_line[0]:
                current_start_position = int(split_line[1])
        elif len(split_line) == 0 and len(split_lines) > 0:
            fd = Fragment(split_lines)
            coordinate_count += len(fd)
            assert( current_start_position )
            if current_start_position not in fragment_position_data:
                fragment_position_data[current_start_position] = FragmentPositionData(current_start_position)
            fragment_position_data[current_start_position].add_fragment( fd, current_start_position )
            split_lines = []

    f.close()

    # Even though there are these nice Python objects with all the data in it,
    # here they are converted to nasty C-style arrays to be shared among processes
    global coord_array
    global starting_position_array
    global fragment_number_array
    global reporter_n
    reporter_n = multiprocessing.Value('L', 0)
    coord_array = multiprocessing.sharedctypes.Array('d', coordinate_count*3, lock=False)
    starting_position_array = multiprocessing.sharedctypes.Array('L', coordinate_count, lock=False)
    fragment_number_array = multiprocessing.sharedctypes.Array('L', coordinate_count, lock=False)
    coord_array_i = 0
    i = 0
    for start_position in sorted(fragment_position_data.keys()):
        for fragment_index, fragment in enumerate(fragment_position_data[start_position]):
            for coords in fragment.get_coords():
                starting_position_array[i] = start_position
                fragment_number_array[i] = fragment_index
                i += 1
                for coord in coords:
                    coord_array[coord_array_i] = coord
                    coord_array_i += 1

    assert( coord_array_i == coordinate_count*3)
    assert( i == coordinate_count )

    return fragment_position_data

def parse_path_data( path_file ):
    path_data = []
    path_length = None
    with open(path_file, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                new_path = [int(x) for x in line.strip().split()] 
                new_path_length = len(new_path)
                if path_length:
                    assert( path_length == new_path_length )
                else:
                    path_length = new_path_length
                path_data.append( [int(x) for x in line.strip().split()] )
    return (path_data, path_length)

def parse_anchor_data( anchor_file ):
    anchor_points = {}
    with open(anchor_file, 'r') as f:
        for line in f:
            split_line = line.strip().split('|')
            if len( split_line ) == 6:
                anchor_points[int(split_line[1])] = (float(split_line[2]), float(split_line[3]), float(split_line[4]))
    return anchor_points

def calc_drms(ref_coords, alt_coords, ref_distances, alt_distances):
    assert( len(ref_coords) == len(alt_coords) )
    if MDAnalysis:
        MDAnalysis.analysis.distances.self_distance_array(ref_coords, result=ref_distances)
        MDAnalysis.analysis.distances.self_distance_array(alt_coords, result=alt_distances)
        return np.sum( np.absolute( ref_distances-alt_distances ) )
    else:
        total_drms = 0
        for atom1_ref, atom1_alt in zip(ref_coords, alt_coords):
            for atom2_ref, atom2_alt in zip(ref_coords, alt_coords):
                total_drms += abs(abs(np.linalg.norm(np.array(atom2_ref)-np.array(atom1_ref)))-abs(np.linalg.norm(np.array(atom2_alt)-np.array(atom1_alt))))
        return total_drms

def calc_rms(ref_coords, alt_coords):
    # print ref_coords, alt_coords
    assert( len(ref_coords) == len(alt_coords) )
    if pyRMSD:
        calculator = pyRMSD.RMSDCalculator.RMSDCalculator(pyrmsd_calc, np.array([ref_coords, alt_coords]))
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
    parser.add_argument('--single_thread',
                        default = False,
                        action = 'store_true',
                        help = 'Do not use multiprocessing')

    args = parser.parse_args()

    assert( os.path.isfile( args.fragment_file ) )
    assert( os.path.isfile( args.path_file) )
    assert( os.path.isfile( args.anchor_file) )

    print 'Loading fragment data'
    fragment_position_data = parse_fragment_data( args.fragment_file )
    print 'Loading path data'
    path_data, path_length = parse_path_data( args.path_file )
    print 'Loading anchor point data'
    anchor_points = parse_anchor_data( args.anchor_file )

    if args.single_thread:
        cpu_count = 1
    else:
        cpu_count = multiprocessing.cpu_count()
        pool = multiprocessing.Pool(cpu_count)
    results_dict = {}

    print 'Starting calculating RMS for all paths vs. all first fragments for each position'
    total_count = len(path_data)
    starting_time = time.time()

    def helper_callback(outer_results):
        for results_tuple in outer_results:
            path_number, rms_results = results_tuple
            results_dict[path_number] = rms_results

    path_nums_for_jobs = [[] for x in xrange(cpu_count)]
    path_coords_for_jobs = [[] for x in xrange(cpu_count)]
    for i, path_coords in enumerate([[anchor_points[x] for x in path] for path in path_data]):
        path_nums_for_jobs[ i%cpu_count ].append( i )
        path_coords_for_jobs[ i%cpu_count ].append( path_coords )

    for path_nums_list, path_coords_list in zip(path_nums_for_jobs, path_coords_for_jobs):
        if args.single_thread:
            helper_callback( paths_against_all_fragments(path_nums_list, path_coords_list, starting_time, total_count, path_length) )
        else:
            pool.apply_async(paths_against_all_fragments, (path_nums_list, path_coords_list, starting_time, total_count, path_length), callback=helper_callback)

    if not args.single_thread:
        pool.close()
        pool.join()

    n = int(reporter_n.value)
    completion_time = time.time()
    print 'Finished! Processed %d %s, took %.3f seconds\n' % (n, 'paths', completion_time-starting_time)

    print 'Dumping json'
    with open('results.json', 'w') as f:
        json.dump(results_dict, f)

def paths_against_all_fragments(path_nums_list, path_coords_list, starting_time, total_count, path_length):
    try:
        outer_results = []
        # Outer loop: loop through all paths assigned to this process
        for path_num, path_coords in zip(path_nums_list, path_coords_list):
            path_coords = [np.float32(x) for x in path_coords]
            rms_results = []

            fragment_coords = []
            last_fragment_number = fragment_number_array[0]
            last_starting_position = starting_position_array[0]
            coord_array_index = 0

            # Pre-allocated arrays for DRMS
            fragment_distances = np.zeros(path_length*(path_length-1)/2, dtype=np.float64)
            path_distances = np.zeros(path_length*(path_length-1)/2, dtype=np.float64)            

            # Inner loop: loop through all fragments to RMS against outer path
            for i in xrange(len(starting_position_array)):
                starting_position = starting_position_array[i]
                fragment_number = fragment_number_array[i]
                if (starting_position != last_starting_position) or (fragment_number != last_fragment_number):
                    # rms = calc_rms(fragment_coords, path_coords)

                    drms = calc_drms(np.array(fragment_coords), np.array(path_coords), fragment_distances, path_distances)

                    rms_results.append( (drms, last_starting_position, last_fragment_number) )
                    fragment_coords = []
                    last_fragment_number = fragment_number
                    last_starting_position = starting_position

                fragment_coords.append((
                    np.float32(coord_array[coord_array_index]),
                    np.float32(coord_array[coord_array_index+1]),
                    np.float32(coord_array[coord_array_index+2])
                ))
                coord_array_index += 3

            rms_results.sort()
            outer_results.append( (path_num, rms_results[:100]) )

            with reporter_n.get_lock():
                reporter_n.value += 1

                n = reporter_n.value
                t = time.time()
                percent_done = float(n)/float(total_count)
                time_now = time.time()
                est_total_time = (time_now-starting_time) * (1.0/percent_done)
                time_remaining = est_total_time - (time_now-starting_time)
                minutes_remaining = math.floor(time_remaining/60.0)
                seconds_remaining = int(time_remaining-(60*minutes_remaining))
                sys.stdout.write("  Processed: "+str(n)+" paths (%.1f%%) %02d:%02d\r"%(percent_done*100.0, minutes_remaining, seconds_remaining) )
                sys.stdout.flush()

        return outer_results

    except Exception as e:
        print('Caught exception in worker thread')

        # This prints the type, value, and stack trace of the
        # current exception being handled.
        traceback.print_exc()

        print()
        raise e

if __name__ == "__main__":
    main()
