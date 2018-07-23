#!/usr/bin/python

import sys
import os
import argparse
import time
import numpy
import datetime

models=["HISEQ_2500_RAPID_DUALBC"]

HISEQ_2500_RAPID_DUALBC=[[1,1500],[24,1],[5,72],[-1,380],[1,2880],[4,1],[-1,380],[1,5100],[4,1],[-1,380],[1,7920],[9,1],[-1,380]]
HISEQ_2500_RAPID_DUALBC_STRUCTURE=["R","B","B","R"]

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="Input Basecalls folder to duplicate over time.")
    parser.add_argument("-o", "--output", help="Output directory for BaseCalls folder."
                        " (default: current directory)", default=os.getcwd())
    #parser.add_argument("-d", "--delay", help="Delay (s) before writing a new cycle.", 
    #                    default=120)
    parser.add_argument("-r", "--reads", help="Read structure (e.g., 100R,8B,100R). Lowercase characters , e.g. '8b', can be used to simulate the time of a segment without copying any files (e.g. if barcode runtime should be simulated without having barcodes in the basecall files). The given structure must fulfill the requirements of the selected mode.")
    parser.add_argument("-m", "--mode", help="Sequencing machine and mode.\n\nAvailable modes: {0}".format(models), default="HISEQ_2500_RAPID_DUALBC")
    parser.add_argument("-s", "--speedup", help="Speed up the sequencing model by the given factor. (default: 1)", default=1)
    parser.add_argument("--getruntime", help="Instead of executing the simulation, show the total runtime for the given parameters.", action='store_true')
    args = parser.parse_args()
    
    return args

def parse_reads(reads):
    read_vector = reads.split(",")
    ### return vector contains elements like [100,'R'] or [8,'B']
    return_vector = [[i[:-1],i[-1:]] for i in read_vector]
    return return_vector

def check_read_structure(read_vector, mode_structure):
    if len(read_vector) != len(mode_structure):
        print "ERROR: The number of read segments is different for the given read structure ({0}) and the selected model ({1}). Exit.".format(len(read_vector), len(mode_structure))
        exit(1);
    for segment in range(0,len(read_vector)-1):
        if read_vector[segment][1].upper() != mode_structure[segment].upper():
            print "ERROR: The given read structure does not match the selected mode (Segment {0} is {1} but should be {2}). Exit.".format(segment+1, read_vector[segment][1].upper(), mode_structure[segment].upper())
            exit(1)

def select_mode(mode, read_vector):

    ### HISEQ 2500 in rapid mode and dual barcode
    if mode=="HISEQ_2500_RAPID_DUALBC":
        check_read_structure(read_vector, HISEQ_2500_RAPID_DUALBC_STRUCTURE)
        return HISEQ_2500_RAPID_DUALBC;

    ### Didn't select a valid mode
    else:
        print "ERROR: The selected mode is not valid. Available modes are described in the help of this tool. Exit."
        exit(1)

def get_delays(mode, read_vector):
    delays = []
    copy = []
    read_vec_pos = 0
    remaining_cycles = int(read_vector[read_vec_pos][0])
    selected_mode = select_mode(mode, read_vector)
    for subarray in selected_mode:
        write = read_vector[read_vec_pos][1].isupper()
        if subarray[0] < 0:
            delays.extend(numpy.repeat(subarray[1], remaining_cycles))
            copy.extend(numpy.repeat(write, remaining_cycles))
            read_vec_pos += 1
            if read_vec_pos >= len(read_vector):
                break
            remaining_cycles = read_vector[read_vec_pos][0]
        else:
            delays.extend(numpy.repeat(subarray[1], subarray[0]))
            copy.extend(numpy.repeat(write, subarray[0]))
            remaining_cycles = int(remaining_cycles) - int(subarray[0])
    return [delays,copy]

def print_runtime(delays):
    total_seconds = int(0)
    for seconds in delays[0]:
        total_seconds += int(seconds)
    hours = total_seconds // 60 // 60
    minutes = (total_seconds - (hours*60*60)) // 60
    seconds = total_seconds % 60
    sys.stdout.write("Total runtime for the selected run: {0}h {1}m {2}s\n".format(hours, minutes, seconds))
    sys.stdout.flush()

def main(argv):
    args = parse_args()

    reads = parse_reads(args.reads)
    delays = get_delays(args.mode,reads)
    root = os.path.abspath(args.input)

    if args.getruntime:
        print_runtime(delays)
        exit(0)

    output_dir = os.path.abspath(args.output) + "/BaseCalls"
    os.mkdir(output_dir)
    
    # descent into the different lane folders
    cycle_dict = {} # one list of paths to be written for every cycle
    for lane_dir in os.listdir(root):
        if lane_dir[0] != "L":
            continue

        lane_root = os.path.join(root, lane_dir)
        lane_files = os.listdir(lane_root)

        # create the lane folder in the BaseCalls output directory
        lane_copy = os.path.join(output_dir, lane_dir)
        os.mkdir(lane_copy)
        
        # get cycle folders and filter file names
        cycle_dirs = [f for f in lane_files 
                      if os.path.isdir(os.path.join(lane_root, f))]
        filter_files = [f for f in lane_files 
                        if os.path.splitext(f)[1] == ".filter"]

        # sort cyle dirs (by cyle number)..
        cycle_dirs = sorted(cycle_dirs, key = lambda s: int(s.split('.')[0][1:]))
        for i, cycle in enumerate(cycle_dirs):
            if (i not in cycle_dict):
                cycle_dict[i] = []

            cycle_dict[i].append((os.path.join(lane_root, cycle), os.path.join(lane_copy, cycle)))

        for f in filter_files:
            os.symlink(os.path.join(lane_root, f), os.path.join(lane_copy, f))

    copy_cycle=0
    sys.stdout.write('[{0}] Start basecalling\n'.format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
    sys.stdout.flush()
    for i in range(0,len(delays[0])):
        time.sleep(float(delays[0][i])/float(args.speedup))
        if delays[1][i]:
            for src,dst in cycle_dict[copy_cycle]:
                os.symlink(src,dst)
            sys.stdout.write('[{0}] Simulated cycle {1}\n'.format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),copy_cycle+1))
            sys.stdout.flush()
            copy_cycle += 1

if __name__=="__main__":
    main(sys.argv)

