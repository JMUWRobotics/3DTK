import math
import random
import os
import subprocess
import sys
import argparse


def main():
    # Evaluate arguments

    parser = argparse.ArgumentParser(prog = 'accuracy_evaluation', description='Compares accuracy of two ICP-Algorithms')
    parser.add_argument("-d", "--directory", required = True, help="output directory")
    parser.add_argument("-s", "--start", required = False, default = 1, type = int,  help = "start at scan number")
    parser.add_argument("-e", "--end", required = False, type = int, help = "end at scan number")

    args = parser.parse_args()
    dir = args.directory
    scan_nr = args.start
    end = args.end


    if not os.path.exists(dir):
        parser.error("Path does not exist")

    # Open scans
    while end is None or scan_nr <= end:
        fp_file = dir + f"/fix_scan{scan_nr:03d}.frames"
        ref_file = dir +f"/scan{scan_nr:03d}.frames"
        try:
            with open(fp_file, "r") as fp:
                fp_frame = fp.readlines()
            with open(ref_file, "r") as ref:
                ref_frame = ref.readlines()
        except FileNotFoundError:
            print(f"One or both frames for scan{scan_nr:03d} missing \n")
            break

        evaluate(fp_frame, ref_frame, scan_nr, dir)
        subprocess.run(["gnuplot", "-c", "Plot_translation_error.plt", str(scan_nr), dir])

        scan_nr += 1



def evaluate(fp_frame, ref_frame, scan_nr, dir):
    fp_t_err = get_translation_errors(fp_frame)
    ref_t_err = get_translation_errors(ref_frame)

    #fp_r_err = get_rotation_errors(fp_frame)
    #ref_r_err = get_rotation_errors(ref_frame)

    # write transformation errors
    with open(f"{dir}/fix_scan{scan_nr:03d}_trans.dat", "w") as fp:
        for i, error in enumerate(fp_t_err):
            fp.write(f"{i}\t{error}\n")

    with open(f"{dir}/ref_scan{scan_nr:03d}_trans.dat", "w") as ref:
        for i, error in enumerate(ref_t_err):
            ref.write(f"{i}\t{error}\n")

    # write rotation errors
    #with open(f"{dir}/fix_scan{scan_nr:03d}_rot.dat", "w") as fp:
    #    for i, error in enumerate(fp_r_err):
    #        fp.write(f"{i}\t{error}\n")

    #with open(f"{dir}/ref_scan{scan_nr:03d}_rot.dat", "w") as ref:
    #    for i, error in enumerate(ref_r_err):
    #        ref.write(f"{i}\t{error}\n")

# def evaluate_median_error

def get_translation_errors(frame):
    errors = []
    for line in frame:
        values = line.split()
        if (values[-1] == "1"):
            values = [float(value) for value in values]
            t = [values[12], values[13], values[14]]
            errors.append(math.sqrt(t[0]*t[0] + t[1]*t[1] + t[2]*t[2]))

    return errors

def get_rotation_errors(frame):
    errors = []
    for line in frame:
        values = line.split()
        if (values[-1] == "1"):
            values = [float(value) for value in values]
            R = [[values[0], values[1], values[2]], [values[4], values[5], values[6]], [values[8], values[9], values[10]]]

    return errors


if __name__ == "__main__":
    main()
