import math
import random
import os
import sys
import argparse
import datetime

def cube(new_scan):
    for i in range(100):
        new_scan.write(str(random.uniform(-50,50)) + "\t" + str(random.uniform(-50,50)) + "\t" + "100" + "\n")
        new_scan.write(str(random.uniform(-50,50)) + "\t" + str(random.uniform(-50,50)) + "\t" + "150" + "\n")
        new_scan.write(str(random.uniform(-50,50)) + "\t" + "-50" + "\t" + str(random.uniform(50, 150)) + "\n")
        new_scan.write(str(random.uniform(-50,50)) + "\t" + "50" + "\t" + str(random.uniform(50, 150)) + "\n")
        new_scan.write("-50" + "\t" + str(random.uniform(-50,50)) + "\t" + str(random.uniform(50, 150)) + "\n")
        new_scan.write("50" + "\t" + str(random.uniform(-50,50)) + "\t" + str(random.uniform(50, 150)) + "\n")

def volume_cube(new_scan):
    for i in range(1200):
        new_scan.write(str(random.uniform(-50, 50)) + "\t" + str(random.uniform(-50, 50)) + "\t" + str(random.uniform(50, 150)) + "\n")

def sphere(new_scan): # https://mathworld.wolfram.com/SpherePointPicking.html
    # Radius
    r = 50
    for i in range(1200):
        u = random.uniform(0, 1)
        v = random.uniform(0, 1)

        phi = 2 * math.pi* u
        theta = math.acos(2*v-1)

        # Polar coordinates to xyz values
        new_scan.write(str(r * math.sin(theta) * math.cos(phi)) + "\t" + str(r * math.sin(theta) * math.sin(phi)) + "\t" + str(r * math.cos(theta)+125) + "\n")

def volume_sphere(new_scan): # https://codepen.io/thebabydino/pen/abLRNqV
    for i in range(1200):
        # Radius
        r = random.uniform(0, 1)
        r = math.pow(r, 1/3) * 50

        u = random.uniform(0, 1)
        v = random.uniform(0, 1)

        phi = 2 * math.pi* u
        theta = math.acos(2*v-1)

        new_scan.write(str(r * math.sin(theta) * math.cos(phi)) + "\t" + str(r * math.sin(theta) * math.sin(phi)) + "\t" + str(r * math.cos(theta)+125) + "\n")



def main():
    # Evaluate arguments

    parser = argparse.ArgumentParser(prog = 'generate', description='Generate point cloud in different shapes')
    parser.add_argument("-d", "--directory", required = True, help="output directory")
    parser.add_argument("-s",  "--shape", required = True, choices = ['cube', 'volume_cube', 'sphere', 'volume_sphere'],  help = "Shape. Possible shapes: \n \t cube \n \t volume_cube \n \t sphere \n \t volume_sphere" )
    parser.add_argument("-n", "--number", type = int, default = 1, help = "Amount of scans to be generated")
    parser.add_argument("-tmin", "--translation_minimum", type = float, nargs = 3, help = "Lower border for translation")
    parser.add_argument("-tmax", "--translation_maximum", type = float, nargs = 3, help = "Upper border for translation")
    parser.add_argument("-t", "--translation", type = float, nargs = 3, default = [10,10, 10], help = "Symmetric border for translation")
    parser.add_argument("-rmin", "--rotation_minimum", type = float, nargs = 3, help = "Lower border for rotation")
    parser.add_argument("-rmax", "--rotation_maximum", type = float, nargs = 3, help = "Upper border for rotation")
    parser.add_argument("-r", "--rotation", type = float, nargs = 3, default = [20, 20, 20], help = "Symmetric rotation border for rotation")


    args = parser.parse_args()
    dir = args.directory
    shape = args.shape
    n = args.number
    # min/max has higher priority than symmetric t / r
    tmin = args.translation_minimum if args.translation_minimum is not None else [-x for x in args.translation]
    tmax = args.translation_maximum if args.translation_maximum is not None else args.translation
    rmin = args.rotation_minimum if args.rotation_minimum is not None else [-x for x in args.rotation]
    rmax = args.rotation_maximum if args.rotation_maximum is not None else args.rotation




    seed = int(datetime.datetime.now().strftime("%Y%m%d%H%M%S"))
    random.seed(seed)
    #Create scan:
    if not os.path.exists(dir):
        parser.error("Path does not exist")

    for i in range(0, n):
        with open(dir + f"/scan{i:03d}.3d", "w", encoding="utf-8") as new_scan:
            match (shape):
                case "cube":
                    cube(new_scan)
                case "volume_cube":
                    volume_cube(new_scan)
                case "sphere":
                    sphere(new_scan)
                case "volume_sphere":
                    volume_sphere(new_scan)
                case _:
                    parser.error("Invalid shape argument")
        with open(dir + f"/scan{i:03d}.pose","w", encoding="utf-8" ) as new_pose:
            if(i == 0):
                translation = [0, 0, 0]
                rotation = [0, 0, 0]
            else:
                translation = [random.uniform(tmin[j], tmax[j]) for j in range(3)]
                rotation = [random.uniform(rmin[j], rmax[j]) for j in range(3)]

            new_pose.write(f"{translation[0]} {translation[1]} {translation[2]}\n")
            new_pose.write(f"{rotation[0]} {rotation[1]} {rotation[2]}\n")

    print("\n" + "Generated files in " + dir)
    print("Random seed: " + str(seed))


if __name__ == "__main__":
    main()
