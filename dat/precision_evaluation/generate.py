import math
import random
import os
import sys
import argparse

def cube(new_file):
    for i in range(200):
        new_file.write(str(random.uniform(-50,50)) + "\t" + str(random.uniform(-50,50)) + "\t" + "100" + "\n")
        new_file.write(str(random.uniform(-50,50)) + "\t" + str(random.uniform(-50,50)) + "\t" + "150" + "\n")
        new_file.write(str(random.uniform(-50,50)) + "\t" + "-50" + "\t" + str(random.uniform(50, 150)) + "\n")
        new_file.write(str(random.uniform(-50,50)) + "\t" + "50" + "\t" + str(random.uniform(50, 150)) + "\n")
        new_file.write("-50" + "\t" + str(random.uniform(-50,50)) + "\t" + str(random.uniform(50, 150)) + "\n")
        new_file.write("50" + "\t" + str(random.uniform(-50,50)) + "\t" + str(random.uniform(50, 150)) + "\n")

def volume_cube(new_file):
    for i in range(1200):
        new_file.write(str(random.uniform(-50, 50)) + "\t" + str(random.uniform(-50, 50)) + "\t" + str(random.uniform(50, 150)) + "\n")

def sphere(new_file): # https://mathworld.wolfram.com/SpherePointPicking.html
    # Radius
    r = 50
    for i in range(1200):
        # Polar coordinates
        #theta = random.uniform(0, math.pi)
        #phi = random.uniform(0, 2*math.pi)
        u = random.uniform(0, 1)
        v = random.uniform(0, 1)

        phi = 2 * math.pi* u
        theta = math.acos(2*v-1)

        # Polar coordinates to xyz values
        new_file.write(str(r * math.sin(theta) * math.cos(phi)) + "\t" + str(r * math.sin(theta) * math.sin(phi)) + "\t" + str(r * math.cos(theta)+125) + "\n")

def volume_sphere(new_file): # https://codepen.io/thebabydino/pen/abLRNqV
    for i in range(1200):
        # Radius
        r = random.uniform(0, 1)
        r = math.pow(r, 1/3) * 50
        # Polar coordinates
        #theta = random.uniform(0, math.pi)
        #phi = random.uniform(0, 2*math.pi)

        u = random.uniform(0, 1)
        v = random.uniform(0, 1)

        phi = 2 * math.pi* u
        theta = math.acos(2*v-1)

        # Polar coordinates to xyz values
        new_file.write(str(r * math.sin(theta) * math.cos(phi)) + "\t" + str(r * math.sin(theta) * math.sin(phi)) + "\t" + str(r * math.cos(theta)+125) + "\n")



def main():
    # Evaluate arguments

    parser = argparse.ArgumentParser(prog = 'generate', description='Generate point cloud in different shapes')
    parser.add_argument("-d", "--directory", required = True, help="output directory")
    parser.add_argument("-s",  "--shape", required = True, choices = ['cube', 'volume_cube', 'sphere', 'volume_sphere'],  help = "Shape. Possible shapes: \n \t cube \n \t volume_cube \n \t sphere \n \t volume_sphere" )

    args = parser.parse_args()
    dir = args.directory
    shape = args.shape


    #Create scan:
    if not os.path.exists(dir):
        parser.error("Path does not exist")
    with open(dir + "/scan000.3d", "w", encoding="utf-8") as new_file:

        match (shape):
            case "cube":
                cube(new_file)
            case "volume_cube":
                volume_cube(new_file)
            case "sphere":
                sphere(new_file)
            case "volume_sphere":
                volume_sphere(new_file)
            case _:
                parser.error("Invalid shape argument")
    print("\n" + "Generated file \"scan000.3d\" in " + dir)


if __name__ == "__main__":
    main()
