from importlib.util import set_package
import sys
import csv
import enum

def main():
    if len(sys.argv) != 4:
        print("Invalid Input!")
        sys.exit()
    k = int(sys.argv[1])
    goal = sys.argv[2]
    file_name = sys.argv[3]
    len = 0
    points = []
    if file_name[:-4] == ".csv":
        with open(file_name, newline='') as csvfile:
            reader = csv.reader(csvfile, delimiter=',')
            for row in reader:
                p = []
                for e in row:
                    p.append(float(e))
                points.append(p)
                len += 1
    elif file_name[:-4] == ".txt":
        with open(file_name,"r") as txt:
            for line in txt.readlines():
                p = []
                for e in line:
                    p.append(float(e))
                points.append(p)
                len += 1
    else:
        print("Invalid Input!")
        sys.exit()



