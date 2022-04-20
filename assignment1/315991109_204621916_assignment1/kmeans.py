__author__ = "Liad Panker, Gal Gorfung"
__copyright__ = "Copyright 2021,Project C"
__version__ = "1.0.0"
__maintainer__ = "Liad Panker"
__email__ = "Liadpanker@mail.tau.ac.il"
__status__ = "Production"
import sys

espillon = (0.001)**2

class cluster:
    def __init__(self, center):
        self.center = center
        self.point_cord_num = len(self.center)
        self.sum = [float(0) for i in range(self.point_cord_num)]
        self.points_len = 0
        self.less_then_espillon = False
        self.distance = float(0)
    def cal_new_center(self):
        if(self.points_len != 0):
            new_center = [self.sum[i]/float(self.points_len) for i in range(0,self.point_cord_num)]
        else:
            new_center = self.center
        euclidean_norm = float(0)
        for i in range(0,self.point_cord_num):
            euclidean_norm += (new_center[i] - self.center[i]) ** 2
        self.less_then_espillon = euclidean_norm < espillon
        self.center = new_center
    def calculate_distance_from_center(self,point):
        self.distance = float(0)
        for i in range(0,self.point_cord_num):
            self.distance += (point[i] - self.center[i])**2
        return self.distance
    def reset_points_array(self):
        self.points_len = float(0)
        self.sum = [float(0) for i in range(0,self.point_cord_num)]
    def add_point(self,point):
        self.points_len += 1
        self.sum = [self.sum[i] + point[i] for i in range(0,self.point_cord_num)]
    def print_center(self,fd):
        for i in range(0,self.point_cord_num):
            if(i != self.point_cord_num-1):
                print("%.4f" % self.center[i], file=fd, end=",")
            else:
                print("%.4f" % self.center[i], file=fd)
#read the array fron file
def read_points():
    points_array = []
    try:
        fd = open(input_file, 'r')
    except:
        print("Invalid Input!")
    lines = fd.readlines()
    for i in range(0,len(lines)):
        points_array.append(lines[i][:-1].split(','))
        points_array[-1] = [float(cord) for cord in points_array[-1]]
    fd.close()
    return points_array

def export_results(cluster_array):
    try:
        fd = open(output_file, 'w')
    except:
        print("Invalid Output!")
    for clust in cluster_array:
        clust.print_center(fd)
    fd.close()


def clasify_to_clusters(cluster_array,points_array):
    [cluster_array[i].reset_points_array() for i in range(0,K)]
    for point in points_array:
        cluster_array[0].calculate_distance_from_center(point)
        min_dis_cluster_i = 0
        for i in range(1,K):
            if(cluster_array[i].calculate_distance_from_center(point) <= cluster_array[min_dis_cluster_i].distance):
                min_dis_cluster_i = i
        cluster_array[min_dis_cluster_i].add_point(point)
    [cluster_array[i].cal_new_center() for i in range(0, K)]



#main
def main():
    args = sys.argv[1:]
    global K
    global max_iter
    global input_file
    global output_file
    if(len(args) > 4 or len(args) < 3):
        print("Invalid Input!")
        exit(1)
    if(len(args) == 4):
        if (not args[0].isdigit()) or (not args[1].isdigit()):
            print("Invalid Input!")
            exit(1)
        K = int(args[0])
        max_iter = int(args[1])
        input_file = args[2]
        output_file = args[3]
    if(len(args) == 3):
        if (not args[0].isdigit()):
            print("Invalid Input!")
            exit(1)
        K = int(args[0])
        max_iter = 200
        input_file = args[1]
        output_file = args[2]
    points_array = read_points()
    cluster_array = [cluster(points_array[i]) for i in range(0,K)]
    num_of_iter = 1
    less_then_espillon = [cluster_array[i].less_then_espillon for i in range(0,K)]
    while((num_of_iter <= max_iter) and not all(less_then_espillon)):
        clasify_to_clusters(cluster_array,points_array)
        less_then_espillon = [cluster_array[i].less_then_espillon for i in range(0, K)]
    export_results(cluster_array)
main()





        
