__author__ = "Liad Panker, Gal gorfung"
__copyright__ = "Copyright 2021,Project C"
__version__ = "1.0.0"
__maintainer__ = "Liad Panker"
__email__ = "Liadpanker@mail.tau.ac.il"
__status__ = "Production"
import sys
import numpy as np
import mykmeanssp as kmeans


#read the array fron file
def read_points():
    global D
#open file .txt or csv
    try:
        fd1 = open(f'{input1_file}.txt', 'r')
    except:
        try:
            fd1 = open(f'{input1_file}.csv', 'r')
        except:
            print("Invalid Input!")
            exit(1)
    try:
        fd2 = open(f'{input2_file}.txt', 'r')
    except:
        try:
            fd2 = open(f'{input2_file}.csv', 'r')
        except:
            print("Invalid Input!")
            exit(1)
    lines_1 = fd1.readlines()
    lines_2 = fd2.readlines()
#reading len
    len_line_1 = len(lines_1)
    len_line_2 = len(lines_2)
    if(len_line_1 != len_line_2):
        print("Invalid Input!")
    if(len_line_1 == 0):
        print("Invalid Input!")
#reading num of points cord
    half_cord_1 = len(lines_1[0][:-1].split(',')) -1 #calculate the half dim
    half_cord_2 = len(lines_1[0][:-1].split(',')) -1 #calculate the half dim
    D = half_cord_1 + half_cord_2
    points = np.array([np.array([0 for j in range(half_cord_1+half_cord_2)]).astype('float') for i in range(len_line_1)]) # define the points array
#reading points to the array
    for i in range(0,len_line_1):
        line_1 = lines_1[i][:-1].split(',')
        line_2 = lines_2[i][:-1].split(',')
        for cord_i in range(half_cord_1):
            points[int(float(line_1[0]))][cord_i] = float(line_1[cord_i+1])
        for cord_i in range(half_cord_2):
            points[int(float(line_2[0]))][cord_i+half_cord_1] = float(line_2[cord_i+1])
    fd1.close()
    fd2.close()
    return points

def calculate_distance_from_center(point,center):
    distance = float(0)
    for i in range(0,len(point)):
        distance += (point[i] - center[i])**2
    return distance


def k_means_pp(points_array):
    global centers_ind
    i = 1
    np.random.seed(0)
    len_point = len(points_array[0])
    len_points_array = len(points_array)
    centers = np.array([np.array([0 for j in range(len_point)]).astype('float') for k in range(K)])
    centers_ind = np.array([-1 for k in range(K)])
    prob_array = np.array([0 for j in range(len_points_array)]).astype('float')
    d_l_array = np.array([0 for j in range(len_points_array)]).astype('float')
    # random point
    center_point_i = np.random.choice(len_points_array, 1)
    while(i <= K):
        centers_ind[i - 1] = center_point_i[0]
        center_point = points_array[center_point_i[0]]
        centers[i - 1] = center_point
        #calculate_dl
        dl_sum = 0
        for point_i in range(0,len_points_array):
            min_distance = calculate_distance_from_center(points_array[point_i],centers[0])
            for j in range(1,i):
                distance_new = calculate_distance_from_center(points_array[point_i],centers[j])
                if(distance_new<min_distance):
                    min_distance = distance_new
            d_l_array[point_i] = min_distance
            dl_sum += min_distance
        #calcualte probability
        for prob_i in range(0,len_points_array):
            prob_array[prob_i] = d_l_array[prob_i]/dl_sum
        i+=1
        # random point
        center_point_i = np.random.choice(len_points_array, 1, p = prob_array)
    return centers


#main
def main():
    args = sys.argv[1:]
    global K
    global max_iter
    global input1_file
    global input2_file
    global output_file
    global D
    if(len(args) > 5 or len(args) < 4):
        print("Invalid Input!")
        exit(1)
    if(len(args) == 5):
        K = int(args[0])
        max_iter = int(args[1])
        eps = float(args[2])
        input1_file = args[3]
        input2_file = args[4]
    if(len(args) == 4):
        K = int(args[0])
        max_iter = 200
        eps = float(args[1])
        input1_file = args[2]
        input2_file = args[3]
    points_array = read_points()
    centers = k_means_pp(points_array)
 #calling kmeans function
    Points = points_array.flatten().astype('float64').tolist()
    Centers = centers.flatten().astype('float64').tolist()
    e = eps
    MAX = max_iter
    lg = K
    N = len(points_array)
    for center_i in centers_ind[:-1]:
        print(center_i,end=',')
    print(centers_ind[-1])
    centers_fin = np.array(kmeans.fit(N, lg, D, MAX, e,Centers,Points))
    centers_fin = centers_fin.reshape(K,D)
    for i in range(K):
        for j in range(D):
            if(j != D-1):
                print("%.4f" % centers_fin[i][j], end=",")
            else:
                print("%.4f" % centers_fin[i][j])



main()





        
