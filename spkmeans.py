__author__ = "Liad Panker, Gal Gorfung"
__copyright__ = "Copyright 2022,Project C"
__version__ = "1.0.0"
__maintainer__ = "Liad Panker"
__email__ = "Liadpanker@mail.tau.ac.il"
__status__ = "Production"
import sys
import numpy as np
import mykmeanssp as kmeans
import spkmeansmodule as spk
#parameters
eps = 0
goal_enum = ["spk","wam","ddg","lnorm","jacobi"]

#read the array fron file
def read_points():
    points = 0
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

def kmeans(v_mat):
    points_array = read_points(v_mat)
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

#main
def main():
    args = sys.argv[1:]
    global K
    global goal
    global input1_file
    if(len(args) > 3 or len(args) < 2):
        print("Invalid Input!")
        exit(1)
    K = int(args[0])
    goal = args[1]
    if(goal not in goal_enum):
        print("Invalid Input!")
    input1_file = args[2]
    mat = spk.spkmeans(goal,input1_file)
    if(goal == "spk"):
        kmeans(mat)
        exit(0)
    else:
        print(mat)
    
main()





        
