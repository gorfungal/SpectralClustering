import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

from sklearn import cluster
from sklearn import datasets

def calculate_distance_from_center(point,center):
    distance = float(0)
    for i in range(0,len(point)):
        distance += (point[i] - center[i])**2
    return distance

def k_means_pp_bonus(points_array,K):
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

def main():
    iris = datasets.load_iris()
    X = iris.data
    ax = np.array([i for i in range(10)])
    values = np.array([0 for i in range(10)])
    for i in range(1,11):
        for_i = cluster.KMeans(n_clusters=i, init=k_means_pp_bonus(X, i), n_init=1, random_state=0)
        kmeans = for_i.fit(X)
        values[i-1]=kmeans.inertia_
    
    figure, axes = plt.subplots()
    plt.annotate('Elbow Point', xy =(2, 80),
                xytext =(3, 1.8), 
                arrowprops = dict(facecolor ='black',
                                  shrink = 0.05),)
    circle1 = Ellipse((2.1, 80),  width=0.5, height=50,linestyle='--' ,color='black',fill=False)
    axes.add_artist(circle1)
    plt.plot(ax, values)
    plt.xlabel('K')
    plt.ylabel('inertia')

    plt.title('inertia / num of clusters')
    plt.savefig(fname='elbow.png')

main()