import dataset
import numpy as np
from pprint import pprint
import kmeans
import sys
import select
from sklearn.cluster import KMeans
from scipy.spatial import distance
import matplotlib.pyplot as plt

import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D


db = dataset.connect('sqlite:///prova.sqlite')
all_points=[]
playerLabels = []

true_k = raw_input("How many groups (clusters) do you want to form?")

for player in db['dbscan']:
    if player['name'] != 'Referee':
        all_points.append([player['successfulPassages'],
                           player['unsuccessfulPassages'],
                           player['unsuccessfulPassagesDefensiveThird'],
                           player['successfulCrosses'],
                           player['unsuccessfulCrosses'],
                           player['shotsOnGoal'],
                           player['shotsFailed'],
                           player['wonTackles'],
                           player['lostTackles']])
        playerLabels.append(player['name'])


D = distance.squareform(distance.pdist(all_points, 'cosine'))

all_points = np.array(all_points)

km = KMeans(n_clusters=int(true_k), init='k-means++', max_iter=300, verbose=1)
km.fit(all_points)

clusters = km.labels_

'''fig = plt.figure()
axPl = fig.add_subplot(111, projection='3d')'''
prettyClusters = [[[] for i in range(0)] for i in range(int(true_k))]
prettyLabels = [[[] for i in range(0)] for i in range(int(true_k))]
for (i, item) in enumerate(clusters):

    prettyClusters[item].append(all_points[i])
    prettyLabels[item].append(playerLabels[i])

    if item == 0:
        color = "b"
    elif item == 1:
        color = "r"
    elif item == 2:
        color = "m"
    elif item == 3:
        color = "y"
    else:
        color = "g"

    #axPl.scatter(all_points[i][0], all_points[i][1], all_points[i][2], c=color, marker='o')

for (i,cluster) in enumerate(prettyClusters):
    print "CLUSTER "+str(i)+" :"+str(prettyLabels[i]) + " MEDIA: " + str(np.mean(cluster, 0))




plt.show()

