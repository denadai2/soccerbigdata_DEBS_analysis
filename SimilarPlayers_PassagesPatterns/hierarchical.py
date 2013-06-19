
import dataset
import numpy as np
import scipy.cluster
from pprint import pprint
from scipy.cluster.hierarchy import dendrogram
import matplotlib.pyplot as plt
from scipy.spatial import distance

db = dataset.connect('sqlite:///prova.sqlite')
all_points=[]
playerLabels = []

for player in db['dbscan']:
    if player['name'] != 'Referee':
        all_points.append([player['unsuccessfulPassages'], player['unsuccessfulPassagesDefensiveThird']])#, player['unsuccessfulCrosses']])
        playerLabels.append(player['name'])


    #    mlabOrder = N.array(1,2,6,7,3,8,9,4,5,10) #the order returned by the original MATLAB code
    # Remeber MATLAB counts from 1, python from 0

all_points = np.array(all_points)
thresh = 4

clusters =  scipy.cluster.hierarchy.fclusterdata(all_points, thresh, criterion='distance', metric='euclidean', depth=2, method='single', R=None)

title = "threshold: %f, number of clusters: %d" % (thresh, len(set(clusters)))
print title

pprint(clusters)
prettyClusters = [[[] for i in range(0)] for i in range(6)]
for (i, item) in enumerate(clusters):
    print "CLUSTER "+str(item)+" :"+str(all_points[i])
    prettyClusters[item-1].append(all_points[i])

pprint(prettyClusters)


def augmented_dendrogram(*args, **kwargs):

    ddata = dendrogram(*args, **kwargs)

    if not kwargs.get('no_plot', False):
        for i, d in zip(ddata['icoord'], ddata['dcoord']):
            x = 0.5 * sum(i[1:3])
            y = d[1]
            plt.plot(x, y, 'ro')
            plt.annotate("%.3g" % y, (x, y), xytext=(0, -8),
                         textcoords='offset points',
                         va='top', ha='center')

    return ddata


distanceMatrix = distance.pdist(all_points,'euclidean') #use euclidean function
distanceSquareMatrix = distance.squareform(distanceMatrix)
linkageMatrix = scipy.cluster.hierarchy.linkage(distanceSquareMatrix)
augmented_dendrogram(linkageMatrix,
           color_threshold=1,
           show_leaf_counts=True,
           labels=playerLabels)
plt.show()