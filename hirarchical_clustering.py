# Hierarchical (agglomerative) clustering 
# As an example of an unsupervised learning problem: 
# classify it in a number of clusters derived from the dendogram
# compare with K-means clustering results

#----------------------------
## written by Sadra Dec 2019
#----------------------------
import numpy as np
import pandas as pd
from numpy import savetxt
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as shc
from sklearn.cluster import AgglomerativeClustering
from sklearn import metrics
from sklearn.metrics import silhouette_score
from sklearn.metrics import davies_bouldin_score
from sadra import pubfig as pf
pf.setup(label_font=24, tick_font=22, axis_width=2, tick_major_width=2,
         tick_minor_width=1.5, tick_major_size=4, tick_minor_size=4, showminorticks=True)
# plt.figure(figsize=(10,4))
qn_path = '/media/data01/26s-S3/atpaseparts/gomodel/2gfp/analysis/editedparam/2gfp-rep-force/repetitve_random-qn-rmsd-rg-2gfp-26s-l51-T300-f100/repetitive-rand_pull-6atpase-nocyc400'
qn_ngfp = []
qn_gfpc = []
# for traj in range(1, 6):
for traj in range(1,67):
    # print(traj)
    qn_address = qn_path +"/qn-rmsd-2gfp-26s-repetitive-l51-f100-traj%d/qn-rept-2gfp.dat.merged" % traj
    qn_file = np.genfromtxt(qn_address, comments='#')
    ### make new file to save every 10 data
    qn_rows = np.shape(qn_file)[0]  # number of rows
    qn_columns = np.shape(qn_file)[1]  # number of columns
    # index of each 10 elements including the first element
    qn_index = range(qn_rows)[0::10]
    new_qn_file = np.ndarray(
        shape=(len(qn_index), qn_columns))  # my new matrix
    # print(new_qn_file)
    for i in range(len(qn_index)):
        # adding each 10th row from qn_file to the new qn_file
        new_qn_file[i] = qn_file[qn_index[i]]
    qn_ngfp_data = new_qn_file[:, 3]
    qn_gfpc_data = new_qn_file[:, 4]
    # print(len(qn_gfpc_data))
    qn_ngfp.extend(qn_ngfp_data)
    qn_gfpc.extend(qn_gfpc_data)
qn_ngfp = np.array(qn_ngfp)
qn_gfpc = np.array(qn_gfpc)
X = np.array(list(zip(qn_ngfp, qn_gfpc)))
# we need to know the clusters that we want our data to be split to: create the dendogram
# we use the single, or the complete, or the complete method for linkage
# plt.figure(figsize=(10, 7))
# plt.title("Customer Dendogram")
"""
#this is for single linkage
dend = shc.dendrogram(shc.linkage(data, method='single'))
plt.savefig('dendrogram_shoppingdata_single.cols3to4.png')
plt.show()
"""
"""
#this is for complete linkage
dend = shc.dendrogram(shc.linkage(data, method='complete'))
plt.savefig('dendrogram_shoppingdata_complete.cols3to4.png')
plt.show()
"""
#this is for complete linkage
# dend = shc.dendrogram(shc.linkage(data, method='complete'))
# plt.savefig('dendrogram_shoppingdata_complete.cols3to4.png')
# plt.show()
plt.figure(figsize=(10, 7))
plt.title("Dendrogram")
data = X
#this is for complete linkage
# dend = shc.dendrogram(shc.linkage(data, method='complete'))
dend = shc.dendrogram(shc.linkage(data, method='complete'),truncate_mode='lastp',p=12, leaf_rotation=45., leaf_font_size=15., show_contracted=True)
plt.savefig('dendrogram.png')
plt.show()
# # use the Calinski-Harabasz score to determine number of clusters
pSF = []
linkage_type = "ward"
for number_clust in range(2, 6):
    cluster = AgglomerativeClustering(n_clusters=number_clust, affinity='euclidean', linkage=linkage_type)
    y = cluster.fit_predict(data)
    labels = cluster.labels_
    pSF.append(metrics.calinski_harabasz_score(data, labels))
plt.plot(range(2, 6), pSF)
plt.plot(range(2, 6), pSF)
plt.scatter(range(2, 6), pSF, c = 'blue')
plt.title('Calinski Harabasz Score vs. clusters (%s Linkage)'%linkage_type)
plt.xlabel('Number of clusters')
plt.ylabel('Calinski-Harabasz')
plt.savefig('CalinskiHarabasz_Agglom-%s_linkage.png'%linkage_type, dpi=300, bbox_inches="tight")
# plt.show()
plt.close()
# use the silhouette score to determine the number of clusters
silh = []
for number_clust in range(2, 6):
    cluster = AgglomerativeClustering(n_clusters=number_clust, affinity='euclidean', linkage=linkage_type)
    y = cluster.fit_predict(data)
    labels = cluster.labels_
    silh.append(metrics.silhouette_score(data, labels, metric='euclidean'))
    message = "For n_clusters = {} The complete silhouette_score is: {}"
    print(message.format(number_clust, round(metrics.silhouette_score(data, labels, metric='euclidean'), 2)))    
plt.plot(range(2, 6), silh)
plt.plot(range(2, 6), silh)
plt.scatter(range(2, 6), silh, c = 'blue')
plt.title('Silhouette Score vs. clusters (%s linkage)'%linkage_type)
plt.xlabel('Number of clusters')
plt.ylabel('Silhouette') 
plt.savefig('AvSilhouette_Agglom-%s_linkage.png'%
            linkage_type, dpi=300, bbox_inches='tight')
# plt.show()
plt.close()
# calculate the Davies Bouldin score (DBI) for agglomerative clustering
scores = []
centers = list(range(2,6))
for center in centers:
    cluster = AgglomerativeClustering(n_clusters=center, affinity='euclidean', linkage=linkage_type)
    y = cluster.fit_predict(data)
    labels = cluster.labels_
    scores.append(metrics.davies_bouldin_score(data, labels))
plt.plot(centers, scores, linestyle='--', marker='o', color='b')
plt.xlabel('clusters')
plt.ylabel('Davies Bouldin score')
#plt.title('Davies Bouldin score vs. clusters (Single linkage)')
#plt.title('Davies Bouldin score vs. clusters (complete linkage)')
plt.title('Davies Bouldin score vs. clusters (%s linkage)'%linkage_type)
#plt.savefig('DBIscore_Agglomsingle_shopping.png')
#plt.savefig('DBIscore_Agglomcomplete_shopping.png')
plt.savefig('DBIscore_Agglom-%s_linkage.png' % linkage_type, dpi=300, bbox_inches='tight')
# plt.show()
plt.close()
# we know the number of clusters for complete linkage for our dataset is 4 or 5 (from most scores it is 5)
# the next step is to group the data points into these 5 clusters
for n_cluster in range(2,6):
    cluster = AgglomerativeClustering(n_clusters=n_cluster, affinity='euclidean', linkage=linkage_type)
    cluster.fit_predict(data)
    #save cluster labels in a file
    savetxt('clusters.csv', cluster.labels_, fmt = '%d', delimiter=',')
    # plot the clusters in the 2D space
    # plt.figure(figsize=(10, 7))
    plt.scatter(data[:,0], data[:,1], c=cluster.labels_, cmap='rainbow')
    plt.xlabel(r'$Q_N$ (GFP-C)')
    plt.ylabel(r'$Q_N$ (N-GFP)')
    plt.title("Agglomerative clustering with %s linkage with %d clusters"%(linkage_type, n_cluster))
    image_name = 'agglomerativeclusters%s-linkage-%d-clusters.png'%(linkage_type, n_cluster)
    plt.savefig(image_name, dpi=300, bbox_inches='tight' )
    # plt.show()
    plt.close()
