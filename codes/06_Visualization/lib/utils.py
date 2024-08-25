import numpy as np
import scipy.sparse
from collections import defaultdict
from scipy import stats
import ncut
import sklearn.metrics.pairwise
import time

import matplotlib.pyplot as plt
import time

import torch
import torch.nn as nn

import psutil




######################################
# Usages: 
# (1) W = construct_knn_graph(X,k,'euclidean')
# (2) W = construct_knn_graph(X,k,'euclidean_zelnik_perona')
# (3) W = construct_knn_graph(X,k,'cosine')
# (4) W = construct_knn_graph(X,k,'cosine_binary')
#
# X = torch.Tensor([[1,2,9],[4,5,2],[7,5,9],[0,11,12]])
# print(X, X.size())
# kNN = 3
# batch_size = 2
# W = construct_knn_graph(X, kNN, 'euclidean', batch_size)
# print(W)
#
# Notations:
#   n = nb_data
#   d = data_dimensionality
#
# Input variables:
#   X = Data matrix. Size = n x d.
#   kNN = Number of nearest neaighbors. 
#   dist = Type of distances: 'euclidean' or 'cosine'
#
# Output variables:
#   W = Adjacency matrix. Size = n x n.
######################################

def construct_knn_graph(X, kNN, dist, batch_size=2000, device='cpu'):

    batch_size = int( batch_size * 16 / (psutil.virtual_memory().total/2**30) ) # ram mem = psutil.virtual_memory().total/2**30
    print('batch_size',batch_size)

    if torch.is_tensor(X) == False:
        print('Convert from numpy to pytorch')
        X = torch.Tensor(X) # convert from numpy to pytorch
        
    X = X.to(device)

    # Compute L2/Euclidean distance between pairs of points with matrix-matrix mulitplications
    # Dist_sqr = |xi-xj|^2 = |xi|^2 + |xj|^2 - 2 xi.xj
    # Dist_sqr = [ D11 D12 Dij ... D1m ]
    #            [ D21 D22     ... D2m ]
    #            [             ...     ]
    #            [ Dn1 Dn2     ... Dnm ]
    #
    #          = [ |x1|^2 |x1|^2 ... |x1|^2 ]
    #            [ |x2|^2 |x2|^2 ... |x2|^2 ]
    #            [               ...        ]
    #            [ |xn|^2 |xn|^2 ... |xn|^2 ]
    #          +
    #            [ |x1|^2 |x2|^2 ... |xn|^2 ]
    #            [ |x1|^2 |x2|^2 ... |xn|^2 ]
    #            [               ...        ]
    #            [ |x1|^2 |x2|^2 ... |xn|^2 ]
    #      - 2 *
    #            [ x1.x1 x1.x2 ... x1.xn ]
    #            [ x2.x1 x2.x2 ... x2.xn ]
    #            [             ...       ]
    #            [ xn.x1 xn.x2 ... xn.xn ]
  
    ######################################
    # Construct a k-NN graph with L2/Euclidean distance
    ######################################
    if dist == 'euclidean':
        
        print('k-NN graph construction with euclidean distance')
        start = time.time()

        # Compute L2 distance between all pairs of points
        X = X - X.mean(dim=0).unsqueeze(dim=0) # centered data, size=[n, d] 
        X2 = (X**2).sum(dim=1) # size=[n] 
        n = X.size(0)
        num_batch = n // batch_size
        if (n//batch_size)*batch_size < n: num_batch += 1
        print('n, batch_size, num_batch, kNN, device:', n, batch_size, num_batch, kNN, device)
        values_knn = []
        idx_knn = []
        idx_current = 0
        for idx in range(num_batch):
            X_sliced = X[idx*batch_size:(idx+1)*batch_size, :] # size=[batch_size, d] 
            batch_size_current = X_sliced.size(0) 
            term1 = X2.repeat(batch_size_current,1) # size=[batch_size, n]
            term2 = X2[idx_current:idx_current+batch_size_current].unsqueeze(1).repeat(1,n) # size=[batch_size, n]
            idx_current += batch_size_current
            term3 = X_sliced @ X.transpose(1,0) # size=[batch_size, n] 
            D = term1 + term2 - 2 * term3 # size=[batch_size, n]
            sorted, indices = torch.sort(D, dim=1) # from smallest to largest values, size=[batch_size, n]
            values_knn.append(sorted[:,:kNN]) # size=[batch_size, kNN]
            idx_knn.append(indices[:,:kNN]) # size=[batch_size, kNN]
            del X_sliced, term1, term2, term3, D, sorted, indices
            print('batch idx:',idx+1,'/',num_batch,' time(sec):',(time.time()-start)/1) 
        values_knn = torch.cat(values_knn)
        idx_knn = torch.cat(idx_knn)

        # Compute Weight matrix
        D = values_knn.sqrt()
        sigma2 = (D[:,-1]).mean()**2 # graph-level scale
        D = (D**2).view(n*kNN)
        W = torch.exp( -D / sigma2 )
        
        # Make W sparse
        row = torch.arange(n).repeat_interleave(kNN).to('cpu')
        col = idx_knn.view(n*kNN).to('cpu')
        data = W.to('cpu')
        W = scipy.sparse.csr_matrix((data, (row, col)), shape=(n, n))
        
        # Make W is symmetric
        bigger = W.T > W
        W = W - W.multiply(bigger) + W.T.multiply(bigger)
        W.setdiag(0) # remove self-connections

        print('Computational time(sec) for k-NN graph:',(time.time()-start)/1) 

    
    ######################################
    # Construct a k-NN graph with Zelnik-Perona technique
    # "Self-Tuning Spectral Clustering", 2005
    ######################################
    if dist == 'euclidean_zelnik_perona':

        print('k-NN graph construction with Zelnik-Perona technique')
        start = time.time()

        # Compute L2 distance between all pairs of points
        X = X - X.mean(dim=0).unsqueeze(dim=0) # centered data, size=[n, d] 
        X2 = (X**2).sum(dim=1) # size=[n] 
        n = X.size(0)
        num_batch = n // batch_size
        if (n//batch_size)*batch_size < n: num_batch += 1
        print('n, batch_size, num_batch, kNN, device:', n, batch_size, num_batch, kNN, device)
        values_knn = []
        idx_knn = []
        idx_current = 0
        for idx in range(num_batch):
            X_sliced = X[idx*batch_size:(idx+1)*batch_size, :] # size=[batch_size, d] 
            batch_size_current = X_sliced.size(0) 
            term1 = X2.repeat(batch_size_current,1) # size=[batch_size, n]
            term2 = X2[idx_current:idx_current+batch_size_current].unsqueeze(1).repeat(1,n) # size=[batch_size, n]
            idx_current += batch_size_current
            term3 = X_sliced @ X.transpose(1,0) # size=[batch_size, n] 
            D = term1 + term2 - 2 * term3 # size=[batch_size, n]
            sorted, indices = torch.sort(D, dim=1) # from smallest to largest values, size=[batch_size, n]
            values_knn.append(sorted[:,:kNN]) # size=[batch_size, kNN]
            idx_knn.append(indices[:,:kNN]) # size=[batch_size, kNN]
            del X_sliced, term1, term2, term3, D, sorted, indices
            print('batch idx:',idx+1,'/',num_batch,' time(sec):',(time.time()-start)/1) 
        values_knn = torch.cat(values_knn)
        idx_knn = torch.cat(idx_knn)

        # Compute Weight matrix
        D = values_knn.sqrt()
        sigma = D[:,-1] 
        sigma_i = sigma[torch.arange(n).repeat_interleave(kNN)] # sigma_i in Zelnik-Perona technique
        sigma_j = sigma[idx_knn.view(n*kNN)] # sigma_j in Zelnik-Perona technique
        D = (D**2).view(n*kNN) / ( sigma_i * sigma_j )
        W = torch.exp(-D)
        
        # Make W sparse
        row = torch.arange(n).repeat_interleave(kNN).to('cpu')
        col = idx_knn.view(n*kNN).to('cpu')
        data = W.to('cpu')
        W = scipy.sparse.csr_matrix((data, (row, col)), shape=(n, n))
        
        # Make W is symmetric
        bigger = W.T > W
        W = W - W.multiply(bigger) + W.T.multiply(bigger)
        W.setdiag(0) # remove self-connections

        print('Computational time(sec) for k-NN graph:',(time.time()-start)/1) 

    
    ######################################
    # Construct a k-NN graph with Cosine distance
    ######################################
    if dist == 'cosine':
        
        print('k-NN graph construction with cosine distance')
        start = time.time()

        # Compute dot products between all pairs of points
        X = X - X.mean(dim=0).unsqueeze(dim=0) # centered data, size=[n, d] 
        X = X / (X**2).sum(dim=1).sqrt().unsqueeze(1) # normalized data / project data on unit sphere, size=[n, d] 
        n = X.size(0)
        num_batch = n // batch_size
        if (n//batch_size)*batch_size < n: num_batch += 1
        print('n, batch_size, num_batch, kNN, device:', n, batch_size, num_batch, kNN, device)
        values_knn = []
        idx_knn = []
        idx_current = 0
        for idx in range(num_batch):
            X_sliced = X[idx*batch_size:(idx+1)*batch_size, :] # size=[batch_size, d] 
            batch_size_current = X_sliced.size(0) 
            idx_current += batch_size_current
            D = X_sliced @ X.transpose(1,0) # size=[batch_size, n] 
            sorted, indices = torch.sort(D, dim=1, descending=True) # from *largest* to smallest values, size=[batch_size, n]
            sorted = torch.arccos(sorted).abs() # positive angles
            values_knn.append(sorted[:,:kNN]) # size=[batch_size, kNN]
            idx_knn.append(indices[:,:kNN]) # size=[batch_size, kNN]
            del X_sliced, D, sorted, indices
            print('batch idx:',idx+1,'/',num_batch,' time(sec):',(time.time()-start)/1) 
        values_knn = torch.cat(values_knn)
        idx_knn = torch.cat(idx_knn)

        # Compute Weight matrix
        D = values_knn.sqrt()
        sigma2 = (D[:,-1]).mean()**2 # graph-level scale
        D = (D**2).view(n*kNN)
        W = torch.exp( -D / sigma2 )

        # Make W sparse
        row = torch.arange(n).repeat_interleave(kNN).to('cpu')
        col = idx_knn.view(n*kNN).to('cpu')
        data = W.to('cpu')
        W = scipy.sparse.csr_matrix((data, (row, col)), shape=(n, n))
        
        # Make W is symmetric
        bigger = W.T > W
        W = W - W.multiply(bigger) + W.T.multiply(bigger)
        W.setdiag(0) # remove self-connections

        print('Computational time(sec) for k-NN graph:',(time.time()-start)/1) 

    
    ######################################
    # Construct a k-NN graph with Cosine distance
    ######################################
    if dist == 'cosine_binary':
        
        print('k-NN graph construction with cosine distance with binary weights')
        start = time.time()

        # Compute dot products between all pairs of points
        X = X - X.mean(dim=0).unsqueeze(dim=0) # centered data, size=[n, d] 
        X = X / (X**2).sum(dim=1).sqrt().unsqueeze(1) # normalized data / project data on unit sphere, size=[n, d] 
        n = X.size(0)
        num_batch = n // batch_size
        if (n//batch_size)*batch_size < n: num_batch += 1
        print('n, batch_size, num_batch, kNN, device:', n, batch_size, num_batch, kNN, device)
        values_knn = []
        idx_knn = []
        idx_current = 0
        for idx in range(num_batch):
            X_sliced = X[idx*batch_size:(idx+1)*batch_size, :] # size=[batch_size, d] 
            batch_size_current = X_sliced.size(0) 
            idx_current += batch_size_current
            D = X_sliced @ X.transpose(1,0) # size=[batch_size, n] 
            sorted, indices = torch.sort(D, dim=1, descending=True) # from *largest* to smallest values, size=[batch_size, n]
            sorted = torch.arccos(sorted).abs() # positive angles
            values_knn.append(sorted[:,:kNN]) # size=[batch_size, kNN]
            idx_knn.append(indices[:,:kNN]) # size=[batch_size, kNN]
            del X_sliced, D, sorted, indices
            print('batch idx:',idx+1,'/',num_batch,' time(sec):',(time.time()-start)/1) 
        values_knn = torch.cat(values_knn)
        idx_knn = torch.cat(idx_knn)

        # Compute Weight matrix, W sparse
        row = torch.arange(n).repeat_interleave(kNN).to('cpu')
        col = idx_knn.view(n*kNN).to('cpu')
        data = torch.ones(n*kNN).to('cpu')
        W = scipy.sparse.csr_matrix((data, (row, col)), shape=(n, n))
        
        # Make W is symmetric
        bigger = W.T > W
        W = W - W.multiply(bigger) + W.T.multiply(bigger)
        W.setdiag(0) # remove self-connections

        print('Computational time(sec) for k-NN graph:',(time.time()-start)/1) 
    
    return W





######################################
# Function that reindexes W according to communities/classes
######################################

######################################
# Usage: 
#   [reindexed_W,reindexed_C] = reindex_W_with_C(W,C)
#
# Notations:
#   n = nb_data
#   nc = nb_communities
#
# Input variables:
#   W = Adjacency matrix. Size = n x n.
#   C = Classes used for reindexing W. Size = n x 1. Values in [0,1,...,nc-1].
#
# Output variables:
#   reindexed_W = reindexed adjacency matrix. Size = n x n.
#   reindexed_C = reindexed classes C. Size = n x 1. Values in [0,1,...,nc-1].
######################################

def reindex_W_with_classes(W,C):
    n = C.shape[0] # nb of vertices
    nc = len(np.unique(C)) # nb of communities
    reindexing_mapping = np.zeros([n]) # mapping for reindexing W
    reindexed_C = np.zeros([n]) # reindexed C
    tot = 0
    for k in range(nc):
        cluster = (np.where(C==k))[0]
        length_cluster = len(cluster)
        x = np.array(range(tot,tot+length_cluster))
        reindexing_mapping[cluster] = x
        reindexed_C[x] = k
        tot += length_cluster
        
    idx_row,idx_col,val = scipy.sparse.find(W)
    idx_row = reindexing_mapping[idx_row]
    idx_col = reindexing_mapping[idx_col]
    reindexed_W = scipy.sparse.csr_matrix((val, (idx_row, idx_col)), shape=(n, n))

    return reindexed_W,reindexed_C



######################################
# Graph Laplacian Operator
######################################

######################################
# Usages: 
#   L = compute_graph_laplacian(W); # compute normalized graph Laplacian
#   L = compute_graph_laplacian(W,False); # compute UNnormalized graph Laplacian
#
# Notations:
#   n = nb_data
#
# Input variables:
#   W = Adjacency matrix. Size = n x n.
#
# Output variables:
#   L = Graph Laplacian. Size = n x n.
######################################

def graph_laplacian(W, normalized=True):
    
    # Degree vector
    d = W.sum(axis=0)

    # Laplacian matrix
    if not normalized:
        D = scipy.sparse.diags(d.A.squeeze(), 0)
        L = D - W
    else:
        d += np.spacing(np.array(0, W.dtype)) # d += epsilon
        d = 1.0 / np.sqrt(d)
        D = scipy.sparse.diags(d.A.squeeze(), 0)
        I = scipy.sparse.identity(d.size, dtype=W.dtype)
        L = I - D * W * D
    return L






######################################
# Graph Gradient Operator
######################################

######################################
# Usage: 
#   G = compute_graph_gradient(W); # compute normalized graph Laplacian
#
# Notations:
#   n = number of nodes
#   m = number of edges
#
# Input variables:
#   W = Adjacency matrix. Size = n x n.
#
# Output variables:
#   G = Graph Gradient Operator. Size = m x n.
######################################


def graph_gradient(W):

    W = W.todense()
    n = W.shape[0] # number of nodes
    Wtri = np.triu(W,1) # Extract upper triangular part of W 
    r,c = np.where(Wtri>0.0) # scipy.sparse.find
    v = Wtri[r,c]    
    ne = len(r)
    Dr = np.arange(0,ne); Dr = np.concatenate([Dr,Dr])
    Dc = np.zeros([2*ne], dtype='int32')
    Dc[:ne] = r
    Dc[ne:2*ne] = c
    Dv = np.zeros([2*ne])
    Dv[:ne] = np.sqrt(v)
    Dv[ne:2*ne] = -np.sqrt(v)
    G = scipy.sparse.csr_matrix((Dv, (Dr, Dc)), shape=(ne, n), dtype='float32')

    return G





######################################
# Visualization technique:
#   Belkin-Niyogi, Laplacian eigenmaps for dimensionality reduction and data representation, 2003
######################################

######################################
# Usage: 
#   X,Y,Z = nldr_visualization(W)
#
# Notations:
#   n = nb_data
#
# Input variables:
#   W = Adjacency matrix. Size = n x n.
#
# Output variables:
#   X = 1st data coordinates in low-dim manifold. Size n x 1.
#   Y = 2nd data coordinates in low-dim manifold. Size n x 1.
#   Z = 3rd data coordinates in low-dim manifold. Size n x 1.
######################################

def sortEVD(lamb, U):
    idx = lamb.argsort() # increasing order
    return lamb[idx], U[:,idx]

def nldr_visualization(W):
    
    # Compute normalized graph Laplacian
    L = graph_laplacian(W)
    
    # Regularization for ill-posed graphs
    L = L + 1e-6* scipy.sparse.identity(L.shape[0], dtype=W.dtype)

    # Compute the first three Laplacian Eigenmaps
    lamb, U = scipy.sparse.linalg.eigsh(L, k=4, which='SM')
    
    # Sort eigenvalue from smallest to largest values
    lamb, U = sortEVD(lamb, U)
    
    # Coordinates of graph vertices in the low-dim embedding manifold
    X = U[:,1]
    Y = U[:,2]
    Z = U[:,3]

    return X,Y,Z





######################################
# Clustering accuracy can be defined with the purity measure, defined here:
#   Yang-Hao-Dikmen-Chen-Oja, Clustering by nonnegative matrix factorization 
#   using graph random walk, 2012.
######################################

######################################
# Usages: 
#   accuracy = compute_clustering_accuracy(C_computed,C_grndtruth,R)
#
# Notations:
#   n = nb_data
#
# Input variables:
#   C_computed = Computed clusters. Size = n x 1. Values in [0,1,...,R-1].
#   C_grndtruth = Ground truth clusters. Size = n x 1. Values in [0,1,...,R-1].
#   R = Number of clusters.
#
# Output variables:
#   accuracy = Clustering accuracy of computed clusters.
######################################

def compute_purity(C_computed,C_grndtruth,R):
    N = C_grndtruth.size
    nb_of_dominant_points_in_class = np.zeros((R, 1))
    w = defaultdict(list)
    z = defaultdict(list)       
    for k in range(R):
        for i in range(N):
            if C_computed[i]==k:
                w[k].append(C_grndtruth[i])
        if len(w[k])>0:
            val,nb = stats.mode(w[k])
            z[k] = int(nb.squeeze()) 
        else:
            z[k] = 0
    sum_dominant = 0
    for k in range(R):
        sum_dominant = sum_dominant + z[k]
    purity = float(sum_dominant) / float(N)* 100.0
    return purity




######################################
# Graph spectral clustering technique NCut:
#   Yu-Shi, Multiclass spectral clustering, 2003 
#   Code available here: http://www.cis.upenn.edu/~jshi/software
######################################

######################################
# Usages: 
#   C,acc = compute_ncut(W,Cgt,R)
#
# Notations:
#   n = nb_data
#
# Input variables:
#   W = Adjacency matrix. Size = n x n.
#   R = Number of clusters.
#   Cgt = Ground truth clusters. Size = n x 1. Values in [0,1,...,R-1].
#
# Output variables:
#   C = NCut solution. Size = n x 1. Values in [0,1,...,R-1].
#   acc = Accuracy of NCut solution.
######################################

def compute_ncut(W, Cgt, R):

    # Apply ncut
    eigen_val, eigen_vec = ncut.ncut( W, R )
    
    # Discretize to get cluster id
    eigenvec_discrete = ncut.discretisation( eigen_vec )
    res = eigenvec_discrete.dot(np.arange(1, R + 1)) 
    C = np.array(res-1,dtype=np.int64)
    
    # Compute accuracy
    computed_solution = C
    ground_truth = Cgt
    acc = compute_purity( computed_solution,ground_truth, R)

    return C, acc












######################################
# Usages: 
#   [PC,PD,EnPD] = compute_pca(X,nb_pca)
#
# Notations:
#   n = nb_data
#   d = data_dimensionality
#
# Input variables:
#   X = Data matrix. Size = n x d.
#   nb_pca = Number of principal components. 
#
# Output variables:
#   PC = Principal components. Size = n x nb_pca.
#   PD = Principal directions. Size = d x nb_pca.
#   EnPD = Energy/variance of the principal directions. Size = np_pca x 1.
######################################

def sortPCA(lamb, U):
    idx = lamb.argsort()[::-1] # decreasing order
    return lamb[idx], U[:,idx]

def compute_pca(X,nb_pca):

    print('nb_pca:',nb_pca)
    
    Xzc = X - np.mean(X,axis=0) # zero-centered data
    
    n,d = X.shape
    if n>d:
        CovX = (Xzc.T).dot(Xzc) 
        dCovX = CovX.shape[0]
        # Compute the eigensolution
        #CovX = scipy.sparse.csr_matrix(CovX)
        if nb_pca<dCovX:
            lamb, U = scipy.sparse.linalg.eigsh(CovX, k=nb_pca, which='LM') # U = d x nb_pca
            lamb, U = sortPCA(lamb, U)  
            #U = U[:,::-1]
            #lamb = lamb[::-1]
        else: # nb_pca==dCovX
            lamb, U = np.linalg.eig(CovX)
            lamb, U = sortPCA(lamb, U)        
        PD = U # Principal directions
        PC = Xzc.dot(U) # Principal components/Projected data on PDs
        EnPD = lamb # Energy of PDs
        
    else:
        #Xzc = scipy.sparse.csr_matrix(Xzc)
        U,S,V = scipy.sparse.linalg.svds(Xzc, k=nb_pca, which='LM')
        U = U[:,::-1]
        S = S[::-1]
        V = V[::-1,:]
        PD = V # Principal directions
        PC = U.dot(np.diag(S)) # Principal components/Projected data on PDs
        EnPD = S**2 # Energy of PDs

    return PC,PD,EnPD








######################################
# Usages: 
#   [C_kmeans,En_kmeans] = compute_kernel_kmeans_EM(K,Ker,Theta,nb_trials)
#
# Note: Code based on Michael Chen (sth4nth@gmail.com)'s code
#
# Notations:
#   n = nb_data
#
# Input variables:
#   K = nb_clusters
#   Ker = Kernel matrix. Size = n x n.
#   Theta = Weight for each data term. Size = n x 1.
#   nb_trials = Number of kmeans runs.
#
# Output variables:
#   C_kmeans = Computed kmeans clusters. Size = n x 1.
#   En_kmeans = Energy of the kmeans partition.
######################################

def compute_kernel_kmeans_EM(nc,Ker,Theta,nb_trials):

    start = time.time()
    n = Ker.shape[0]
    Theta = np.diag(np.ones(n)) # Equal weight for each data
    Ones = np.ones((1,n))
    En_kmeans = 1e10; # Init energy value

    # Run K-Means "nb_trials" times and keep the solution that best minimizes 
    # the K-Means energy.
    for k in range(nb_trials):

        # Initialization
        C = np.random.randint(nc,size=n) # random initialization
        Cold = np.ones([n])
        diffC = np.linalg.norm(C-Cold)/np.linalg.norm(Cold)
    
        # Loop
        k = 0
        while (k<30) & (diffC>1e-2):
    
            # Update iteration
            k += 1
            #print(k)
    
            # Distance Matrix D
            row = np.array(range(n))
            col = C
            data = np.ones(n)
            F = scipy.sparse.csr_matrix((data, (row, col)), shape=(n, nc)).todense()
            O = np.diag( np.array( 1./ (Ones.dot(Theta).dot(F) + 1e-6) ).squeeze() )
            T = Ker.dot(Theta.dot(F.dot(O)))
            D = - 2* T + np.repeat( np.diag(O.dot((F.T).dot(Theta.dot(T))))[None,:] ,n,axis=0)
    
            # Extract clusters
            C = np.array(np.argmin(D,1)).squeeze()
                
            # L2 difference between two successive cluster configurations
            diffC = np.linalg.norm(C-Cold)/np.linalg.norm(Cold)
            Cold = C
        
        # K-Means energy
        En = np.multiply( (np.repeat(np.diag(Ker)[:,None],nc,axis=1) + D) , F)
        En = np.sum(En)/n


        # Check energy and keep the smallest one
        if En < En_kmeans:
            En_kmeans = En
            C_kmeans = C

    # Computational time
    #print('Computational time for Kernel K-Means with EM approach (sec): ',time.time() - start);

    return C_kmeans, En_kmeans





######################################
# Usages: 
#   [C_kmeans,En_kmeans] = compute_kernel_kmeans_spectral(K,Ker,Theta,nb_trials)
#
# Notations:
#   n = nb_data
#
# Input variables:
#   K = nb_clusters
#   Ker = Kernel matrix. Size = n x n.
#   Theta = Weight for each data term. Size = n x 1.
#   nb_trials = Number of standard kmeans runs.
#
# Output variables:
#   C_kmeans = Computed kmeans clusters. Size = n x 1.
#   En_kmeans = Energy of the kmeans partition.
######################################

def compute_kernel_kmeans_spectral(nc,Ker,Theta,nb_trials):


    start = time.time()
    n = Ker.shape[0]
    Theta = np.diag(Theta) # Weight for each data
    Ones = np.ones((1,n))

    # Eigenvalue Decomposition (EVD)
    A = (pow(Theta,0.5)).dot( Ker.dot(pow(Theta,0.5)) )
    lamb, U = scipy.sparse.linalg.eigsh(A, k=nc, which='LM') 
    U = U[:,::-1] # largest = index 0
    lamb = lamb[::-1]

    # Pre-process embedding coordinates Y
    Y = U - np.mean(U,axis=0) # zero-centered data
    Y = ( Y.T / np.sqrt(np.sum(Y**2,axis=1)+1e-10) ).T # L2 normalization of rows of Y

    # Run Standard/Linear K-Means on embedding coordinates Y
    Ker = construct_kernel(Y,'linear')
    C, En = compute_kernel_kmeans_EM(nc,Ker,Theta,10)

    # Outputs
    C_kmeans = C
    En_kmeans = En

    # Computational time
    #print('Computational time for Kernel K-Means with Spectral approach (sec):',time.time() - start);

    return C_kmeans, En_kmeans







######################################################
# Incremental Reseeding Algorithm for Clustering
# Xavier Bresson, Huiyi Hu, Thomas Laurent, Arthur Szlam, James von Brecht
# arXiv:1406.3837
######################################################

def plant_seeds(C, n, R, seeds_per_class):

    # Init
    F = np.zeros((n,R))
    Bal = np.zeros(R)

    # Construct indicators of seeds
    for k in range(R):
        idx_in_class_k = [ i for i,item in enumerate(C) if C[i]==k ]
        size_class_k = len(idx_in_class_k)
        if size_class_k > n:
            size_class_k = int(n)
        seed_idx_k = np.random.permutation(range(size_class_k))
        seed_idx_k = seed_idx_k[0:int(seeds_per_class)]
        seed_idx_k = [ idx_in_class_k[item] for i,item in enumerate(seed_idx_k) ]
        F[seed_idx_k,k] = 1

    return F

def compute_pcut(W,Cgt,R,speed=5.0,max_nb_iter=50,display=False):

    # Compute Random Walk Transition Matrix
    D = scipy.sparse.csc_matrix.sum(W,axis=1) # sum along rows
    Dinv = 1/D
    Dinv = Dinv.squeeze()
    n = D.shape[0]
    invD = scipy.sparse.spdiags(Dinv,0,n,n)
    RW = W * invD

    # Incremental Clustering
    speed = float(speed)
    maxiter = int(max_nb_iter)
    delta_seeds = np.array(speed*10**-4*n/R)
    nb_of_seeds = np.array(1) # initial nb of seeds 
    seeds_per_class = nb_of_seeds
    eps = np.array(10**-6)
    seed_idx = np.array(np.random.permutation(range(n)),dtype=np.int64)
    seed_idx = seed_idx[0:nb_of_seeds*R]
    C = np.random.randint(R,size=n)
    purity_tab = np.zeros(maxiter)

    # Loop
    for iter in range(maxiter):

        # Plant Seeds
        F = plant_seeds(C,n,R,seeds_per_class)

        # Random Walk
        k = 0
        while (np.min(F)<eps) & (k<50):
            F = RW.dot(F)
            k += 1

        # Thresholding
        C = np.argmax(F,axis=1)

        # Update seeds
        nb_of_seeds = nb_of_seeds + delta_seeds
        seeds_per_class = round(nb_of_seeds)

        # Purity
        if display & (not iter%10):
            print('iter=', iter, ', accuracy=', compute_purity(C,Cgt,R))
        
    return C, compute_purity(C,Cgt,R)






######################################
# 
# Usages: 
# (i)    Ker = construct_kernel(X,'linear');              # Ker = <Xi,Xj>
# (ii)   Ker = construct_kernel(X,'polynomial',[a,b,c]);  # Ker = ( a* <Xi,Xj> + b )^c
# (iii)  Ker = construct_kernel(X,'gaussian');            # Ker = exp( -|Xi-Xj|_2^2 / a ), automatic a
# (iv)   Ker = construct_kernel(X,'gaussian',a);          # Ker = exp( -|Xi-Xj|_2^2 / a ), given a
# (v)    Ker = construct_kernel(X,'sigmoid',[a,b]);       # Ker = tanh( a* <Xi,Xj> + b )
# (vi)   Ker = construct_kernel(X,'kNN_gaussian',k);
# (vii)  Ker = construct_kernel(X,'kNN_cosine',k);
# (viii) Ker = construct_kernel(X,'kNN_cosine_binary',k);
#
# Notations:
#   n = nb_data
#   d = data_dimensionality
#   k = nb_nearest_neighbors
#
# Input variables:
#   X = Data matrix. Size = n x d.
#   type_kernel = Type of kernels: 'linear', 'polynomial', 'gaussian', 'sigmoid', 'kNN'.
#
# Output variables:
#   Ker = Kernel matrix. Size = n x n.
######################################

def construct_kernel(X,type_kernel,parameters=None):
#def construct_kernel(X,type_kernel):

    start = time.time()
    n = X.shape[0]


    ####################################################
    # Linear Kernel = <Xi,Xj>
    # Ker = <Xi,Xj>
    ####################################################

    if type_kernel=='linear':
    
        print('Construct Linear Kernel')
    
        if isinstance(X, np.ndarray)==False:
            X = X.toarray()
    
        # Compute D = <Xi,Xj> for all pairwise data points size(D) =  n x n
        
        Ker = X.dot(X.T)
        

    ####################################################
    # Gaussian Kernel
    # Ker = exp( -|Xi-Xj|_2^2 / a ), automatic a
    ####################################################

    if type_kernel=='gaussian':
    
        print('Construct Gaussian Kernel')

        # Compute L2/Euclidean distance between all pairs of points
        if isinstance(X, np.ndarray)==False:
            X = X.toarray()        
        #Xzc = X - np.mean(X,axis=0) # zero-centered data
        D = sklearn.metrics.pairwise.pairwise_distances(X, metric='euclidean', n_jobs=1)
        Ddist = 1.0* D

        # Parameters
        if parameters==None:
            
            k = round(n/100.)
            #print('k=',k)

            # Sort Distance matrix
            idx = np.argsort(D)[:,:k] # indices of k nearest neighbors
            D.sort() # sort D from smallest to largest values
            D = D[:,:k]

            # Compute Ker matrix
            sigma2 = np.mean(D[:,-1])**2 # kernel scale
            #print('sigma2c',sigma2)

        else:

            sigma = parameters
            sigma2 = sigma**2
            #print('sigma2c',sigma2)

        Ker = np.exp(- Ddist**2 / sigma2)

        # Make Ker is symmetric
        bigger = Ker.T > Ker
        Ker = Ker - Ker.dot(bigger) + Ker.T.dot(bigger)



    ####################################################
    # Polynomial Kernel
    # Ker = ( a* <Xi,Xj> + b )^c
    ####################################################

    if type_kernel=='polynomial':
    
        print('Construct Polynomial Kernel')

        # Parameters
        a = parameters[0]
        b = parameters[1]
        c = parameters[2]

        # Compute Cosine distance between all pairs of points
        if isinstance(X, np.ndarray)==False:
            X = X.toarray()
        
        Xzc = X - np.mean(X,axis=0) # zero-centered data
        Xl2proj = ( Xzc.T / np.sqrt(np.sum(Xzc**2,axis=1)+1e-10) ).T # Projection on the sphere, i.e. ||x_i||_2 = 1
        D = Xl2proj.dot(Xl2proj.T)

        # Polynomial Kernel
        Ker = pow( np.abs( a* D + b ) , c )
    
        # Make Ker is symmetric
        bigger = Ker.T > Ker
        Ker = Ker - Ker.dot(bigger) + Ker.T.dot(bigger)



    ####################################################
    # kNN Gaussian Kernel
    ####################################################

    if type_kernel=='kNN_gaussian':
    
        print('Construct kNN Gaussian Kernel')
    
        # Parameters
        k = parameters

        # Compute L2/Euclidean distance between all pairs of points
        if isinstance(X, np.ndarray)==False:
            X = X.toarray()        
        #Xzc = X - np.mean(X,axis=0) # zero-centered data
        D = sklearn.metrics.pairwise.pairwise_distances(X, metric='euclidean', n_jobs=1)

        # Sort Distance matrix
        idx = np.argsort(D)[:,:k] # indices of k nearest neighbors
        D.sort() # sort D from smallest to largest values
        D = D[:,:k]

        # Compute Ker matrix
        sigma2 = np.mean(D[:,-1])**2 # graph scale
        Ker_val = np.exp(- D**2 / sigma2)

        # Make Ker sparse
        n = X.shape[0]
        row = np.arange(0, n).repeat(k)
        col = idx.reshape(n*k)
        data = Ker_val.reshape(n*k)
        Ker = scipy.sparse.csr_matrix((data, (row, col)), shape=(n, n))

        # Make Ker is symmetric
        bigger = Ker.T > Ker
        Ker = Ker - Ker.multiply(bigger) + Ker.T.multiply(bigger)
        Ker = Ker.todense()



    ####################################################
    # kNN Cosine Kernel
    ####################################################

    if type_kernel=='kNN_cosine':
    
        print('Construct kNN Cosine Kernel')
    
        # Parameters
        k = parameters
        
        # Compute Cosine distance between all pairs of points
        if isinstance(X, np.ndarray)==False:
            X = X.toarray()
        
        Xzc = X - np.mean(X,axis=0) # zero-centered data
        Xl2proj = ( Xzc.T / np.sqrt(np.sum(Xzc**2,axis=1)+1e-10) ).T # Projection on the sphere, i.e. ||x_i||_2 = 1
        D = Xl2proj.dot(Xl2proj.T)

        # Sort D according in descending order     
        idx = np.argsort(D)[:,::-1][:,:k] # indices of k nearest neighbors
        D.sort(axis=1)
        D[:] = D[:,::-1]

        # Cosine distance
        Dcos = np.abs(np.arccos(D-1e-10))
        D = Dcos
        D = D[:,:k]

        # Compute Ker matrix
        sigma2 = np.mean(D[:,-1])**2 # graph scale
        Ker_val = np.exp(- D**2 / sigma2)

        # Make Ker sparse
        n = X.shape[0]
        row = np.arange(0, n).repeat(k)
        col = idx.reshape(n*k)
        data = Ker_val.reshape(n*k)
        Ker = scipy.sparse.csr_matrix((data, (row, col)), shape=(n, n))

        # Make W is symmetric
        bigger = Ker.T > Ker
        Ker = Ker - Ker.multiply(bigger) + Ker.T.multiply(bigger)
        Ker = Ker.todense()


    ####################################################
    # kNN Cosine Binary Kernel
    ####################################################

    if type_kernel=='kNN_cosine_binary':
    
        print('Construct kNN Cosine Binary Kernel')
    
        # Parameters
        k = parameters

        # Compute Cosine distance between all pairs of points
        if isinstance(X, np.ndarray)==False:
            X = X.toarray()
        
        Xzc = X - np.mean(X,axis=0) # zero-centered data
        Xl2proj = ( Xzc.T / np.sqrt(np.sum(Xzc**2,axis=1)+1e-10) ).T # Projection on the sphere, i.e. ||x_i||_2 = 1
        D = Xl2proj.dot(Xl2proj.T)

        # Sort D according in descending order     
        idx = np.argsort(D)[:,::-1][:,:k] # indices of k nearest neighbors
        D.sort(axis=1)
        D[:] = D[:,::-1]

        # Cosine distance
        Dcos = np.abs(np.arccos(D-1e-10))
        D = Dcos
        D = D[:,:k]

        # Compute Ker matrix
        sigma2 = np.mean(D[:,-1])**2 # graph scale
        Ker_val = np.exp(- D**2 / sigma2)

        # Make Ker sparse
        n = X.shape[0]
        row = np.arange(0, n).repeat(k)
        col = idx.reshape(n*k)
        data = np.ones([n*k])
        Ker = scipy.sparse.csr_matrix((data, (row, col)), shape=(n, n))

        # Make W is symmetric
        bigger = Ker.T > Ker
        Ker = Ker - Ker.multiply(bigger) + Ker.T.multiply(bigger)
        Ker = Ker.todense()






    #print('Computational time for Kernel construction (sec):',time.time() - start)
    
    return Ker







######################################
# 
# Usage: 
#   [alpha,f_test,C_test,acc_test] = compute_SVM(Xtrain,l,Xtest,Cgt_test,gammaG,type_graph,parameters_graph,gamma,type_kernel,parameters)
#
# Notations:
#   n_train = nb of training data
#   n_test = nb of test data
#   d = data_dimensionality
#   k = nb_nearest_neighbors
#
# Input variables:
#   Xtrain = Training data matrix. Size = n_train x 1.
#   l = Label vector. Size = n_train x 1. Values = +/-1.
#   Xtest = Test data matrix. Size = n_test x 1.
#   Cgt_test = Ground truth clusters of the test dataset. Size = n_test x 1.
#   gammaG = Graph regularization paramater.
#   type_graph = Type of graphs: 'no_graph', 'euclidean', 'euclidean_zelnik_perona', 'cosine', 'cosine_binary'
#   gamma = SVM regularization paramater.
#   type_kernel = Type of kernels: 'linear', 'polynomial', 'gaussian', 'sigmoid', 'kNN'.
#   parameters = parameters used for Kernel construction
#
# Output variables:
#   alpha = Coefficient vector of classification function. Size = n_train x 1.
#   f_test = Continuous decision function. Size = n_test x 1.
#   C_test = Binary decision function. Size = n_test x 1.
#   acc_test = Accuracy of SVM classification on test dataset.
######################################

def compute_SVM(Xtrain,Cgt_train,l_train,type_svm,param_svm=None,Xtest=None,Cgt_test=None,plot=None,type_graph=None,param_graph=None):


    # Parameters
    n = Xtrain.shape[0]
    nc = len(np.unique(Cgt_train))
    gamma = param_svm[0]
    if not (plot==False):
        fig_id = plot[0]
        size_vertex_plot = plot[1]

    # Type of SVM technique
    if type_graph==None:
        if type_svm=='soft_linear':
            print('Run Linear SVM')
        
        if type_svm=='gaussian_kernel':
            print('Run Kernel SVM')

    if not (type_graph==None):
        if type_svm=='gaussian_kernel':
            print('Run Graph SVM with Gaussian Kernel')

    # Compute kernel
    if type_svm=='soft_linear':
        Ker = construct_kernel(Xtrain,'linear')

    if type_svm=='gaussian_kernel':
        if len(param_svm)<2: Ker = construct_kernel(Xtrain,'gaussian') # TODO
        sigma = param_svm[1]
        Ker = construct_kernel(Xtrain,'gaussian',sigma) 

    # Compute graph
    if not (type_graph==None):
        kNN = param_graph[0]
        gammaG = param_graph[1]
        W = construct_knn_graph(Xtrain,kNN,type_graph)
        Lap = graph_laplacian(W).todense()

    # Compute test kernel for classification error
    if type_svm=='soft_linear':
        KXtest = Xtrain.dot(Xtest.T)

    if type_svm=='gaussian_kernel':
        Ddist = sklearn.metrics.pairwise.pairwise_distances(Xtrain,Xtest, metric='euclidean', n_jobs=1)
        sigma2 = sigma**2
        KXtest = np.exp(- Ddist**2 / sigma2)

    # Compute Indicator function of labels
    H = np.zeros([n])
    H[np.abs(l_train)>0.0] = 1
    H = np.diag(H)

    # Compute L, Q
    L = np.diag(l_train)
    l = l_train
    T = np.eye(n) 
    if not (type_graph==None): 
        T += gammaG* Lap.dot(Ker) 

    Tinv = np.linalg.inv(T)
    Q = L.dot(H.dot(Ker.dot(Tinv.dot(H.dot(L)))))

    # Time steps
    sigma = 1./ np.linalg.norm(L,2)
    tau = 1./ np.linalg.norm(Q,2)

    # For conjuguate gradient
    Acg = tau* Q + np.eye(n)

    # Initialization
    x = np.zeros([n])
    y = 0.0
    xold = x

    # Plot
    if not (plot==None):
        fig = plt.figure(fig_id)
        fig.canvas.draw()
        plt.show(block=False)

    # Loop
    k = 0
    diffx = 1e6
    while (diffx>1e-3) & (k<1000):
            
        # Update iteration
        k += 1
            
        # Update x
        # Approximate solution with conjuguate gradient
        b = x + tau - tau* l* y
        x,_ = scipy.sparse.linalg.cg(Acg, b, x0=x, tol=1e-3, maxiter=50)   
        x[x<0.0] = 0 # Projection on [0,gamma]
        x[x>gamma] = gamma 

        # Update y
        y = y + sigma* l.T.dot(x)
            
        # Stopping condition
        diffx = np.linalg.norm(x-xold)
        xold = x
            
        # Plot
        if not(k%50):
                   
            # Lagrangian multipliers
            alpha = x

            # a vector
            a = Tinv.dot(H.dot(L.dot(alpha)))

            # b offset
            idx_unlabeled_data = np.where( np.abs(l)<1./2 )
            alpha_labels = alpha; alpha_labels[idx_unlabeled_data] = 0
            idx = np.where( np.abs(alpha_labels)>0.25* np.max(np.abs(alpha_labels)) )
            Isv = np.zeros([n]); Isv[idx] = 1 # Indicator function of Support Vectors
            nb_sv = len(Isv.nonzero()[0])
            if nb_sv > 1:
                b = (Isv.T).dot( l - Ker.dot(a)  )/ nb_sv
            else:
                b = 0

            # Continuous decision function
            f_test = a.T.dot(KXtest) + b 

            # Binary decision function
            C_test = np.array( 1./2* ( 1+ np.sign(f_test) ) , dtype=np.int64) # decision function in {0,1}
            acc_test = compute_purity(C_test,Cgt_test,nc)

            # Plot
            if not (plot==None):
                plt.subplot(121)
                plt.scatter(Xtest[:,0], Xtest[:,1], s=size_vertex_plot*np.ones(n), c=f_test)
                plt.title('Continuous decision function')
                plt.subplot(122)
                plt.scatter(Xtest[:,0], Xtest[:,1], s=size_vertex_plot*np.ones(n), c=C_test)
                plt.title('Decision function, iter=' + str(k) + ', acc=' + str(round(acc_test,3)) )
                plt.tight_layout()
                plt.show()
                fig.canvas.draw()
                time.sleep(0.01)
            

    # Final operations
    # Lagrangian multipliers
    alpha = x
    # a vector
    a = Tinv.dot(H.dot(L.dot(alpha)))
    
    # b offset
    idx_unlabeled_data = np.where( np.abs(l)<1./2 )
    alpha_labels = alpha; alpha_labels[idx_unlabeled_data] = 0
    idx = np.where( np.abs(alpha_labels)>0.25* np.max(np.abs(alpha_labels)) )
    Isv = np.zeros([n]); Isv[idx] = 1 # Indicator function of Support Vectors
    nb_sv = len(Isv.nonzero()[0])
    if nb_sv > 1:
        b = (Isv.T).dot( l - Ker.dot(a)  )/ nb_sv
    else:
        b = 0
    #print('b=',b)

    # Continuous decision function
    f_test = a.T.dot(KXtest) + b
    # Binary decision function
    C_test = np.array( 1./2* ( 1+ np.sign(f_test) ) , dtype=np.int64) # decision function in {0,1}
    acc_test = compute_purity(C_test,Cgt_test,nc)

    # Plot
    if not (plot==None):
        plt.subplot(121)
        plt.scatter(Xtest[:,0], Xtest[:,1], s=size_vertex_plot*np.ones(n), c=f_test)
        plt.title('Continuous decision function')
        plt.subplot(122)
        plt.scatter(Xtest[:,0], Xtest[:,1], s=size_vertex_plot*np.ones(n), c=C_test)
        plt.title('Decision function, iter=' + str(k) + ', acc=' + str(round(acc_test,3)) )
        plt.tight_layout()
        plt.show()
        fig.canvas.draw()
        time.sleep(0.01)


    return alpha, f_test, C_test, acc_test





######################################
# Soft shrinkage operator
######################################

def shrink(x,mu):

    s = np.sqrt(x**2)
    ss = s - mu
    ss = ss* ( ss>0.0 )
    s = s + ( s<mu )
    ss = ss/ s
    res = ss* x

    return res




