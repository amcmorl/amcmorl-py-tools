import numpy as np

def pca(M):
    "Perform PCA on M, return eigenvectors and eigenvalues, sorted."
    M = np.asarray(M)
    
    T, N = M.shape
    P = M - M.mean(axis=0)[np.newaxis,...]

    # if there are less rows T than columns N, use
    # snapshot method
    if T < N:
        C = np.dot(P, np.transpose(P)) / float(T-1)
        evals, evecsC = np.linalg.eig(C)
        # HACK: make sure evals are all positive
        evals = np.where(evals < 0, 0, evals)
        evecs = 1./np.sqrt(evals) * np.dot(np.transpose(P),
                                           np.transpose(evecsC))
    else:
        # calculate covariance matrix
        K = 1./float(T-1) * np.dot(np.transpose(P), P)
        evals, evecs = np.linalg.eig(K)
    # sort the eigenvalues and eigenvectors, decending order
    order = (np.argsort(evals)[::-1])
    evecs = np.take(evecs, order, 1)
    evals = np.take(evals, order)
    return evals, np.transpose(evecs)
