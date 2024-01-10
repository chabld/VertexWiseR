"""Python function to smooth surface images
   Modified from BrainStat toolbox https://github.com/MICA-MNI/BrainStat/blob/master/brainstat/mesh/data.py 
   -removed third dimension
   -removed printing of messages (messes up the R console)
   -replace surf(dict, BSPolyData) input with ndarray edgelists that are loaded automatically depending on the number of columns in the CT data
   -remove edges in the medial wall
   -converted mesh units to mm
"""
import numpy as np
import requests
import io

def mesh_smooth(
    Y: np.ndarray,  FWHM: float
) -> np.ndarray:
    """Smooths surface data by repeatedly averaging over edges.

    Parameters
    ----------
    Y : numpy.ndarray
        Surface data of shape (n,v). v is the number of vertices,
        n is the number of observations.
    FWHM : float
       Gaussian smoothing filter in mesh units.

    Returns
    -------
    numpy.ndarray
        Smoothed surface data of shape (n,v).
    """
    if Y.shape[1]==20484:
        url = requests.get('https://raw.githubusercontent.com/CogBrainHealthLab/VertexWiseR/main/data/edgelistfs5.npy')
        edg=np.load(io.BytesIO(url.content))
        FWHM=FWHM/3.5
        
    elif Y.shape[1]==81924:
        url = requests.get('https://raw.githubusercontent.com/CogBrainHealthLab/VertexWiseR/main/data/edgelistfs6.npy')
        edg=np.load(io.BytesIO(url.content))
        FWHM=FWHM/2
   
    niter = int(np.ceil(pow(FWHM, 2) / (2 * np.log(2))))
    if isinstance(Y, np.ndarray):
        Y = np.array(Y, dtype="float")
        n, v = np.shape(Y)
        isnum = True
    
    agg_1 = np.bincount(edg[:, 0], minlength=(v + 1)) * 2
    agg_2 = np.bincount(edg[:, 1], minlength=(v + 1)) * 2
    Y1 = (agg_1 + agg_2)[1:]

    n10 = np.floor(n / 10)

    for i in range(0, n):
          if isnum:
        
              Ys = Y[i, :]

              for itera in range(1, niter + 1):
                  Yedg = Ys[edg[:, 0] - 1] + Ys[edg[:, 1] - 1]
                  agg_tmp1 = np.bincount(edg[:, 0], Yedg, (v + 1))[1:]
                  agg_tmp2 = np.bincount(edg[:, 1], Yedg, (v + 1))[1:]
                  with np.errstate(invalid="ignore"):
                      Ys = (agg_tmp1 + agg_tmp2) / Y1
                        
              Y[i, :] = Ys
                
    return Y


