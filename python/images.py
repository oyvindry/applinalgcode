import matplotlib.pyplot as plt
from numpy import *

# Generic image functionality

def imread( filename, format = ''):
    img = plt.imread(filename)
    if img.dtype == dtype('float32'):
        img *= 255
        img = img.astype(uint8)
    return img

# Assume that data is uint8
def imwrite(data, filename, format = ''):
    if ndim(data) == 3:
        plt.imsave(filename, data, origin='upper')
    else:
        data = data/255.
        data = data.astype(float32)
        plt.imsave(filename, data, origin='upper', cmap = 'gray')

def imshow(data):
    if ndim(data) == 3:
        plt.imshow(data)
    else :
        data = data/255.
        data = data.astype(float32)
        plt.imshow(data, cmap = 'gray')
    plt.show()
    


# Functions for image manipulation
        
def mapto01(X):
    minval, maxval = X.min(), X.max()
    X -= minval
    X /= (maxval-minval)
    
def contrastadjust(X,epsilon):
    """
    Assumes that the values are in [0,255]
    """
    X /= 255.
    X += epsilon
    log(X, X) 
    X -= log(epsilon)
    X /= (log(1+epsilon)-log(epsilon))
    X *= 255
  
def contrastadjust0(X,n):
    """
    Assumes that the values are in [0,255]
    """
    X /= 255.
    X -= 1/2.
    X *= n
    arctan(X, X)
    X /= (2*arctan(n/2.)) 
    X += 1/2.0
    X *= 255 # Maps the values back to [0,255]
    
def combineimages(imgs):
    N = len(imgs[0])
    M = len(imgs)
    sz = shape(imgs[0][0])
    ind = sz[0]/20 + 1
    newsz = [i for i in sz]
    newsz[0] *= M
    newsz[0] += (M-1)*ind
    newsz[1] *= N
    newsz[1] += (N-1)*ind
    newimg = 255*ones(newsz)
    for m in range(M):
        for n in range(N):
            newimg[(m*(sz[0]+ind)):((m+1)*sz[0]+m*ind), (n*(sz[1]+ind)):((n+1)*sz[1]+n*ind)] = imgs[m][n][:]
    return newimg

def CreateExcerpt():
    img = double(imread('images/lena.png','png'))
    return img[128:384,128:384,:]