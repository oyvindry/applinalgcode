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
        
def map_to_01(X):
    minval, maxval = X.min(), X.max()
    X -= minval
    X /= (maxval-minval)
# End map_to_01

def contrast_adjust(X,epsilon):
    # Assumes that the values are in [0,255]
    X /= 255.
    X += epsilon
    log(X, X) 
    X -= log(epsilon)
    X /= (log(1+epsilon)-log(epsilon))
    X *= 255
# End contrast_adjust