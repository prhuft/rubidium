"""
General helper functions for simulations

Preston Huft
"""
import csv
import numpy as np
from random import random as rand
import matplotlib.pyplot as plt

## reading to and writing from files

def soln_to_csv(fname, data, labels):
    """
    Args:
        fname: myfile.csv
        data: a list or array of equal length lists or arrays of data
        labels: a label describing each list in data
    
        e.g., 
        data = [array([1,2,3]), array([2,4,6], array([1,4,9]]
        labels = ['x', '2x', 'x^2']
    """
    
    if type(data) == np.ndarray:
    # listify to avoid malforming strings
        data = [d for d in data]
    
    with open(fname, 'w', newline='') as f:
        writer = csv.writer(f, delimiter=',')

        for d,l in zip(data, labels):
            writer.writerow([l] + list(d))
            
    print(f"wrote data to {fname}")

def soln_from_csv(fname):
    """
    Args:
        fname: myfile.csv
    Returns:
        data: an array of equal length arrays of data
        labels: a label describing each array in data
    
        e.g. 
            data = [array([1,2,3]), array([2,4,6], array([1,4,9]]
            labels = ['x', '2x', 'x^2']
    """
    labels = []
    data = []
    
    with open(fname, 'r', newline='') as f:
        reader = csv.reader(f, delimiter=',')

        for row in reader:
            labels.append(row[0])
            try:
                data.append(np.array([complex(x) for x in row[1:]]))
            except (ValueError,TypeError) as e:
                print(e)
                print(f"problematic row: {row}")
                break
            
    return data, labels
    
## data checks and sanitization

def is_complexable(x):
    try:
        complex(x)
    except ValueError:
        return False
    else:
        return True

def is_numeric(x):
    test = lambda x: np.array([is_complexable(a) for a in x]).prod() 
    try:
        l = len(x)
        x = np.array(x).flatten()
        return bool(test(x))
    except TypeError:
        return bool(test([x]))
        
def is_scalar_2D(xarr):
    """check if 2D array passed in has non-list like elements"""
    
    dimx,dimy = np.array(xarr).shape
    for i in range(dimx):
        for j in range(dimy):
            try:
                len(xarr[i,j])
                print(f"non-scalar like at elem ({i},{j}): {xarr[i,j]}")
                return False
            except Exception as e:
                pass
    return True

def sample_dist(fx,fmax,domain,n,integral=False,showplot=False):
    """
    Get random samples from a probability distribution function (PDF)
    
    Args:
        fx: callable, the 1D continuous PDF to be sampled, i.e. f(x)
        fmax: the maximum of f(x) on a domain [x1,x2]
        n: number of samples to take
        domain: list, domain over which to sample f(x), i.e. [x1,x2]
        integral: set to True if fx is an integer distribution, i.e. only valid for x an integer (e.g. a Poisson distribution)
        showplot: bool, an optional graphical check of this function
    Return:
        x_dist: np array (float), the samples we took, shape (n,).
    """
    x1,x2 = domain
    y_dist = np.empty(n) 
    f_dist = np.empty(n) 
    x_dist = np.empty(n) # this is the distribution we want
    j = 0 # dist index
    
    if not integral:
        xsamp = lambda: (x2-x1)*rand()+x1
    else:
        xsamp = lambda: int((x2-x1)*rand()+x1 + 0.5)
    
    while j < n:
        x = xsamp() # rand val on domain of f(x)
        f = fx(x)
        y = rand()*fmax # rand val on range of f(x)
        if y <= f:
            y_dist[j]=y
            f_dist[j]=f
            x_dist[j]=x
            j+=1

    # plot distribution as a check:
    if showplot is not False:
        plt.scatter(x_dist,y_dist,c='red',s=10)
        plt.scatter(x_dist,f_dist,c='blue',s=10)
        plt.show()

    return x_dist


    
