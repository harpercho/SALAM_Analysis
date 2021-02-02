import numpy as np
import matplotlib.pyplot as plt
import sys, os, glob, pickle


def load_halos_pickle(pickle_path):
    
    '''
    Returns pickle as dictionary
    '''
        
    data = pickle.load( open( pickle_path , "rb" ))
    
    output = dict([(str(k),np.zeros(len(data))) for k in data[0]])
    
    for i in range(len(data)):
        gal_dict = data[i]
        
        if gal_dict is not None:
            for key, value in output.items():
                try:
                    value[i] = gal_dict[str(key)] + 1e-12
                except:
                    pass
    return output


def density_estimate(x, y):
    
    '''
    Calculates the plot density given function m2(m1).
    Returns X, Y, Z
    '''
    
    from scipy import stats
    xmin = x.min()
    xmax = x.max()
    ymin = y.min()
    ymax = y.max()

    X, Y = np.mgrid[xmin:xmax:1000j, ymin:ymax:1000j]
    positions = np.vstack([X.ravel(),Y.ravel()])
    values = np.vstack([x, y])

    kernel = stats.gaussian_kde(values)

    Z = np.reshape(kernel(positions).T, X.shape)
    
    return X, Y, Z

def scatter_to_contour(ax, x, y):
    
    X, Y, Z = density_estimate(x, y)

    # Show denisty
    ax.imshow(np.rot90(Z), cmap='Blues', extent=[xmin, xmax, ymin, ymax])

    # Add contour lines
    plt.contour(X, Y, Z, alpha = 0.1)

    #ax.plot(x, y, 'k.', markersize=2,alpha=0.1)
    ax.set_xlim([xmin, xmax])
    ax.set_ylim([ymin, ymax])

    return ax

def get_plot_coords(n_plots, per_row):

    y = int(idx/per_row)
    x = idx % per_row

    return x, y
