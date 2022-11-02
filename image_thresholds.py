import numpy as np
import matplotlib.pyplot as plt
import pickle
from matplotlib import patches

def return_com_array(array: np.array):
    """short helper function to return the center of mass of an array in pixles

    Args:
        array (np.array): array of which to return the center of mass

    Returns:
        tuple: tuple(y, x) center of mass
    """
    height, width = array.shape[0], array.shape[1]
    X, Y = np.meshgrid(range(height), range(width))
    total = np.sum(array)
    y_mean, x_mean = np.sum(Y*array)/total, np.sum(X*array)/total
    return y_mean, x_mean

def return_thresholds(array, thresholds=None):
    """For a given list of thresholds and an intensity map,
    a list of (x, y) are returned such that the intensity contained by the centered rectangle with height 2*x and width 2*y
    surpasses the respective threshold.

    Args:
        array (np.array): array of intensities
        thresholds (list, optional): List of thresholds default is ten equidistant steps from 0.1 to 1. Defaults to None.

    Returns:
        list: list of tuple (x, y) for array-centered rectangle with heigth 2*x and width 2*y
    """
    if thresholds is None:
        thresholds = np.linspace(0.1, 0.9, 9)
    threshold_xy = []
    shape = array.shape
    all = np.sum(array)
    midpoint = int(shape[0]/2), int(shape[1]/2)
    xmin = 0
    y_divided_x = round(shape[1]/shape[0]) # calc the y steps from x
    print(y_divided_x)
    for ind, threshold in enumerate(thresholds):
        for x in range(xmin, midpoint[0]):
            y = round(y_divided_x*x)
            try:
                data = array[midpoint[0]-x: midpoint[0]+x+1, midpoint[1]-y: midpoint[1]+y+1]
            except IndexError:
                break
            if np.sum(data)/all >= threshold:
                threshold_xy.append((x, y))
                xmin = x
                break
    return threshold_xy

def plot_thresholds(array, fig=None, ax=None, thresholds=None, extent=None, circle=None):
    """plots the calculated thresholds in a rectangular map

    Args:
        array (np.array): data array for wich the thresholds are calculated
        fig (plt.fig, optional): Figure in which to plot if none is given a new figure is created. Defaults to None.
        ax (plt.ax, optional): ax in which to plot if none is give a new one is created. Defaults to None.
        thresholds (list, optional): list of thresholds defaults to (0.1, 0.2, ..., 0.9). Defaults to None.
        circle (center, radius): if a center and radius are given

    Returns:
        (fig, ax): tuple of figure and ax
    """
    if thresholds is None:
        thresholds = np.linspace(0.1, 0.9, 9)
    threshold_xy = return_thresholds(array, thresholds)
    midpoint = int(array.shape[0]/2), int(array.shape[1]/2)
    if extent==None:
        extent=(0, array.shape[0], 0, array.shape[1])
    height = extent[1]-extent[0]
    width = extent[3]-extent[2]

    def px2datapoints(px, py):
        x = (px/array.shape[0])*(height)+ extent[0]
        y = (py/array.shape[1])*(width)+ extent[2]
        return x, y

    if fig==None or ax==None:
        fig, ax = plt.subplots()
        im = ax.imshow(array[::-1], interpolation='none', extent=extent)

    for ind, (x, y) in enumerate(threshold_xy):
        lowerleft_x, lowerleft_y = px2datapoints(midpoint[0]-x, midpoint[1]-y)
        rect = patches.Rectangle((lowerleft_x, lowerleft_y), 2*x*height/array.shape[0], 2*y*width/array.shape[1], alpha=1, fc='none', edgecolor='red')
        ax.add_patch(rect)
        ax.text(lowerleft_x+0.2, lowerleft_y+0.2, round(thresholds[ind], 2), color='red')
    #fig.colorbar(im)
    return fig, ax


def return_thresholds_circle(array: np.array, xmid_px, ymid_px, extent: tuple=None, thresholds: list=None) -> np.array:
    """finds and returns the radii of circles crossing certain thresholds

    Args:
        array (np.array): array of data to find thresholds for
        extent (tuple): extent of the array. Defaults to (0,xpix,0,ypix)
        thresholds (list): list of thresholds. Defaults to np.linspace(0.1, 0.9, 9)

    Returns:
        np.array: list(ind, threshold, r[px]) list for all thresholds
    """
    if thresholds is None:
        thresholds = np.linspace(0.1, 0.9, 9)
    threshold_radii = []
    total = np.sum(array)
    radius = 0
    for ind, threshold in enumerate(thresholds):
        for r in range(radius, max([array.shape[0], array.shape[1]])):
            try:
                X, Y = np.meshgrid(range(array.shape[0]), range(array.shape[1]))
                mask = ((X-xmid_px)**2 + (Y-ymid_px)**2 < r**2)
                if np.sum(array[mask])/total > threshold:
                    threshold_radii.append((ind, threshold, r))
                    radius = r
                    break
            except IndexError:
                break
        else:
            print('no suc', threshold)
    return threshold_radii

def plot_thresholds_circle(array: np.array, xmid_px=None, ymid_px=None, extent: tuple=None, thresholds: list=None, figax=None):
    if extent is None:
        extent = (-0.5, array.shape[0]-0.5, -0.5, array.shape[1]-0.5)
    if xmid_px is None:
        xmid_px = array.shape[0]//2
    if ymid_px is None:
        ymid_px = array.shape[1]//2
    midpointx, midpointy = extent[0] + (xmid_px/(array.shape[0]-1))*(extent[1]-extent[0]), extent[2] + (ymid_px/(array.shape[1]-1))*(extent[3]-extent[2])
    print(midpointx, midpointy)
    if figax is None:
        fig, ax = plt.subplots(1)
        im = ax.imshow(array[::-1], extent=extent, interpolation='none')
    else:
        fig, ax = figax
    # transform the coordinates to indices
    thresholds=return_thresholds_circle(array,  xmid_px=xmid_px, ymid_px=ymid_px, extent=extent, thresholds=thresholds)
    for ind, (i, threshold, r_px) in enumerate(thresholds):
        r_ext = r_px*(extent[1]-extent[0])/array.shape[0]
        print(r_ext, threshold)
        circ = patches.Circle((midpointx, midpointy), r_ext, fc='none', edgecolor='red')
        ax.add_patch(circ)
        ax.text(xmid_px/(array.shape[0]-1)+r_px/(array.shape[0]-1), ymid_px/(array.shape[1]-1),'$I_p$ = {:.2}'.format(threshold), color='red', transform=ax.transAxes)
    return fig, ax


def radius_vs_content(array: np.array, mid_x: float=None, mid_y: float=None, radii=None):
    """takes an array and returns the radius of the cirle vs the incorporated content

    Args:
        array (np.array): array with data
        mid_x (float, optional): x-value of the midpoint. Defaults to None.
        mid_y (float, optional): y-value of the midpoint. Defaults to None.

    Returns:
        (np.array): [(rad, int)]
    """
    height = array.shape[0]
    width = array.shape[1]
    total = np.sum(array)
    rad_int = []
    if mid_x is None:
        mid_x = (height)/2
    if mid_y is None:
        mid_y = (width)/2
    if radii is None:
        radii = range(0, max([height, width]))
    X, Y = np.meshgrid(range(height), range(width))
    for radius in radii:
        #fig, ax = plt.subplots(1)

        circle = ((X-mid_x)**2 + (Y-mid_y)**2 < radius**2)

        #ax.imshow(array, interpolation ='none')

        #print(np.sum(array[circle]))
        fraction = np.sum(array[circle])/total
        rad_int.append((radius, fraction))
    return rad_int
#test_array = np.ones((100, 100))
#print(radius_vs_content(test_array, radii = [1,5, 10, 15], mid_x=60, mid_y=49.5))

#with open('data/data_finite_g_length160finite_beamspot.txt', 'rb') as myFile:
#    data_dict = pickle.load(myFile)
#with open('data/meta_data_finite_g_length160finite_beamspot.txt', 'rb') as myFile:
#    meta_dict = pickle.load(myFile)
##X, Y = np.meshgrid(range(100), range(100))
##circ = ((X-2)**2 + (Y-34)**2 <= 30**2)
##test_array[circ] += 1
#
#extent = (-21.8/2, 21.8/2, -21.8/2, 21.8/2)
#print(data_dict.keys())
##test_array = data_dict['f_psd.dat'][:100]
##X, Y = np.meshgrid(range(100), range(100))
##circ = ((X-40)**2 + (Y-12)**2 <= 30**2)
##test_array[circ] *= 20
##
##circ = ((X-1)**2 + (Y-12)**2 <= 5**2)
##test_array[circ] *= 20000
#fig, ax = plt.subplots(1)
#im = ax.imshow(test_array[::], interpolation='none')
#midy, midx = return_com_array(test_array)
##midy, midx = midy/99*21.8-21.8/2, midx/99*21.8-21.8/2
#print(midy, midx)
#plot_thresholds_circle(test_array, xmid_px=midx, ymid_px=midy, extent=extent, thresholds=[0.30, 0.5, 0.70])
#print('center of mass', return_com_array(test_array))
##fig, ax = plot_thresholds_circle(test_array)
##print(return_thresholds_circle(test_array))
##plt.show()
##print('done')
#plt.show()
def calc_efficiency(array, height, axis, p=0, **kwargs):
    """returns the intensity integrated from 
        HENLO
    Args:
        array (np.array): array with the data
        height (float): height of the roi area in pixels
        axis (int): axis along which to integrate the data 0 if vertical 1 if horizontal

    Returns:
        float: intensity within the roi
    """
    #first integrated the efficiency along one axis
    integrated = np.sum(array, axis=axis)
    
    length = len(integrated)
    #determine the center of mass
    center_of_mass = 1/np.sum(integrated)*sum([k*ind for ind, k in enumerate(integrated)])
    center_of_mass = np.argmax(integrated)
    #better to run through all of the possible windows and check the max
    max_int = 0
    max_center = 0
    for center_of_mass in range(len(array)):
        mask = (np.array(range(length))-center_of_mass)**2 <= (height/2)**2
        if np.sum(integrated[mask]) > max_int:
            max_center = center_of_mass
            max_int = np.sum(integrated[mask])
    #only for testing delete after
    #max_center=500#WARNING DELETE if moving window is wanted
    print('max center = {}'.format(max_center))
    mask = (np.array(range(length))-max_center)**2 <= (height/2)**2

    if p:
        fig, ax = plt.subplots(1)
        ax.plot(integrated)
        ax.plot(mask)
        ax.set_xticklabels((float(k)-500)*0.5 for k in list(ax.get_xticks()))
        ax.set_xlabel('y (mm)')
        ax.set_ylabel('intensit (arb. u.)')
        ax.set_yscale('log')
        ax.tick_params(which='both', direction='in')
        title= ''
        for key in kwargs:
            title += "{} = {} ".format(key, kwargs[key])
        ax.set_title(title)
    return np.sum(integrated[mask])