from __future__ import print_function
import lightkurve as lk
import numpy as np
import matplotlib.pyplot as plt
import pandas
print('lk version:',lk.__version__)

# Fetching planet information
crossmatch_planets = pandas.read_csv('data/crossmatch_planets.csv')
crossmatch_candidates = pandas.read_csv('data/crossmatch_candidates.csv')
combined = pandas.concat([crossmatch_planets,crossmatch_candidates],sort=False)

def get_star_name(tic=None,epic=None):
    """Return the name of a host star based on its TIC ID."""
    if tic:
        index = combined.index[combined['TIC'] == tic].tolist()
    if epic:
        index = combined.index[combined['EPIC'] == epic].tolist()
    if index ==[]:
        return None
    name = combined['pl_hostname'][index].tolist()[0]
    return name

def get_tic_from_epic(epic):
    index = combined.index[combined['EPIC'] == epic].tolist()
    tic = combined['TIC'][index].tolist()[0]
    return tic

def get_epic_from_tic(tic):
    index = combined.index[combined['TIC'] == tic].tolist()
    epic = combined['EPIC'][index].tolist()[0]
    return epic

def get_tic_from_name(hostname):
    index = combined.index[combined['pl_hostname'] == hostname].tolist()
    tic = combined['TIC'][index].tolist()[0]
    return tic

def make_lc(tic=None,epic=None,hostname=None,mission='tess',verbose=False):
    """Return a LightCurve object from TESS or K2 planet.

    Keyword arguments:
    Need one of:
        tic -- TIC ID # of the host star
        epic -- EPIC ID # of the host star
        hostname -- host star name
    mission -- which mission to get the data from (default 'tess')
    verbose -- whether to print status updates (default False)
    """
    if tic == None:
        if epic:
            tic = get_tic_from_epic(epic)
        elif hostname:
            tic = get_tic_from_name(hostname)

    if verbose: print('Fetching TPF...')
    if mission.lower() == 'tess':
        tpf = lk.search_targetpixelfile(tic,mission='tess').download()
    elif mission.lower() == 'k2':
        epic = get_epic_from_tic(tic)
        tpf = lk.search_targetpixelfile(epic,mission='k2').download()
    else:
        raise ValueError('Mission must be \'tess\' or \'k2\'.')

    if verbose: print('De-trending...')
    if mission.lower() == 'tess':
        corrector = lk.PLDCorrector(tpf)
        lc = corrector.correct()
    elif mission.lower() == 'k2':
        lc = tpf.to_lightcurve(aperture_mask=tpf.pipeline_mask)
        lc.remove_outliers(sigma=6).flatten()    


    if verbose: print('Success')    
    return lc

# For pretty plotting
colors = ['dodgerblue','hotpink','darkmagenta','navy','mediumseagreen']

def test_colors():
    """Test the color scheme."""
    fig, ax = plt.subplots()
    for c in colors:
        x = np.random.uniform(-3, 3, size=20)
        y = np.random.uniform(-3, 3, size=20)
        ax.scatter(x,y,color=c)
    plt.show()

def plot_lc(lc=None,tic=None,epic=None,hostname=None,mission='tess',transit_dict=None,verbose=False):
    """Plot the light curve of a host star with predicted transit overlaid.

    Keyword arguments:
    Must have one:
        lc -- LightCurve object
        tic -- TIC ID # of the host star
        epic -- EPIC ID # of the host star
        hostname -- host star name
    mission -- which mission to get the data from (default 'tess')
    transit_dict -- predicted transits, of the form {'b':[JD1,JD2]}
    verbose -- whether to print status updates (default False)
    """
    if lc == None:
        lc = make_lc(tic=tic,epic=epic,hostname=hostname,mission=mission,verbose=verbose)

    if hostname:
        title = hostname
    else:
        if tic == None and epic == None:
            title = 'Unknown'
        elif tic:
            title = get_star_name(tic=tic)
        elif epic:
            title = get_star_name(epic=epic)

    title = title+' '+mission

    if mission.lower() == 'tess':
        offset = 2457000
    elif mission.lower() == 'k2':
        offset = 2454833

    ax = lc.scatter(s=0.1)
    
    if transit_dict:
        j = 0
        color_index = j % len(colors)
        for p in transit_dict:
            for time in transit_dict[p]:
                plt.axvline(x=time-offset, color=colors[color_index], alpha=0.5, linewidth = 3)
        j += 1
    
    plt.title(title) 
    plt.savefig('images/'+title+'.png')
    plt.show()


